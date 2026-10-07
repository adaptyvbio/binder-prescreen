"""Verdict layer: turn the two arms' numbers into flags a reviewer can act on.

The flags are ordered by how strongly they disqualify a submission, and the first one
that fires becomes the ``verdict``; every flag that fires is also listed in ``flags`` so
nothing is hidden by the short-circuit. This mirrors the evaluation order of
``novelty_scales``, where the short-circuit is what makes the levels partition cleanly.

What each flag means
--------------------
``known_target_binder``
    Whole-sequence match to a curated known binder of the target at or above
    ``tnf_known``. The submission is an existing anti-target molecule.
``target_region_match``
    The binding region (CDR3 / paratope) matches a known binder's region at or above
    ``region_match`` while the rest of the sequence may differ. This is the grafted-
    paratope case, and it is the flag the current novelty service cannot raise.
``verbatim_target_fragment``
    The submission is an exact substring of a known binder (or contains one). Catches
    copies that are too short for a significant alignment.
``existing_design``
    Whole-sequence match at or above ``known_sequence`` in the Proteinbase design
    corpus: someone already submitted or published this design. Ranked above the generic
    ``known_sequence`` because a hit among designs is the more actionable finding.
``known_sequence``
    The same whole-sequence match, against any other public reference (PDB, SwissProt,
    PLAbDab, THPdb, patents). Not target-specific: the sequence is simply not new.
``target_binder_homolog``
    Whole-sequence identity to a known binder between ``tnf_homolog`` and ``tnf_known``.
    A relative of a known binder rather than a copy.
``near_known_sequence``
    Prior-art identity between ``near_known_sequence`` and ``known_sequence``.
``batch_duplicate``
    Another submission in the same batch clusters with this one.
``no_reference_signal``
    Nothing was found on either arm *and* the sequence is short enough that the search
    may simply be blind — absence of evidence, reported as such rather than as a pass.
``pass``
    Searched, nothing found.

Context observations (reported, never flagged)
----------------------------------------------
These are the cases a whole-sequence filter gets wrong, and the reason this filter
compares binding regions at all. For antibody formats and alternative scaffolds the
framework is most of the sequence, so whole-sequence identity to a public reference is
usually framework reuse — which these competitions allow:

``known_framework_new_paratope``
    Whole-sequence match to a public reference at or above ``known_sequence``, but the
    hit's own binding region differs (below ``prior_art_region_known``). A published
    framework carrying a new paratope. Not a flag.
``design_framework_reused``
    The same, where the match is in the Proteinbase design corpus.
``shared_framework_with_target_binder``
    Whole-sequence similarity to a known target binder in the homolog band, with no
    region match. Typically "built on the same scaffold as a known anti-target agent".
``framework_similar_to_public_sequence``
    Prior-art similarity in the near-known band for a framework-dominated category.

Deliberately absent: any judgement of whether the scaffold is de novo. Framework
identity to the category reference is reported as the category's confidence, and is never
scored.
"""

from __future__ import annotations

from .config import Config

FLAG_ORDER = (
    "known_target_binder",
    "target_region_match",
    "verbatim_target_fragment",
    "existing_design",
    "known_sequence",
    "target_binder_homolog",
    "near_known_sequence",
    "batch_duplicate",
    "no_reference_signal",
    "pass",
)

# A short sequence that finds nothing may be a true negative or a blind search; this is
# the length below which we refuse to call silence a pass.
BLIND_LENGTH = 25


def evaluate(record: dict, cfg: Config) -> dict:
    """Return ``{"verdict": str, "flags": [str], "reasons": {flag: evidence}}``.

    Args:
        record: a per-sequence result carrying ``length``, ``prior_art`` (best/best_design),
            ``target`` (whole/region/substring) and ``batch`` (cluster) sections.
        cfg: thresholds.
    """
    flags: list = []
    reasons: dict = {}

    pa = record.get("prior_art") or {}
    best = pa.get("best")
    best_design = pa.get("best_design")
    tgt = record.get("target") or {}
    whole = tgt.get("whole")
    region = tgt.get("region") or {}
    substring = tgt.get("substring")
    batch = record.get("batch") or {}
    length = record.get("length") or 0
    region_kind = record.get("region_kind")

    context: dict = {}
    pa_region = pa.get("region") or {}
    # None means "no region could be compared", which is not evidence of anything. Keep it
    # distinct from 0.0 ("compared, and different") all the way to the suppression rule.
    _pa_region_raw = pa_region.get("region_identity")
    pa_region_verified = _pa_region_raw is not None
    pa_region_id = float(_pa_region_raw or 0.0)
    # For antibody formats and alternative scaffolds the framework dominates the
    # sequence, so a whole-sequence match to a public reference (or to a known binder)
    # is usually framework sharing, which these competitions allow. Those categories are
    # therefore judged on the binding region, and whole-sequence similarity alone is
    # recorded as context instead of raising a flag.
    region_aware = region_kind != "whole"

    whole_fident = float(whole.get("fident") or 0.0) if whole else 0.0
    whole_qcov = float(whole.get("qcov") or 0.0) if whole else 0.0
    whole_sim = float(whole["similarity_check"]) if whole else 0.0
    region_id = float(region.get("region_identity") or 0.0)
    region_exact = bool(region.get("region_exact"))
    region_len = len(region.get("region_query_seq") or "")
    pa_sim = float(best["similarity_check"]) if best else 0.0
    design_sim = float(best_design["similarity_check"]) if best_design else 0.0

    region_matches = region_exact or region_id >= cfg.region_match
    # Whole-sequence "already a known binder". For antibody and scaffold formats the
    # paratope is the binding determinant, so a known-binder call requires BOTH a high
    # whole-sequence match (identity and coverage past the measured cut) AND the binding
    # region to match. A known binder resubmitted verbatim clears both; a known framework
    # carrying a NEW paratope clears the whole-sequence cut but not the region, and must
    # pass — reusing a published scaffold is allowed. For non-antibody formats
    # (miniprotein, peptide) the whole sequence IS the binder, so the composite suffices.
    if region_aware:
        known_binder = (
            whole_fident >= cfg.tnf_known_identity
            and whole_qcov >= cfg.tnf_known_coverage
            and region_matches
        )
    else:
        known_binder = whole_sim >= cfg.tnf_known
    if whole and known_binder:
        flags.append("known_target_binder")
        reasons["known_target_binder"] = {
            "reference": whole["target"],
            "identity": round(whole_fident, 4),
            "query_coverage": round(whole_qcov, 4),
            "similarity_check": round(whole_sim, 4),
            "binder_status": whole.get("binder_status"),
            "named_agent": whole.get("named_agent"),
        }
    # Known anti-target paratope on any framework: exact CDR3 match, or CDR3 identity at
    # or above region_match on a matching chain. Only for region-bearing categories; for a
    # whole-sequence category the region restates the whole-sequence number.
    if (
        region_aware
        and region_len >= cfg.region_min_length
        and (region_exact or region_id >= cfg.region_match)
    ):
        flags.append("target_region_match")
        reasons["target_region_match"] = {
            "reference": region.get("region_reference"),
            "region_identity": round(region_id, 4),
            "exact": region_exact,
            "region": region.get("region_query_seq"),
            "reference_region": region.get("region_reference_seq"),
            "chain": region.get("region_chain"),
            "binder_status": region.get("binder_status"),
            "named_agent": region.get("named_agent"),
            "citation": region.get("citation"),
            "supporting_exact_regions": region.get("supporting_exact_regions"),
        }
    # Verbatim-copy check. For antibody and scaffold formats the exact-CDR / paratope arm
    # already catches a copied binding region, and a whole-sequence substring there only
    # matches shared framework (a short framework-derived peptide in the corpus is a
    # substring of any same-class framework), so the flag is reserved for whole-sequence
    # categories — peptides and miniproteins, where a verbatim copy has no region arm to
    # catch it. The match is still recorded as context for region-aware categories.
    if substring:
        if region_aware:
            context["framework_contains_reference_fragment"] = substring
        else:
            flags.append("verbatim_target_fragment")
            reasons["verbatim_target_fragment"] = substring

    # Suppressing a public-sequence match because the paratope looks new is only
    # defensible when the paratope was actually examined. When the region could not be
    # computed there is no finding to suppress with, so the flag stands and the evidence
    # records that it is unverified — a verbatim copy of a published binder whose CDRs
    # happened not to extract must not be reported as a pass.
    paratope_known = (
        (not region_aware)
        or not pa_region_verified
        or pa_region_id >= cfg.prior_art_region_known
    )
    if design_sim >= cfg.known_sequence:
        evidence = {
            "reference": best_design["target"],
            "db": best_design["db"],
            "similarity_check": round(design_sim, 4),
            "region_identity": round(pa_region_id, 4) if pa_region_verified else None,
            "paratope_verified": pa_region_verified if region_aware else None,
        }
        if paratope_known:
            flags.append("existing_design")
            reasons["existing_design"] = evidence
        else:
            context["design_framework_reused"] = evidence
    if pa_sim >= cfg.known_sequence:
        evidence = {
            "reference": best["target"],
            "db": best["db"],
            "identity": round(best["fident"], 4),
            "query_coverage": round(best["qcov"], 4),
            "similarity_check": round(pa_sim, 4),
            "region_identity": round(pa_region_id, 4) if pa_region_verified else None,
            "paratope_verified": pa_region_verified if region_aware else None,
            "region_hit": pa_region.get("hit"),
        }
        if paratope_known:
            flags.append("known_sequence")
            reasons["known_sequence"] = evidence
        else:
            # The framework is public but the paratope is not: exactly the case these
            # competitions permit, so it is context and never a flag.
            context["known_framework_new_paratope"] = evidence
    if cfg.tnf_homolog <= whole_sim < cfg.tnf_known:
        evidence = {
            "reference": whole["target"],
            "similarity_check": round(whole_sim, 4),
            "region_identity": round(region_id, 4) if region_aware else None,
        }
        if not region_aware:
            flags.append("target_binder_homolog")
            reasons["target_binder_homolog"] = evidence
        elif region_id < cfg.region_match:
            # Framework-dominated category, similar overall to a known binder but with a
            # different binding region: scaffold reuse, which is allowed. When the region
            # *does* match, target_region_match above already carries the finding.
            context["shared_framework_with_target_binder"] = evidence
    if cfg.near_known_sequence <= pa_sim < cfg.known_sequence:
        evidence = {
            "reference": best["target"],
            "db": best["db"],
            "similarity_check": round(pa_sim, 4),
        }
        if region_aware:
            context["framework_similar_to_public_sequence"] = evidence
        else:
            flags.append("near_known_sequence")
            reasons["near_known_sequence"] = evidence
    if (batch.get("cluster_size") or 1) > 1:
        flags.append("batch_duplicate")
        reasons["batch_duplicate"] = {
            "cluster_representative": batch.get("cluster_representative"),
            "cluster_size": batch.get("cluster_size"),
        }
    if not flags:
        if pa_sim == 0.0 and whole_sim == 0.0 and length <= BLIND_LENGTH:
            flags.append("no_reference_signal")
            reasons["no_reference_signal"] = {
                "length": length,
                "note": "no hit on either arm at a length where k-mer search is unreliable",
            }
        else:
            flags.append("pass")

    verdict = next(f for f in FLAG_ORDER if f in flags)
    return {
        "verdict": verdict,
        "flags": flags,
        "reasons": reasons,
        "context": context,
    }
