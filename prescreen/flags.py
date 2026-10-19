"""Verdict layer: turn the two arms' numbers into flags a reviewer can act on.

The flags are ordered by how strongly they disqualify a submission, and the first one
that fires becomes the ``verdict``; every flag that fires is also listed in ``flags`` so
nothing is hidden by the short-circuit. This mirrors the evaluation order of
``novelty_scales``, where the short-circuit is what makes the levels partition cleanly.

What each flag means
--------------------
``known_target_binder``
    Whole-sequence match to a curated known binder of the target at or above
    ``target_known``. The submission is an existing anti-target molecule.
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
    Whole-sequence identity to a known binder between ``target_homolog`` and ``target_known``.
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

#: Whole-sequence identity AND coverage at which a hit is the reference itself. Used only
#: to keep a verbatim copy flagged when its region match is exact but germline-common.
VERBATIM_IDENTITY = 0.99


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
    if pa.get("carrier_explained"):
        context["known_carrier_prior_art"] = {
            "hits": pa["carrier_explained"],
            "note": "the wild-type Fc this submission is fused to is public; its hits "
            "are carrier, not the binder, and were left out of the prior-art call",
        }
    # For antibody formats and alternative scaffolds the framework dominates the
    # sequence, so a whole-sequence match to a public reference (or to a known binder)
    # is usually framework sharing, which these competitions allow. Those categories are
    # therefore judged on the binding region, and whole-sequence similarity alone is
    # recorded as context instead of raising a flag.
    region_aware = region_kind != "whole"

    whole_fident = float(whole.get("fident") or 0.0) if whole else 0.0
    whole_qcov = float(whole.get("qcov") or 0.0) if whole else 0.0
    whole_sim = float(whole["similarity_check"]) if whole else 0.0
    # None means the target region arm did not run at all (no extractable focus region, or
    # no reference of a comparable category). It cannot corroborate a known-binder call,
    # and it must not be read as a clean region either.
    region_compared = region.get("region_compared", True) and region.get(
        "region_identity"
    ) is not None
    region_id = float(region.get("region_identity") or 0.0)
    region_exact = bool(region.get("region_exact"))
    region_len = len(region.get("region_query_seq") or "")
    pa_sim = float(best["similarity_check"]) if best else 0.0
    design_sim = float(best_design["similarity_check"]) if best_design else 0.0

    # False when the region matched is germline-common - carried by more than
    # ``region_corpus_max`` sequences in a general antibody corpus - so the match says the
    # query shares a framework, not a paratope. Still reported, just not load-bearing.
    # See :func:`prescreen._best_region`.
    region_decisive = bool(region.get("region_decisive", True))
    region_matches = region_decisive and (region_exact or region_id >= cfg.region_match)
    if region_aware and region_compared and not region_decisive and (
        region_exact or region_id >= cfg.region_match
    ):
        context["target_region_match_not_decisive"] = {
            "reference": region.get("region_reference"),
            "region_identity": round(region_id, 4),
            "chain": region.get("region_chain"),
            "corpus_occurrences": region.get("corpus_occurrences"),
            "note": "matched region is germline-common in a general antibody corpus, so "
            "it indicates a shared framework rather than a known paratope; reported as "
            "context only",
        }
    # An antibody-shaped sequence the numbering refused. The paratope arm did not run,
    # so "no region match" here means "not looked at", and a reviewer who reads this row
    # as a clean paratope is reading a check that never happened. Context, not a flag: it
    # changes no verdict, it says what the verdict does not cover.
    numbering_reason = (record.get("category_evidence") or {}).get(
        "antibody_numbering_reason"
    )
    if numbering_reason:
        context["antibody_numbering_failed"] = {
            "reason": numbering_reason,
            "category": record.get("category"),
            "note": "sequence was not recognised as an antibody domain, so the paratope "
            "arm did not run; it was judged on whole-sequence evidence only",
        }
    if region_aware and not region_compared:
        # Surfaced as context rather than a flag: it changes no verdict, but a reviewer
        # must not read this record's silence on the binding region as a negative result.
        context["target_region_not_compared"] = {
            "region_kind": region_kind,
            "note": "no comparable reference region; whole-sequence evidence only",
        }
    # Whole-sequence "already a known binder". For antibody and scaffold formats the
    # paratope is the binding determinant, so a known-binder call requires BOTH a high
    # whole-sequence match (identity and coverage past the measured cut) AND the binding
    # region to match. A known binder resubmitted verbatim clears both; a known framework
    # carrying a NEW paratope clears the whole-sequence cut but not the region, and must
    # pass — reusing a published scaffold is allowed. For non-antibody formats
    # (miniprotein, peptide) the whole sequence IS the binder, so the composite suffices.
    verbatim = False
    if region_aware:
        # Requiring the region to corroborate is right when the region arm had a fair
        # chance to. It did not when the matched reference's class contributes no CDRs to
        # the paratope index, so a byte-identical resubmission of such a binder scores low
        # on the region arm and could not be called a known binder. That is the arm having
        # nothing to say, not evidence of a new paratope.
        #
        # The waiver is scoped to antibody formats because they are the only ones the CDR
        # index serves. An affibody, monobody or DARPin is compared by projected paratope
        # against same-class references (the scaffold branch of
        # ``TargetReference.compare_region``), which never consults that index — its
        # classes are absent from it by construction, not by a gap in curation. Waiving
        # the paratope requirement for them would flag a legitimately reused scaffold
        # carrying a new binding surface on whole-sequence similarity alone, which is the
        # one thing this rule exists to prevent.
        #
        # Deliberately NOT an identity cut: a framework-reuse chimera swaps one CDR3, which
        # in a 432-aa Fab still leaves 97% identity, so any fixed identity threshold high
        # enough to mean "verbatim" for a nanobody flags legitimate reuse in a Fab.
        uncorroborable = (
            region_kind == "antibody_cdrs"
            and whole is not None
            and not whole.get("region_reference_available", True)
        )
        # A germline-common region cannot tell a known binder from a framework that
        # shares it - except when the submission IS the reference: the whole sequence is
        # identical and so is the region. Requiring the region to match exactly keeps
        # the framework-reuse case above out of reach - a swapped CDR3 is never exact.
        verbatim = (
            region_compared
            and region_exact
            and whole_fident >= VERBATIM_IDENTITY
            and whole_qcov >= VERBATIM_IDENTITY
        )
        known_binder = (
            whole_fident >= cfg.target_known_identity
            and whole_qcov >= cfg.target_known_coverage
            and (region_matches or uncorroborable or verbatim)
        )
    else:
        known_binder = whole_sim >= cfg.target_known
    # A bare-Fc query whose hit shares no engineered CH3 patch with the reference - Fc is
    # ~95% conserved, so any Fc-bearing reference matches it: the similarity is the
    # carrier, not the binder. Reported as context, like every other framework-reuse
    # observation.
    if whole is not None and whole.get("carrier_explained"):
        known_binder = False
        context["known_carrier_new_paratope"] = {
            "reference": whole["target"],
            "identity": round(whole_fident, 4),
            "query_coverage": round(whole_qcov, 4),
            "carrier_patch_query": whole.get("carrier_patch_query"),
            "carrier_patch_reference": whole.get("carrier_patch_reference"),
            "carrier_patch_identity": whole.get("carrier_patch_identity"),
            "note": "matched an Fc-bearing reference on carrier sequence; the "
            "submission and the reference share no engineered CH3 patch",
        }
    if whole and known_binder:
        flags.append("known_target_binder")
        reasons["known_target_binder"] = {
            "reference": whole["target"],
            "identity": round(whole_fident, 4),
            "query_coverage": round(whole_qcov, 4),
            "similarity_check": round(whole_sim, 4),
            "binder_status": whole.get("binder_status"),
            "named_agent": whole.get("named_agent"),
            # False when the flag rests on the whole-sequence match alone because the
            # matched reference's class contributes no CDRs to the index. The finding is
            # real but uncorroborated, and a reviewer should see which it is.
            "paratope_corroborated": bool(region_matches),
            # True when the call rests on the submission being the reference verbatim,
            # its exact region match being germline-common and so not decisive alone.
            "verbatim_copy": bool(verbatim and not region_matches),
            "reference_class": whole.get("reference_class"),
        }
    # Known anti-target paratope on any framework: exact CDR3 match, or CDR3 identity at
    # or above region_match on a matching chain. Only for region-bearing categories; for a
    # whole-sequence category the region restates the whole-sequence number.
    if (
        region_aware
        and region_decisive
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
        # The demotion below is right for a short framework-derived scrap, which is a
        # substring of any same-class framework and would fire on every antibody. It is
        # wrong for a verbatim copy of an attributed binder: a known GLP-1 agonist fused
        # to a carrier VHH types as a nanobody, so the region arm compares the VHH's CDRs,
        # finds nothing, and vetoes the known-binder call - while the agonist sits inside
        # the construct, already found, recorded and ignored. Measured: exendin-4 joined
        # to a VHH by a (G4S)3 linker passed both arms.
        promote = (
            substring.get("attributed")
            and (substring.get("matched_length") or 0) >= cfg.verbatim_fragment_min_length
        )
        if region_aware and not promote:
            context["framework_contains_reference_fragment"] = substring
        else:
            flags.append("verbatim_target_fragment")
            reasons["verbatim_target_fragment"] = substring

    # Suppressing a public-sequence match because the paratope looks new is only
    # defensible when the paratope was actually examined. When the region could not be
    # computed there is no finding to suppress with, so the flag stands and the evidence
    # records that it is unverified — a verbatim copy of a published binder whose CDRs
    # happened not to extract must not be reported as a pass.
    # ...and only on the paratope of the HIT BEING JUDGED. ``best`` and ``best_design``
    # are frequently different molecules, and the design hit is often not the one that
    # supplied the strongest region match, so a single best-over-hits number would let an
    # unrelated arm's paratope decide this hit's verdict in both directions — promoting a
    # design whose own paratope differs, or excusing one whose paratope was never looked
    # at. ``per_hit`` holds one identity per hit examined, so each flag reads its own.
    #
    # Note 0.0 is a real answer ("compared, and different") and must suppress; only a hit
    # ABSENT from ``per_hit`` was never examined, and there the flag stands with the
    # evidence recording the paratope as unverified.
    per_hit = pa_region.get("per_hit") or {}

    def paratope_known(hit: dict | None) -> tuple[bool, float | None]:
        """``(suppressible, identity_for_this_hit)`` for one prior-art hit."""
        if not region_aware:
            return True, None
        value = per_hit.get((hit or {}).get("target"))
        if value is None:
            return True, None
        return value >= cfg.prior_art_region_known, value

    if design_sim >= cfg.known_sequence:
        known, hit_region_id = paratope_known(best_design)
        evidence = {
            "reference": best_design["target"],
            "db": best_design["db"],
            "similarity_check": round(design_sim, 4),
            "region_identity": None if hit_region_id is None else round(hit_region_id, 4),
            "paratope_verified": (hit_region_id is not None) if region_aware else None,
        }
        if known:
            flags.append("existing_design")
            reasons["existing_design"] = evidence
        else:
            context["design_framework_reused"] = evidence
    if pa_sim >= cfg.known_sequence:
        known, hit_region_id = paratope_known(best)
        evidence = {
            "reference": best["target"],
            "db": best["db"],
            "identity": round(best["fident"], 4),
            "query_coverage": round(best["qcov"], 4),
            "similarity_check": round(pa_sim, 4),
            "region_identity": None if hit_region_id is None else round(hit_region_id, 4),
            "paratope_verified": (hit_region_id is not None) if region_aware else None,
            "region_hit": pa_region.get("hit"),
        }
        if known:
            flags.append("known_sequence")
            reasons["known_sequence"] = evidence
        else:
            # The framework is public but the paratope is not: exactly the case these
            # competitions permit, so it is context and never a flag.
            context["known_framework_new_paratope"] = evidence
    if whole and whole_sim >= cfg.target_homolog:
        evidence = {
            "reference": whole["target"],
            "similarity_check": round(whole_sim, 4),
            "region_identity": round(region_id, 4) if region_aware else None,
        }
        if not region_aware:
            # A band, deliberately: above target_known the whole-sequence rule already fired.
            if whole_sim < cfg.target_known:
                flags.append("target_binder_homolog")
                reasons["target_binder_homolog"] = evidence
        elif region_id < cfg.region_match and "known_target_binder" not in flags:
            # Framework-dominated category, similar overall to a known binder but with a
            # different binding region: scaffold reuse, which is allowed. When the region
            # *does* match, target_region_match above already carries the finding.
            #
            # No upper bound here. Bounding it at target_known assumed anything above that cut
            # would be caught by known_target_binder, but that flag also needs the region to
            # match — so a submission 90%+ identical to a known binder with an unmatched
            # region fell out of both rules and the target arm said nothing at all. The
            # report went quiet exactly where it should be loudest.
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
