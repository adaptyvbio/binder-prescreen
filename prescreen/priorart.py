"""The prior-art arm: is this sequence already public, and where.

The search itself is delegated to :mod:`prescreen.refdb`, the reference-database layer,
which resolves the eight arms (``pdb``, ``swissprot``, ``plabdab``, ``plabdab_nano``,
``therasabdab``, ``thpdb``, ``proteinbase_public``, ``uspto``), routes every query by
length into a tuned parameter set, and raises before searching if an arm is absent —
because a missing database reads exactly like a sequence with no prior art, which is the
one answer this filter must never give by accident.

This module adds the two things the prescreen needs on top of a hit list:

1. **A separate design channel.** A hit among published designs
   (``proteinbase_public``) is a different finding from a hit in SwissProt or a patent,
   so it is reported on its own.
2. **A paratope-aware check.** Each hit comes back with its own sequence (``tseq``), so
   the hit can be re-numbered and its binding region compared with the query's. Without
   this, every designed nanobody built on a published framework looks like prior art.

Scoring follows the deployed convention so numbers stay comparable with ``db_check``:
``similarity_check = fident * qcov`` and ``distance = qlen * (1 - similarity_check)``.
Threshold on that composite, never on ``fident`` alone: for a short query the aligner
often trims mutated termini and reports ``fident = 1.0`` with ``qcov < 1``.
"""

from __future__ import annotations

from pathlib import Path

from . import mmseqs, refdb
from .chains import _relabel
from .config import Config


def arms_for(cfg: Config) -> list:
    """Arms to search: every public arm, the patent arm included.

    The patent arm is 10.2M sequences and costs seconds per query, but a prior-art screen
    that cannot see patents is not a prior-art screen: a binder claimed in a granted patent
    and never deposited anywhere else is invisible without it. The cost is amortised over
    the batch, which is why submissions are screened in one pass rather than one at a time.
    ``Config.include_patent_arm = False`` (``--no-patent``) drops it.
    """
    return [a for a in refdb.PUBLIC_ARMS if a != "uspto" or cfg.include_patent_arm]


def search(
    records: dict,
    cfg: Config,
    workdir: str | Path,
    db_names: list | None = None,
) -> dict:
    """Search a batch of sequences against the prior-art and design databases.

    Args:
        records: ``{query_id: sequence}``.
        cfg: prescreen configuration (database root, thresholds, mmseqs binary).
        workdir: scratch directory for FASTA and tabular output.
        db_names: override the arm list; defaults to :func:`arms_for`.

    Returns:
        ``{query_id: {"best": row|None, "best_design": row|None, "per_db": {arm: row},
        "top": [row]}}`` where each row is an MMseqs2 hit augmented with ``db`` (the arm
        name), ``profile`` (which length-routed parameter set found it),
        ``similarity_check`` and ``distance``. ``top`` holds the best
        ``cfg.prior_art_top_hits`` hits overall, kept for the region-level check.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    arms = list(db_names or arms_for(cfg))
    roots = cfg.db_root.split(":") if cfg.db_root else None

    present = refdb.available(roots)
    searched = [a for a in arms if present.get(a)]
    missing = [a for a in arms if not present.get(a)]
    # refdb.search_fasta raises for an arm it cannot resolve, but it is only ever handed
    # the arms that already resolved, so the invariant this module's docstring states has
    # to be enforced here or not at all. Silently searching seven arms when eight were
    # asked for reports "no prior art" on evidence that was never gathered.
    if missing and not cfg.allow_missing_arms:
        raise refdb.ReferenceDbMissing(
            f"reference database(s) not found for arm(s): {', '.join(missing)}. "
            f"Build them with scripts/build_reference_dbs.sh, point PRESCREEN_DB_ROOT at "
            f"them, or pass --allow-missing-arms to screen against the rest and accept a "
            f"weaker 'no prior art' result."
        )

    results = {
        qid: {"best": None, "best_design": None, "per_db": {}, "top": []}
        for qid in records
    }
    if not searched:
        for slot in results.values():
            slot["arms_searched"] = []
            slot["arms_missing"] = missing
        return results

    query_fasta = mmseqs.write_fasta(records, workdir / "queries.fasta")
    hits = refdb.search_fasta(
        query_fasta,
        arms=searched,
        roots=roots,
        out_dir=workdir / "m8",
        threads=cfg.threads or 8,
        profile="auto",
    )
    design = {
        name for name, arm in refdb.ARMS.items() if arm.arm == "designs"
    }
    for row in hits.to_dict("records"):
        qid = row.get("query")
        if qid not in results:
            continue
        row["db"] = row.get("arm")
        row["similarity_check"] = float(row.get("similarity") or 0.0)
        row["distance"] = round(
            float(row.get("qlen") or 0) * (1 - row["similarity_check"]), 2
        )
        slot = results[qid]
        slot["per_db"][row["db"]] = row
        key = "best_design" if row["db"] in design else "best"
        prev = slot[key]
        if prev is None or row["similarity_check"] > prev["similarity_check"]:
            slot[key] = row
        slot["top"].append(row)

    keep = max(1, cfg.prior_art_top_hits)
    for slot in results.values():
        slot["top"] = sorted(
            slot["top"], key=lambda r: r["similarity_check"], reverse=True
        )[:keep]
        slot["arms_searched"] = searched
        slot["arms_missing"] = missing
    return results


def merge_chains(per_chain: dict, index, cfg: Config) -> dict:
    """Collapse per-chain prior-art slots into one slot per submission.

    ``search`` is chain-agnostic — it takes opaque query ids and returns the same ids — so
    the chain bookkeeping lives here. Every row is relabelled with its submission id and
    the chain that found it before being merged, so the synthetic query ids never reach the
    report.
    """
    keep = max(1, cfg.prior_art_top_hits)
    merged: dict = {}
    for sid, qids in index.by_submission.items():
        slot = {
            "best": None,
            "best_design": None,
            "per_db": {},
            "top": [],
            "per_chain": {},
            "arms_searched": [],
            "arms_missing": [],
        }
        for qid in qids:
            chain, domain = index.chain_of(qid), index.domain_of(qid)
            part = per_chain.get(qid)
            if part is None:
                continue
            # All chains of a batch share one refdb.available() call, so these agree.
            slot["arms_searched"] = part.get("arms_searched", [])
            slot["arms_missing"] = part.get("arms_missing", [])
            for row in part.get("top") or []:
                _relabel(row, sid, chain, domain)
            for key in ("best", "best_design"):
                row = part.get(key)
                if row is not None:
                    _relabel(row, sid, chain, domain)
                # Strict ``>``: walking chains in order, a tie goes to the lowest chain
                # number, which the template's convention makes the heavy chain.
                if row is not None and (
                    slot[key] is None
                    or row["similarity_check"] > slot[key]["similarity_check"]
                ):
                    slot[key] = row
            for arm, row in (part.get("per_db") or {}).items():
                _relabel(row, sid, chain, domain)
                prev = slot["per_db"].get(arm)
                if prev is None or row["similarity_check"] > prev["similarity_check"]:
                    slot["per_db"][arm] = row
            slot["top"].extend(part.get("top") or [])
            # String key: report.json is plain json.dumps, which cannot take a tuple key.
            slot["per_chain"][f"{chain}.{domain}"] = {
                "best": part.get("best"),
                "best_design": part.get("best_design"),
                "n_hits": len(part.get("top") or []),
            }
        # Sort on the composite ALONE. Python's sort is stable and each chain's list
        # arrived already sorted by it, so a single-chain merge returns the identical
        # order; adding db or chain to the key would reorder ties within one chain.
        slot["top"] = sorted(
            slot["top"], key=lambda r: r["similarity_check"], reverse=True
        )[:keep]
        slot["best_chain"] = (slot["best"] or {}).get("chain")
        slot["best_design_chain"] = (slot["best_design"] or {}).get("chain")
        merged[sid] = slot
    return merged


def best_region_vs_hits(
    chain_regions: list,
    category: str,
    region_kind: str,
    hits: list,
    cache: dict | None = None,
) -> dict:
    """Best :func:`region_identity_vs_hits` result across a submission's input chains.

    Same reason as the target arm: a multivalent construct numbered as one string has its
    second domain mislabelled, so each submitted chain is compared in its own right. A
    single-chain submission does exactly one comparison, as before.
    """
    best = None
    # ``per_hit`` is merged across chains rather than taken from the winning one: a
    # submission's heavy chain and light chain are compared against the same hit list, and
    # the question each flag asks is "what is this hit's paratope identity", not "which
    # chain won". The strongest comparison for a given hit, from whichever chain produced
    # it, is the answer.
    merged: dict = {}
    for regions in chain_regions or []:
        got = region_identity_vs_hits(regions, category, region_kind, hits, cache=cache)
        for target, value in (got.get("per_hit") or {}).items():
            prev = merged.get(target)
            if prev is None or value > prev:
                merged[target] = value
        if best is None or (got.get("region_identity") or -1.0) > (
            best.get("region_identity") or -1.0
        ):
            best = got
    if best is None:
        best = region_identity_vs_hits({}, category, region_kind, hits, cache=cache)
    best = dict(best)
    best["per_hit"] = merged
    return best


def region_identity_vs_hits(
    query_regions: dict,
    category: str,
    region_kind: str,
    hits: list,
    cache: dict | None = None,
) -> dict:
    """Compare the query's binding region with the binding regions of its prior-art hits.

    This is the check that separates "my scaffold is a known framework" from "my binder
    is a known binder". A designed nanobody built on a published framework matches
    PLAbDab at high whole-sequence identity by construction; only if the *paratope* also
    matches is the molecule itself prior art.

    Args:
        query_regions: regions of the query, from :func:`prescreen.regions.extract`.
        category: query category slug.
        region_kind: query region kind; ``whole`` returns an empty result because the
            region and the sequence are the same thing.
        hits: prior-art hit rows carrying ``tseq`` (the hit's own sequence).
        cache: optional ``{hit_sequence: regions}`` map reused across queries, which
            matters because popular framework hits recur throughout a batch.

    Returns:
        ``{"region_identity": float|None, "hit": str|None, "db": str|None,
        "hit_region": str|None, "per_hit": {target_id: float}}``.

        ``per_hit`` carries one identity per hit examined, which is what lets a flag be
        judged on the paratope of the hit it is about. A single best-over-hits number
        cannot do that: ``best`` and ``best_design`` are usually different molecules, and
        the hit supplying the maximum region identity is rarely the one with the highest
        whole-sequence similarity, so one shared number either suppresses a flag using
        evidence about something else or fails to suppress at all. A hit absent from
        ``per_hit`` was never compared — its paratope is unknown, not new.

        ``region_identity`` is ``None`` when no region could be compared at all — the query's focus region would not extract, or no hit carried a
        usable one. That is NOT the same as 0.0, which means a region WAS compared and
        differs, and the two must not be conflated: a caller that reads "could not compute"
        as "the paratope is new" turns a failed check into a clean bill of health for a
        verbatim copy of a published binder.
    """
    from .regions import FOCUS, extract, focus_region, identity_fn_for

    empty = {
        "region_identity": None,
        "hit": None,
        "db": None,
        "hit_region": None,
        "chain": None,
        "per_hit": {},
    }
    if region_kind == "whole":
        return dict(empty)
    _chain, q_focus = focus_region(query_regions)
    if not q_focus:
        return dict(empty)
    identity = identity_fn_for(category)
    cache = cache if cache is not None else {}
    best = dict(empty)
    per_hit: dict = {}
    compared = False
    for hit in hits:
        tseq = (hit.get("tseq") or "").strip().upper()
        if not tseq or len(tseq) < len(q_focus):
            continue
        if tseq in cache:
            hit_foci = cache[tseq]
        else:
            # The hit misassigns its own chains exactly as a query would: a bivalent
            # reference numbered as one string files its second domain under L. Compare
            # against each of the hit's variable domains, so a multivalent reference is
            # not reduced to whichever domain antpack happened to label H.
            from .chains import variable_domains

            units = variable_domains(tseq) or [tseq]
            hit_foci = []
            for unit in units:
                _hc, h_focus = focus_region(extract(unit, category, region_kind))
                if h_focus:
                    hit_foci.append(h_focus)
            cache[tseq] = hit_foci
        if not hit_foci:
            continue
        h_focus = max(hit_foci, key=lambda f: float(identity(q_focus, f)))
        value = float(identity(q_focus, h_focus))
        compared = True
        target = hit.get("target")
        if target is not None and value > (per_hit.get(target, -1.0)):
            per_hit[target] = value
        if value > (best["region_identity"] or 0.0):
            best = {
                "region_identity": round(value, 4),
                "hit": hit.get("target"),
                "db": hit.get("db"),
                "hit_region": h_focus,
                # Which input chain's hit supplied this region, so a reviewer can see
                # whether a light-chain hit is what demoted a heavy-chain match.
                "chain": hit.get("chain"),
            }
    if compared and best["region_identity"] is None:
        # Every comparison scored 0.0, so the loop never beat the initial None. The region
        # WAS compared, and the answer is zero.
        best["region_identity"] = 0.0
    best["per_hit"] = per_hit
    return best
