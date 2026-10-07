"""The prior-art arm: is this sequence already public, and where.

The search itself is delegated to :mod:`prescreen.refdb`, the reference-database layer,
which resolves the nine arms (``pdb``, ``swissprot``, ``plabdab``, ``plabdab_nano``,
``therasabdab``, ``thpdb``, ``proteinbase_public``, ``proteinbase_internal``, ``uspto``),
routes every query by length into a tuned parameter set, and raises before searching if
an arm is absent — because a missing database reads exactly like a sequence with no prior
art, which is the one answer this filter must never give by accident.

This module adds the three things the prescreen needs on top of a hit list:

1. **Visibility discipline.** ``proteinbase_internal`` holds other entrants' unpublished
   submissions. It may raise an internal duplicate flag for the organisers and must never
   appear as evidence shown to a competitor; it is excluded unless
   ``Config.organiser_mode`` is set.
2. **A separate design channel.** A hit among published designs
   (``proteinbase_public``) is a different finding from a hit in SwissProt or a patent,
   so it is reported on its own.
3. **A paratope-aware check.** Each hit comes back with its own sequence (``tseq``), so
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
from .config import Config


def arms_for(cfg: Config) -> list:
    """Arms to search: public arms, minus the patent arm unless asked for.

    The patent arm is 10.2M sequences and costs seconds per query, so it belongs in an
    asynchronous batch rather than the submission path. ``proteinbase_internal`` is added
    only in organiser mode.
    """
    arms = [a for a in refdb.PUBLIC_ARMS if a != "uspto" or cfg.include_patent_arm]
    if cfg.organiser_mode and "proteinbase_internal" in refdb.ARMS:
        arms.append("proteinbase_internal")
    return arms


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
        ``{"region_identity": float, "hit": str|None, "db": str|None,
        "hit_region": str|None}``; identity 0.0 when nothing comparable was found.
    """
    from .regions import FOCUS, extract, focus_region, identity_fn_for

    empty = {"region_identity": 0.0, "hit": None, "db": None, "hit_region": None}
    if region_kind == "whole":
        return empty
    _chain, q_focus = focus_region(query_regions)
    if not q_focus:
        return empty
    identity = identity_fn_for(category)
    cache = cache if cache is not None else {}
    best = dict(empty)
    for hit in hits:
        tseq = (hit.get("tseq") or "").strip().upper()
        if not tseq or len(tseq) < len(q_focus):
            continue
        if tseq in cache:
            hit_regions = cache[tseq]
        else:
            hit_regions = extract(tseq, category, region_kind)
            cache[tseq] = hit_regions
        _hc, h_focus = focus_region(hit_regions)
        if not h_focus:
            continue
        value = float(identity(q_focus, h_focus))
        if value > best["region_identity"]:
            best = {
                "region_identity": round(value, 4),
                "hit": hit.get("target"),
                "db": hit.get("db"),
                "hit_region": h_focus,
            }
    return best
