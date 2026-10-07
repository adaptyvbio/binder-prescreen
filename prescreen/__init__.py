"""Sequence-level prescreen for protein-design competition submissions.

The filter answers two questions per submitted sequence, both at the sequence level and
both fast enough to run at submission time:

1. **Is this sequence already known?** Batched MMseqs2 search against the public corpora
   (PDB, SwissProt, PLAbDab, THPdb, USPTO patents) and against the Proteinbase design
   corpus, with length-aware settings so short peptides are not silently missed.
2. **Is it similar to a known binder of the target?** Whole-sequence and *binding-region*
   comparison against a curated set of known binders for the target — the region arm
   being what catches a known paratope grafted onto a new framework.

It deliberately does **not** ask whether the scaffold is de novo. Reusing a published
framework is allowed in these competitions, so the category (nanobody, DARPin, affibody,
monobody, miniprotein, peptide, ...) is used only to decide which binding region to
compare, and framework identity is reported as context rather than scored.

Typical use::

    from prescreen import Config, screen
    results = screen({"sub_001": "QVQLVESGG..."}, Config.from_env(),
                     target_fasta="data/tnfa_binders.fasta")

or from the shell::

    prescreen submissions.fasta --target data/tnfa_binders.fasta -o report
"""

from __future__ import annotations

import os
import tempfile
import time
from pathlib import Path

from . import cluster as _cluster
from . import flags as _flags
from . import mmseqs, priorart, refdb
from .chains import merge_whole, split_submissions
from .classify import classify, compare_declared
from .config import Config, resolve_db_root
from .regions import extract, focus_region
from .target import TargetReference

__all__ = [
    "Config",
    "TargetReference",
    "screen",
    "screen_fasta",
    "classify",
    "compare_declared",
    "resolve_db_root",
    "__version__",
]

__version__ = "0.1.0"


def _best_region(reference, ann: dict) -> dict:
    """Best target-region comparison across a submission's input chains.

    Each chain is compared in its own right and the strongest result wins — an exact match
    beats any inexact one. ``supporting_exact_regions`` is summed across chains, so a Fab
    scored chain-by-chain corroborates exactly as it does scored whole.
    """
    best, support = None, 0
    for regions in ann["chain_regions"]:
        hit = reference.compare_region(
            regions, ann["category"], ann["region_kind"]
        ).as_dict()
        support += hit.get("supporting_exact_regions") or 0
        if best is None:
            best = hit
            continue
        rank = (hit.get("region_exact"), hit.get("region_identity") or -1.0)
        prev = (best.get("region_exact"), best.get("region_identity") or -1.0)
        if rank > prev:
            best = hit
    if best is None:
        return {}
    best["supporting_exact_regions"] = support
    return best


def _annotate(
    records: dict,
    declared: dict | None = None,
    chains: dict | None = None,
    chain_seqs: dict | None = None,
) -> dict:
    """Classify and region-annotate a batch of submissions."""
    declared = declared or {}
    chains = chains or {}
    chain_seqs = chain_seqs or {}
    out: dict = {}
    for qid, seq in records.items():
        cl = classify(seq)
        regions = extract(seq, cl.category, cl.region_kind)
        chain, focus = focus_region(regions)
        # Regions extracted from each input chain on its own, as well as from the
        # concatenation. antpack numbers a two-domain string by labelling one domain H and
        # the other L, which is right for a Fab and wrong for anything multivalent: a
        # bivalent VHH has its second domain filed under L and compared against the
        # light-chain CDR index, where no VHH CDR exists. Comparing each submitted chain in
        # its own right costs one extra numbering per chain and removes that whole class of
        # misassignment. For a single-chain submission this is the same single extraction.
        parts = chain_seqs.get(qid) or []
        per_chain = (
            [regions]
            if len(parts) < 2
            else [extract(c, cl.category, cl.region_kind) for c in parts]
        )
        # ``parts`` are the units actually searched — the submitted chains, each already
        # broken into its variable domains where it had more than one. Comparing the same
        # units keeps the two arms talking about the same thing.
        out[qid] = {
            "category": cl.category,
            "region_kind": cl.region_kind,
            "category_confidence": cl.confidence,
            "length": cl.length,
            "category_evidence": cl.evidence,
            "regions": regions,
            "chain_regions": per_chain,
            "focus_chain": chain,
            "focus_region": focus,
            "n_chains": chains.get(qid, 1),
            "declared_class": declared.get(qid),
            "declared_class_match": compare_declared(declared.get(qid), cl.category),
        }
    return out


def screen(
    records: dict,
    cfg: Config | None = None,
    target_fasta: str | Path | None = None,
    target_metadata: str | Path | None = None,
    workdir: str | Path | None = None,
    skip_prior_art: bool = False,
    declared: dict | None = None,
) -> dict:
    """Screen a batch of submitted sequences.

    Args:
        records: ``{submission_id: sequence}``. A multi-chain submission joins its
            chains with ``:``, as in the Proteinbase submission template; the chains are
            concatenated before screening (no linker is inserted, so a Fab still
            classifies as a Fab rather than as an scFv).
        cfg: thresholds and database locations; ``Config.from_env()`` if omitted.
        target_fasta: curated known-binder FASTA for the target. Defaults to
            ``cfg.target_fasta`` or the packaged ``data/tnfa_binders.fasta``.
        target_metadata: optional CSV of per-reference annotations (class, source, ...).
        workdir: scratch directory; a temporary one is used if omitted.
        skip_prior_art: run only the target arm (useful when the public databases are not
            mounted, e.g. in a unit test).
        declared: optional ``{submission_id: molecule_class}`` as declared by the
            submitter. Cross-checked against the classifier and reported; never used in
            its place.

    Returns:
        ``{"results": {id: record}, "timing": {...}, "config": {...}, "target": {...}}``
        where each record carries the category, both arms' evidence, and the verdict.
    """
    cfg = cfg or Config.from_env()
    if cfg.db_root is None:
        cfg.db_root = resolve_db_root()
    # ``:`` separates the chains of a multi-chain submission, and the two halves of this
    # pipeline want different things from them. Classification reads the concatenation, so
    # a paired H+L still types as a Fab rather than as two loose chains; the searches run
    # per chain, because every reference is one chain and a concatenated query dilutes
    # coverage by the fraction of the molecule that is not the matching chain. See
    # :mod:`prescreen.chains`.
    index = split_submissions(records)
    records = index.concat
    chain_counts = {k: len(v) for k, v in index.by_submission.items()}
    if not records:
        return {"results": {}, "timing": {}, "config": cfg.as_dict(), "target": {}}

    tmp_holder = None
    if workdir is None:
        tmp_holder = tempfile.TemporaryDirectory(prefix="prescreen-run-")
        workdir = tmp_holder.name
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    timing: dict = {}
    t = time.time()
    chain_seqs = {
        sid: [index.queries[q] for q in qids]
        for sid, qids in index.by_submission.items()
    }
    annotations = _annotate(records, declared, chain_counts, chain_seqs)
    timing["classify_regions"] = round(time.time() - t, 2)

    data_dir = Path(__file__).parent / "data"
    target_path = target_fasta or cfg.target_fasta or (data_dir / "tnfa_binders.fasta")
    cdr_path = cfg.target_cdr_csv or (data_dir / "tnfa_binder_cdrs.csv")
    reference = TargetReference(
        target_path,
        metadata=target_metadata or cfg.target_metadata,
        cdr_csv=cdr_path if Path(cdr_path).exists() else None,
    )
    t = time.time()
    reference.ensure_annotations()
    timing["reference_annotation"] = round(time.time() - t, 2)

    t = time.time()
    # One record per CHAIN from here to the merge: the query ids are synthetic, so they
    # survive a submission name containing whitespace, which MMseqs2 would otherwise
    # truncate into a key that matches nothing.
    query_fasta = mmseqs.write_fasta(index.queries, workdir / "queries.fasta")
    (workdir / "query_map.tsv").write_text(
        "query_id\tsubmission_id\tchain\tdomain\tlength\n"
        + "".join("\t".join(str(c) for c in r) + "\n" for r in index.map_rows())
    )
    target_whole = merge_whole(
        reference.search_whole(
            query_fasta,
            workdir / "target",
            mmseqs_bin=cfg.mmseqs_bin,
            threads=cfg.threads,
        ),
        index,
    )
    timing["target_whole"] = round(time.time() - t, 2)

    t = time.time()
    target_region = {
        qid: _best_region(reference, ann)
        for qid, ann in annotations.items()
    }
    timing["target_region"] = round(time.time() - t, 2)

    t = time.time()
    substrings = {qid: reference.substring_match(seq) for qid, seq in records.items()}
    timing["target_substring"] = round(time.time() - t, 2)

    prior = {}
    if not skip_prior_art and not cfg.db_root and not cfg.allow_missing_arms:
        # Falling through here would run the target arm alone and report every submission
        # as having no prior art, which is indistinguishable from a real clean result.
        raise refdb.ReferenceDbMissing(
            "no reference database root found: set PRESCREEN_DB_ROOT (or pass --db-root) "
            "to the directory built by scripts/build_reference_dbs.sh, or run with "
            "--skip-prior-art to screen against the target set alone."
        )
    if not skip_prior_art and cfg.db_root:
        t = time.time()
        prior = priorart.merge_chains(
            priorart.search(index.queries, cfg, workdir / "priorart"), index, cfg
        )
        timing["prior_art"] = round(time.time() - t, 2)
        # Paratope-aware prior art: re-number the top hits' own sequences and compare
        # binding regions, so a shared framework is not mistaken for a known molecule.
        t = time.time()
        hit_cache: dict = {}
        for qid, slot in prior.items():
            ann = annotations[qid]
            slot["region"] = priorart.best_region_vs_hits(
                ann["chain_regions"],
                ann["category"],
                ann["region_kind"],
                slot.get("top") or [],
                cache=hit_cache,
            )
        timing["prior_art_regions"] = round(time.time() - t, 2)

    t = time.time()
    batches = _cluster.cluster_batch(records, cfg, workdir / "cluster")
    timing["batch_cluster"] = round(time.time() - t, 2)

    results: dict = {}
    for qid, seq in records.items():
        ann = annotations[qid]
        record = {
            "id": qid,
            "sequence": seq,
            "length": ann["length"],
            "category": ann["category"],
            "region_kind": ann["region_kind"],
            "category_confidence": ann["category_confidence"],
            "category_evidence": ann["category_evidence"],
            "focus_chain": ann["focus_chain"],
            "focus_region": ann["focus_region"],
            "n_chains": ann["n_chains"],
            "chain_lengths": index.lengths.get(qid, [ann["length"]]),
            "declared_class": ann["declared_class"],
            "declared_class_match": ann["declared_class_match"],
            "prior_art": prior.get(qid, {}),
            "target": {
                "whole": target_whole.get(qid),
                "region": target_region.get(qid),
                "substring": substrings.get(qid),
            },
            "batch": batches.get(qid, {}),
        }
        record.update(_flags.evaluate(record, cfg))
        results[qid] = record

    if tmp_holder is not None:
        tmp_holder.cleanup()
    return {
        "results": results,
        # Submissions whose sequence was empty once cleaned (":", whitespace, "***").
        # Named rather than dropped in silence: a submission missing from the report is
        # indistinguishable from one that passed.
        "dropped": index.dropped,
        "timing": timing,
        "config": cfg.as_dict(),
        "target": {
            "name": cfg.target_name,
            "fasta": str(target_path),
            "references": len(reference.sequences),
            "checksum": reference.checksum,
        },
    }


def screen_fasta(path: str | Path, **kwargs) -> dict:
    """Screen every sequence in a FASTA file. Keyword arguments go to :func:`screen`."""
    return screen(mmseqs.read_fasta(path), **kwargs)


def to_rows(screened: dict) -> list:
    """Flatten :func:`screen` output into one flat dict per submission, for CSV export."""
    rows = []
    for rec in screened["results"].values():
        pa_best = rec["prior_art"].get("best") or {}
        pa_design = rec["prior_art"].get("best_design") or {}
        pa_region = rec["prior_art"].get("region") or {}
        whole = (rec["target"].get("whole") or {})
        region = (rec["target"].get("region") or {})
        sub = rec["target"].get("substring") or {}
        rows.append(
            {
                "id": rec["id"],
                "length": rec["length"],
                "category": rec["category"],
                "category_confidence": rec["category_confidence"],
                "declared_class": rec["declared_class"],
                "declared_class_match": rec["declared_class_match"],
                "n_chains": rec["n_chains"],
                "verdict": rec["verdict"],
                "flags": ";".join(rec["flags"]),
                "context": ";".join(rec.get("context") or {}),
                "prior_art_similarity": round(pa_best.get("similarity_check", 0.0), 4),
                "prior_art_db": pa_best.get("db"),
                "prior_art_hit": pa_best.get("target"),
                "prior_art_region_identity": pa_region.get("region_identity"),
                "prior_art_arms_missing": ";".join(
                    rec["prior_art"].get("arms_missing") or []
                ),
                "prior_art_region_hit": pa_region.get("hit"),
                "design_similarity": round(pa_design.get("similarity_check", 0.0), 4),
                "design_hit": pa_design.get("target"),
                "target_similarity": round(whole.get("similarity_check", 0.0), 4),
                "target_hit": whole.get("target"),
                "target_region_identity": region.get("region_identity"),
                "target_region_compared": region.get("region_compared", True),
                "target_region_exact": region.get("region_exact", False),
                "target_region_hit": region.get("region_reference"),
                "target_binder_status": region.get("binder_status")
                or whole.get("binder_status"),
                "target_named_agent": region.get("named_agent")
                or whole.get("named_agent"),
                "focus_region": rec["focus_region"],
                "verbatim_fragment_of": sub.get("reference"),
                "cluster_representative": rec["batch"].get("cluster_representative"),
                "cluster_size": rec["batch"].get("cluster_size"),
                # Which INPUT chain (1-based, the n-th ``:``-separated part) produced each
                # hit. Not to be confused with target_region_chain, which is an IMGT chain
                # letter (H/L/S/W) within a single numbered domain.
                "chain_lengths": ";".join(str(n) for n in rec.get("chain_lengths") or []),
                "prior_art_chain": pa_best.get("chain"),
                "prior_art_domain": pa_best.get("domain"),
                "design_chain": pa_design.get("chain"),
                "prior_art_region_chain": pa_region.get("chain"),
                "target_chain": whole.get("chain"),
                "target_domain": whole.get("domain"),
                "target_region_chain": region.get("region_chain"),
            }
        )
    return rows


os.environ.setdefault("KMP_AFFINITY", "disabled")
