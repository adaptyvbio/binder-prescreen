#!/usr/bin/env python3
"""Add the binding regions the CDR index is missing, from the curated binder set itself.

``tnfa_binder_cdrs.csv`` was built only from the classes that look like conventional
antibody chains — ``fv_domain_vh/vl``, ``vhh``, ``vhh_putative``, ``igg_chain_heavy/light``,
``fab``, ``igg_fv``. Four curated classes contribute no rows at all, and the region arm
compares a submission's CDR3 against this index and nothing else. A binding region that is
not in it cannot be recognised however exactly it is copied: the arm reports a low identity
because it has nothing to match, and that silence reads as "the paratope is new".

Two kinds of reference are missing, and they need different treatment:

``peptide_or_cdr`` (1,626 entries, 8-50 aa, almost all from USPTO sequence listings)
    These *are* binding regions — claimed loops and peptides, with no parent chain and so
    no chain or CDR-number to assign. They are emitted as wildcard entries (``*``/``*``)
    that match a query's region whatever chain and position it was numbered at, which is
    the right semantics for a bare claimed loop: it is prior art wherever it reappears.
    Length-normalised identity keeps this honest — a 40-aa entry simply cannot score highly
    against a 14-aa CDR3.

``vnar``, ``other``, ``designed_other``
    These number as antibody domains, so their real CDR1/2/3 and chain types are extracted
    exactly as the original index's were.

Existing rows are never touched: this appends, and skips any sequence already present.

Usage: augment_binder_cdrs.py [--binders FASTA] [--metadata CSV] [--cdrs CSV] [--out CSV]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from prescreen import mmseqs  # noqa: E402
from prescreen.classify import classify  # noqa: E402
from prescreen.regions import extract  # noqa: E402

FIELDS = [
    "cdr_id", "cdr_seq", "cdr_region", "chain_type", "length", "n_parent_sequences",
    "best_status", "named_agents", "parent_classes", "source_dbs", "parent_ids",
    "example_citation",
]

#: Bare claimed loops shorter than this are too generic to be evidence on their own.
MIN_REGION = 8
WILDCARD_CLASSES = {"peptide_or_cdr"}
NUMBERED_CLASSES = {"vnar", "other", "designed_other"}


def covered_classes(cdr_rows: list) -> set:
    out = set()
    for row in cdr_rows:
        for parent in (row.get("parent_classes") or "").split(";"):
            if parent.strip():
                out.add(parent.strip())
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parents[2] / "prescreen" / "data"
    ap.add_argument("--binders", default=str(here / "tnfa_binders.fasta"))
    ap.add_argument("--metadata", default=str(here / "tnfa_binders.csv"))
    ap.add_argument("--cdrs", default=str(here / "tnfa_binder_cdrs.csv"))
    ap.add_argument("--out", default=None, help="default: overwrite --cdrs")
    args = ap.parse_args()

    seqs = {h.split("|")[0]: s for h, s in mmseqs.read_fasta(args.binders).items()}
    meta = {r["id"]: r for r in csv.DictReader(open(args.metadata))}
    existing = list(csv.DictReader(open(args.cdrs)))
    have = {(r.get("cdr_seq") or "").strip().upper() for r in existing}
    covered = covered_classes(existing)
    print(f"existing rows: {len(existing)}  covering classes: {sorted(covered)}")

    added, skipped_dupe, skipped_short, failed = [], 0, 0, 0
    n = 0
    for rid, m in meta.items():
        cls = m.get("class") or ""
        if cls in covered or rid not in seqs:
            continue
        seq = seqs[rid].strip().upper()
        base = {
            "best_status": m.get("binder_status") or "",
            "named_agents": m.get("named_agent") or "",
            "parent_classes": cls,
            "source_dbs": m.get("source_db") or "",
            "parent_ids": rid,
            "example_citation": m.get("citation") or "",
            "n_parent_sequences": 1,
        }
        regions = []
        if cls in WILDCARD_CLASSES:
            # No parent chain, so no chain or CDR number can be assigned honestly.
            regions = [("*", "*", seq)]
        elif cls in NUMBERED_CLASSES:
            try:
                cl = classify(seq)
                got = extract(seq, cl.category, cl.region_kind)
            except Exception:
                failed += 1
                continue
            if cl.region_kind != "antibody_cdrs":
                continue
            for chain, labels in got.items():
                for label, value in labels.items():
                    if label.startswith("cdr") and value:
                        regions.append((label.upper(), chain, value))
        for region, chain, value in regions:
            value = value.strip().upper()
            if len(value) < MIN_REGION:
                skipped_short += 1
                continue
            if value in have:
                skipped_dupe += 1
                continue
            have.add(value)
            n += 1
            added.append(
                {**base, "cdr_id": f"TNFCDRX{n:05d}", "cdr_seq": value,
                 "cdr_region": region, "chain_type": chain, "length": len(value)}
            )

    out_path = Path(args.out or args.cdrs)
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(existing)
        w.writerows(added)
    print(f"added {len(added)} rows "
          f"(duplicates skipped {skipped_dupe}, too short {skipped_short}, "
          f"not numberable {failed})")
    import collections
    print("by class:", dict(collections.Counter(r["parent_classes"] for r in added)))
    print("wrote", out_path, "->", len(existing) + len(added), "rows")


if __name__ == "__main__":
    main()
