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

With ``--cdrs`` pointing at a file that does not exist the same machinery builds an index
from scratch, which is how a newly curated target gets one. Widen ``--numbered-classes``
in that case: the default set assumes the conventional antibody classes are already
covered by an existing index, and from scratch nothing is.

Usage: augment_binder_cdrs.py [--binders FASTA] [--metadata CSV] [--cdrs CSV] [--out CSV]
                              [--prefix TNFCDRX] [--numbered-classes a,b,c]
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

#: Residues a linker is made of. A region drawn entirely from these is a (G4S)n-style
#: spacer that the numbering mistook for a CDR, not a paratope: scFv and multi-domain
#: constructs in patent listings put their linker exactly where a CDR is expected. Such
#: an entry is pure false-positive surface - it matches any poly-glycine control and
#: nothing meaningful - so it never enters the index.
LINKER_RESIDUES = frozenset("GS")


def is_linker(seq: str) -> bool:
    """Whether this region is a glycine/serine spacer rather than a binding region."""
    return bool(seq) and set(seq.upper()) <= LINKER_RESIDUES
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
    here = Path(__file__).resolve().parents[2] / "prescreen" / "data" / "targets" / "tnfa"
    ap.add_argument("--binders", default=str(here / "tnfa_binders.fasta"))
    ap.add_argument("--metadata", default=str(here / "tnfa_binders.csv"))
    ap.add_argument("--cdrs", default=str(here / "tnfa_binder_cdrs.csv"))
    ap.add_argument("--out", default=None, help="default: overwrite --cdrs")
    ap.add_argument(
        "--prefix",
        default="TNFCDRX",
        help="id prefix for the rows this run mints; must not collide with the ids "
        "already in --cdrs, or two different regions end up sharing one id",
    )
    ap.add_argument(
        "--numbered-classes",
        default=",".join(sorted(NUMBERED_CLASSES)),
        help="classes whose CDRs are extracted by numbering",
    )
    ap.add_argument(
        "--wildcard-classes",
        default=",".join(sorted(WILDCARD_CLASSES)),
        help="classes emitted as chain-agnostic wildcard entries",
    )
    args = ap.parse_args()

    numbered = {c for c in args.numbered_classes.split(",") if c}
    wildcard = {c for c in args.wildcard_classes.split(",") if c}

    seqs = {h.split("|")[0]: s for h, s in mmseqs.read_fasta(args.binders).items()}
    meta = {r["id"]: r for r in csv.DictReader(open(args.metadata))}
    # Absent is not an error: that is the from-scratch case for a new target.
    existing = list(csv.DictReader(open(args.cdrs))) if Path(args.cdrs).exists() else []
    purged = [r for r in existing if is_linker((r.get("cdr_seq") or "").strip())]
    if purged:
        existing = [r for r in existing if r not in purged]
        print(f"purged {len(purged)} linker row(s) already in the index: "
              + ", ".join(f"{r['cdr_id']}={r['cdr_seq']}" for r in purged))
    have = {(r.get("cdr_seq") or "").strip().upper() for r in existing}
    covered = covered_classes(existing)
    print(f"existing rows: {len(existing)}  covering classes: {sorted(covered)}")

    added, skipped_dupe, skipped_short, skipped_linker, failed = [], 0, 0, 0, 0
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
        if cls in wildcard:
            # No parent chain, so no chain or CDR number can be assigned honestly.
            regions = [("*", "*", seq)]
        elif cls in numbered:
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
            if is_linker(value):
                skipped_linker += 1
                continue
            if value in have:
                skipped_dupe += 1
                continue
            have.add(value)
            n += 1
            added.append(
                {**base, "cdr_id": f"{args.prefix}{n:05d}", "cdr_seq": value,
                 "cdr_region": region, "chain_type": chain, "length": len(value)}
            )

    clash = {r.get("cdr_id", "") for r in existing} & {r["cdr_id"] for r in added}
    if clash:
        sys.exit(
            f"--prefix {args.prefix} collides with {len(clash)} existing id(s); "
            "pick another or the index silently maps one id onto two regions"
        )

    out_path = Path(args.out or args.cdrs)
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(existing)
        w.writerows(added)
    print(f"added {len(added)} rows "
          f"(duplicates skipped {skipped_dupe}, too short {skipped_short}, "
          f"linkers {skipped_linker}, not numberable {failed})")
    import collections
    print("by class:", dict(collections.Counter(r["parent_classes"] for r in added)))
    print("wrote", out_path, "->", len(existing) + len(added), "rows")


if __name__ == "__main__":
    main()
