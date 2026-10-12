#!/usr/bin/env python3
"""Normalise a curated known-binder FASTA into a target set the prescreen can load.

A curation track hands over one FASTA whose headers already carry provenance::

    >GLP1B00203|antibody_fab_chain|confirmed_binder|PDB|6LN2|Fab7F38_light_chain

Two things stand between that and a usable target set.

**The class vocabulary.** ``prescreen.target`` matches the second header field against a
fixed set of slugs (``TargetReference.ANTIBODY_REF_CLASSES``) to decide whether a
reference is represented by the CDR index or has to be region-annotated itself. A class
it does not recognise takes the wrong branch silently. Curation vocabularies differ per
target, so this script maps them over.

**What the mapping should be driven by.** Not the upstream label, which describes the
*molecule* ("antibody_or_fc_fusion_chain" covers both a real IgG chain and a peptide-Fc
fusion with no CDRs at all), but whether the sequence actually numbers as an antibody
domain — because that is the question the class is consulted to answer. So a header
class in ``--classify-classes`` is resolved by running the classifier over the sequence
and reading the chain type off the numbering. The upstream label is preserved verbatim
in the ``binder_class_raw`` column, so nothing is lost.

Writes ``<slug>_binders.fasta`` (headers rewritten to the mapped class) and
``<slug>_binders.csv`` beside each other, both ready for ``targets.json``.

Usage::

    target_set_from_fasta.py --fasta curated.fasta --slug mytarget --out-dir prescreen/data/targets/mytarget
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from prescreen import mmseqs  # noqa: E402
from prescreen.classify import classify  # noqa: E402
from prescreen.regions import extract  # noqa: E402
from prescreen.targets import MIN_SEARCHABLE_LENGTH  # noqa: E402

#: The columns downstream code actually reads. The curated TNF-alpha CSV carries 24; the
#: rest are provenance a converter cannot invent, so they are written empty rather than
#: filled with plausible-looking values.
FIELDS = [
    "id", "class", "binder_class_raw", "binder_status", "named_agent", "length",
    "chain_type", "source_db", "accession", "note", "sequence",
]

#: Straight renames: the upstream label already says what the molecule is.
DEFAULT_RENAME = {"peptide": "peptide_or_cdr"}

#: Labels too coarse to trust - resolved by numbering the sequence instead.
DEFAULT_CLASSIFY = ("antibody_fab_chain", "antibody_or_fc_fusion_chain")

#: Where a numbered chain lands, by the chain antpack assigned it.
BY_CHAIN = {"H": "igg_chain_heavy", "K": "igg_chain_light", "L": "igg_chain_light"}

#: A sequence in --classify-classes that does not number as an antibody. Deliberately
#: NOT an antibody class: it must be region-annotated like any other scaffold, and
#: putting it in ANTIBODY_REF_CLASSES would skip it as "already in the CDR index" when
#: it contributes no CDRs at all.
NOT_ANTIBODY = "other"


def _chain_of(seq: str) -> tuple[str | None, str]:
    """The antibody chain this sequence numbers as, plus the classifier's category."""
    cl = classify(seq)
    if cl.region_kind != "antibody_cdrs":
        return None, cl.category
    got = extract(seq, cl.category, cl.region_kind)
    chains = [c for c, labels in got.items() if any(v for k, v in labels.items() if k.startswith("cdr"))]
    return (chains[0] if chains else None), cl.category


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="curated binder FASTA")
    ap.add_argument("--slug", required=True, help="target slug, e.g. 'mytarget'")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument(
        "--classify-classes",
        default=",".join(DEFAULT_CLASSIFY),
        help="header classes resolved by numbering rather than by rename",
    )
    ap.add_argument(
        "--rename",
        default=",".join(f"{k}={v}" for k, v in DEFAULT_RENAME.items()),
        help="comma-separated old=new header class renames",
    )
    args = ap.parse_args()

    rename = dict(p.split("=", 1) for p in args.rename.split(",") if p)
    by_numbering = {c for c in args.classify_classes.split(",") if c}
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows, headers, mapped = [], [], Counter()
    short, ambiguous = [], []
    for header, seq in mmseqs.read_fasta(args.fasta).items():
        parts = header.split("|")
        rid = parts[0]
        raw = parts[1] if len(parts) > 1 else ""
        fields = {
            "binder_status": parts[2] if len(parts) > 2 else "",
            "source_db": parts[3] if len(parts) > 3 else "",
            "accession": parts[4] if len(parts) > 4 else "",
            "note": parts[5] if len(parts) > 5 else "",
        }
        chain = ""
        if raw in by_numbering:
            got, category = _chain_of(seq)
            if got is None:
                cls = NOT_ANTIBODY
                ambiguous.append(f"{rid} ({raw}, {fields['note']}) -> {cls}: numbers as {category}")
            else:
                cls = BY_CHAIN.get(got, NOT_ANTIBODY)
                chain = got
        else:
            cls = rename.get(raw, raw)
        mapped[f"{raw} -> {cls}"] += 1
        if len(seq) < MIN_SEARCHABLE_LENGTH:
            short.append(f"{rid} ({len(seq)} aa)")
        headers.append((f"{rid}|{cls}|" + "|".join(
            [fields["binder_status"], fields["source_db"], fields["accession"], fields["note"]]), seq))
        rows.append({
            "id": rid, "class": cls, "binder_class_raw": raw, "length": len(seq),
            "chain_type": chain, "named_agent": fields["note"], "sequence": seq, **fields,
        })

    fasta_out = out_dir / f"{args.slug}_binders.fasta"
    with fasta_out.open("w") as fh:
        for header, seq in headers:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    csv_out = out_dir / f"{args.slug}_binders.csv"
    with csv_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} references -> {fasta_out}")
    for label, n in sorted(mapped.items()):
        print(f"  {label}: {n}")
    for line in ambiguous:
        print(f"  numbering: {line}")
    if short:
        print(
            f"  WARNING: {len(short)} reference(s) under {MIN_SEARCHABLE_LENGTH} aa cannot be "
            f"matched by the MMseqs2 short profile or the substring fallback: {', '.join(short)}"
        )
    x = [r["id"] for r in rows if "X" in r["sequence"]]
    if x:
        print(f"  note: {len(x)} reference(s) contain X; they match only as wildcards: {', '.join(x[:8])}"
              + (" ..." if len(x) > 8 else ""))


if __name__ == "__main__":
    main()
