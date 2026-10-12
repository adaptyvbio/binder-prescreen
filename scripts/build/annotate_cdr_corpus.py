#!/usr/bin/env python3
"""Annotate a CDR index with how common each region is in a general antibody corpus.

The region arm's premise is that a CDR3 match means a shared paratope. That holds only
for regions that are actually rare. Kappa CDR-L3 in particular is germline-encoded and
9-10 residues long, so the same handful of endings recurs across antibodies against
completely unrelated targets - and an index that treats such a match as evidence will
call any antibody using that light chain a known binder of whatever the index is for.

Measured before this existed: ``QQRSNWPLT`` occurs in 358 of 350k PLAbDab antibodies, and
all 358 were called known TNF-alpha binders.

So each row gets a ``corpus_occurrences`` count, and :mod:`prescreen.cdrref` refuses to
let a region above ``region_corpus_max`` decide a verdict. The corpus should be a broad,
target-agnostic antibody set - PLAbDab is what this was measured on - precisely because
its *lack* of relation to the target is what makes a high count mean "germline".

Usage::

    annotate_cdr_corpus.py --cdrs <index.csv> --corpus <antibodies.fasta> [--jobs 8]
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

FIELD = "corpus_occurrences"


def corpus_lines(corpus: Path) -> Path:
    """Sequence lines only, so a header can never contribute a match."""
    tmp = Path(tempfile.mkstemp(prefix="cdrcorpus-", suffix=".txt")[1])
    with tmp.open("w") as out:
        subprocess.run(["grep", "-v", "^>", str(corpus)], stdout=out, check=True)
    return tmp


def count(pattern: str, path: Path) -> int:
    """Sequences containing this region.

    One grep per pattern rather than one multi-pattern pass: ``grep -o -F -f`` consumes
    overlapping matches, so short patterns swallow the text that longer ones would have
    matched. That undercounted ``QQRSNWPLT`` by 8x when this was first measured.

    An ``X``-containing region is a Markush claim and is counted with the ``X`` expanded
    to "any residue", because that is how :mod:`prescreen.cdrref` matches it. Counting it
    literally returns 0 - no real sequence contains ``X`` - which would mark the broadest
    entries in the index as the most distinctive ones. ``QQXXXXXXT`` matches any 9-mer
    kappa CDR-L3 of the form QQ...T; literally it scores 0, expanded it scores tens of
    thousands.
    """
    if "X" in pattern:
        proc = subprocess.run(
            ["grep", "-c", "-E", "--", pattern.replace("X", "."), str(path)],
            capture_output=True, text=True,
        )
    else:
        proc = subprocess.run(
            ["grep", "-c", "-F", "--", pattern, str(path)],
            capture_output=True, text=True,
        )
    return int((proc.stdout or "0").strip() or 0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cdrs", required=True, help="CDR index CSV to annotate in place")
    ap.add_argument("--corpus", required=True, help="general antibody FASTA")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", default=None, help="default: overwrite --cdrs")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.cdrs)))
    fields = [f for f in rows[0] if f != FIELD] + [FIELD]
    patterns = sorted({(r["cdr_seq"] or "").strip().upper() for r in rows if r["cdr_seq"].strip()})
    lines = corpus_lines(Path(args.corpus))
    try:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            counts = dict(zip(patterns, pool.map(lambda p: count(p, lines), patterns)))
    finally:
        lines.unlink(missing_ok=True)

    for r in rows:
        r[FIELD] = counts.get((r["cdr_seq"] or "").strip().upper(), 0)

    out = Path(args.out or args.cdrs)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    hot = sorted(rows, key=lambda r: -int(r[FIELD]))[:5]
    never = sum(1 for r in rows if int(r[FIELD]) == 0)
    print(f"annotated {len(rows)} rows -> {out}")
    print(f"  never seen in corpus: {never}  ({never * 100 // len(rows)}%)")
    print("  most common regions:")
    for r in hot:
        print(f"    {r['cdr_id']:14s} {r['cdr_seq']:16s} {r[FIELD]}")


if __name__ == "__main__":
    main()
