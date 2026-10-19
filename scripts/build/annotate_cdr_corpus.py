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

The same counting applies to a *binder* set, where the problem has the same shape: USPTO
sequence listings enumerate framework segments as separate entries, so a curated set can
end up holding ``GTLVTVSS`` (heavy-chain FR4, in 86,712 of 350,350 PLAbDab antibodies) as
if it were a binder. Pass ``--column sequence`` for that file.

An exact count is not enough for CDR-L3. The region arm matches at ``region_match``
identity (0.80), not exactly, so a reference CDR-L3 one or two edits from a germline ending
is as germline as the ending itself while its exact count stays low: the TfR1 set's
``QQRSNWPPIT`` occurs 30 times, sits at 0.80 from ``QQRSNWPLT``, and called 345 of 358
PLAbDab light chains carrying that ending known TfR1 binders. ``--near-identity 0.80``
counts, for light-chain CDR3 rows, every corpus antibody whose CDR-L3 is at least that
identical to the row, and records the larger of the two counts. Light chain only: CDR-H3
is VDJ-diverse, and the germline-ending problem is a kappa/lambda one.

Usage::

    annotate_cdr_corpus.py --cdrs <index.csv> --corpus <antibodies.fasta> [--jobs 8]
    annotate_cdr_corpus.py --cdrs <index.csv> --corpus <antibodies.fasta> --near-identity 0.80
    annotate_cdr_corpus.py --cdrs <binders.csv> --column sequence --corpus <antibodies.fasta>
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

FIELD = "corpus_occurrences"
#: CDR-L3 as IMGT bounds it: after the conserved Cys, before the J-region FGxG.
CDR_L3 = re.compile(r"C([^C]{5,14}?)FG.G")
LIGHT = {"K", "L"}


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


def corpus_light_cdr3(path: Path) -> Counter:
    """CDR-L3 of every corpus sequence that has one, counted once per sequence."""
    out: Counter = Counter()
    with path.open() as fh:
        for line in fh:
            m = CDR_L3.search(line.strip())
            if m:
                out[m.group(1)] += 1
    return out


_NEAR: dict = {}


def _near_init(corpus_l3: dict, cut: float) -> None:
    _NEAR["l3"], _NEAR["cut"] = corpus_l3, cut


def near_count(pattern: str) -> int:
    """Corpus sequences whose CDR-L3 is at least ``cut`` identical to ``pattern``.

    Identity is the region arm's own metric (``cdr_identity``), so "common" means common
    by the same yardstick the arm matches with. ``X`` is not expanded here; a Markush CDR
    keeps its exact (already wildcard-expanded) count.
    """
    from prescreen._proteintyper.cdr_novelty import cdr_identity

    cut = _NEAR["cut"]
    n = len(pattern)
    # identity >= cut needs the lengths within (1 - cut) of the longer one
    slack = int((1 - cut) * (n / cut)) + 1
    return sum(
        c
        for l3, c in _NEAR["l3"].items()
        if abs(len(l3) - n) <= slack and cdr_identity(pattern, l3) >= cut
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cdrs", required=True, help="CDR index CSV to annotate in place")
    ap.add_argument("--corpus", required=True, help="general antibody FASTA")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument(
        "--column", default="cdr_seq",
        help="column holding the sequence to count ('sequence' for a binder set)",
    )
    ap.add_argument(
        "--max-length", type=int, default=0,
        help="skip rows longer than this (0 = no limit). A reference too long to sit "
        "inside a submission cannot cause a spurious containment, so counting it is "
        "wasted work on a 5,000-row binder set.",
    )
    ap.add_argument(
        "--near-identity", type=float, default=0.0,
        help="also count corpus CDR-L3s at least this identical to a light-chain CDR3 "
        "row, and keep the larger count (0 = exact only). Use the region_match cut.",
    )
    ap.add_argument("--out", default=None, help="default: overwrite --cdrs")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.cdrs)))
    col = args.column
    if col not in rows[0]:
        raise SystemExit(f"{args.cdrs} has no {col!r} column; columns: {list(rows[0])}")
    fields = [f for f in rows[0] if f != FIELD] + [FIELD]
    def _seq(row: dict) -> str:
        return (row.get(col) or "").strip().upper()
    patterns = sorted({
        _seq(r) for r in rows
        if _seq(r) and (not args.max_length or len(_seq(r)) <= args.max_length)
    })
    lines = corpus_lines(Path(args.corpus))
    try:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            counts = dict(zip(patterns, pool.map(lambda p: count(p, lines), patterns)))
        near = {}
        if args.near_identity:
            l3_rows = sorted({
                _seq(r) for r in rows
                if r.get("cdr_region") == "CDR3"
                and (r.get("chain_type") or "").upper() in LIGHT
                and _seq(r) and "X" not in _seq(r)
            })
            corpus_l3 = corpus_light_cdr3(lines)
            print(f"near-identity {args.near_identity}: {len(l3_rows)} CDR-L3 rows vs "
                  f"{len(corpus_l3)} distinct corpus CDR-L3s", file=sys.stderr)
            with ProcessPoolExecutor(
                max_workers=args.jobs, initializer=_near_init,
                initargs=(dict(corpus_l3), args.near_identity),
            ) as pool:
                near = dict(zip(l3_rows, pool.map(near_count, l3_rows, chunksize=4)))
    finally:
        lines.unlink(missing_ok=True)
    raised = 0
    for pattern, n in near.items():
        if n > counts.get(pattern, 0):
            counts[pattern] = n
            raised += 1

    for r in rows:
        # A row skipped by --max-length keeps whatever it had rather than being reset to
        # 0, which would read as "measured, never seen" instead of "not measured".
        if _seq(r) in counts:
            r[FIELD] = counts[_seq(r)]
        else:
            r.setdefault(FIELD, "")

    out = Path(args.out or args.cdrs)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    scored = [r for r in rows if str(r.get(FIELD, "")).strip() != ""]
    hot = sorted(scored, key=lambda r: -int(r[FIELD]))[:5]
    never = sum(1 for r in scored if int(r[FIELD]) == 0)
    ident = "cdr_id" if "cdr_id" in rows[0] else "id"
    print(f"annotated {len(scored)} of {len(rows)} rows -> {out}")
    if near:
        print(f"  raised by near-identity counting: {raised} of {len(near)} CDR-L3 rows")
    if scored:
        print(f"  never seen in corpus: {never}  ({never * 100 // len(scored)}%)")
    print("  most common:")
    for r in hot:
        print(f"    {r[ident]:14s} {_seq(r)[:34]:36s} {r[FIELD]}")


if __name__ == "__main__":
    main()
