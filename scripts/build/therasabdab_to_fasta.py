#!/usr/bin/env python3
"""Thera-SAbDab CSV -> FASTA, one record per chain.

The therapeutic antibodies, including every approved anti-TNF biologic. Each row is one
*agent*; a submission is a *chain*, so every chain is emitted separately.

The dump carries four chain columns, not two: a bispecific agent's second Fv lives in
`HeavySequence(ifbispec)` / `LightSequence(ifbispec)`. Reading only the first pair loses
~300 chains, several of them anti-TNF, so all four are read. Absent chains are the string
`na` rather than an empty cell.

Usage: therasabdab_to_fasta.py IN.csv OUT.fasta
"""

import csv
import sys

NAME_COLUMNS = ("Therapeutic", "therapeutic", "Name", "name")
# label -> candidate column names, first present wins
CHAIN_COLUMNS = {
    "H": ("HeavySequence", "Heavy Sequence", "heavy_sequence", "VH"),
    "L": ("LightSequence", "Light Sequence", "light_sequence", "VL"),
    "H2": ("HeavySequence(ifbispec)", "HeavySequence (ifbispec)"),
    "L2": ("LightSequence(ifbispec)", "LightSequence (ifbispec)"),
}
MIN_LENGTH = 50          # below this it is a placeholder, not a variable domain
AA = set("ACDEFGHIKLMNPQRSTVWYXBZU")


def pick(row, candidates):
    for c in candidates:
        if c in row:
            return c
    return None


def main(src, dst):
    csv.field_size_limit(1 << 24)
    written, agents = 0, 0
    with open(src, newline="", encoding="utf-8-sig") as fh, open(dst, "w") as out:
        reader = csv.DictReader(fh)
        cols = None
        for i, row in enumerate(reader):
            if cols is None:
                cols = {k: pick(row, v) for k, v in CHAIN_COLUMNS.items()}
                cols["name"] = pick(row, NAME_COLUMNS)
                if not any(cols[k] for k in CHAIN_COLUMNS):
                    sys.exit(f"no chain columns in {src}; saw {sorted(row)[:12]}")
            name = (row.get(cols["name"]) or f"thera{i + 1}").strip().replace(" ", "_")
            agents += 1
            for label in CHAIN_COLUMNS:
                col = cols.get(label)
                if not col:
                    continue
                seq = (row.get(col) or "").strip().upper().replace("-", "")
                if len(seq) < MIN_LENGTH or set(seq) - AA:
                    continue
                out.write(f">therasabdab_{name}_{label}\n{seq}\n")
                written += 1
    print(f"therasabdab: {written} chains from {agents} agents -> {dst}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
