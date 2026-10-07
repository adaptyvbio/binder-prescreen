#!/usr/bin/env python3
"""PLAbDab-nano CSV -> FASTA of unique single-domain sequences.

PLAbDab's paired dump has a heavy AND a light chain for every entry, so it contains no
single-domain antibodies at all. The nano dump is the VHH / VNAR / sdAb prior art, and it
is a separate arm for that reason.

Usage: plabdab_nano_to_fasta.py IN.csv OUT.fasta
"""

import csv
import sys

# The dump has moved column names between releases; take the first that is present.
SEQ_COLUMNS = ("sequence", "sequence_alignment_aa", "heavy_sequence", "vh")
ID_COLUMNS = ("id", "entry_id", "name", "pdb")


def pick(row, candidates):
    for c in candidates:
        if row.get(c):
            return c
    return None


def main(src, dst):
    csv.field_size_limit(1 << 24)
    seen = {}
    with open(src, newline="") as fh:
        reader = csv.DictReader(fh)
        seq_col = id_col = None
        for i, row in enumerate(reader):
            if seq_col is None:
                seq_col = pick(row, SEQ_COLUMNS)
                id_col = pick(row, ID_COLUMNS)
                if seq_col is None:
                    sys.exit(
                        f"no sequence column in {src}; saw {sorted(row)[:12]}"
                    )
            seq = (row.get(seq_col) or "").strip().upper().replace("-", "")
            if not seq:
                continue
            # one record per unique sequence: the dump repeats a sequence across entries
            seen.setdefault(seq, (row.get(id_col) or f"nano{i + 1}").strip())

    with open(dst, "w") as out:
        for n, (seq, name) in enumerate(seen.items(), 1):
            out.write(f">{name or f'nano{n}'}\n{seq}\n")
    print(f"plabdab_nano: {len(seen)} unique sequences -> {dst}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
