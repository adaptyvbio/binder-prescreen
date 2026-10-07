#!/usr/bin/env python3
"""THPdb TSV -> FASTA of FDA-approved peptide and protein therapeutics.

The THPdb web host (webs.iiitd.edu.in) serves a ready-made FASTA but is frequently
unreachable, so the build pulls the Figshare deposit instead:

    https://figshare.com/articles/dataset/5198005   (file 8868913, thpdb.txt)

That file is a 39-column TSV, not a FASTA, and needs three things doing to it.

**Rows repeat.** 852 rows cover 239 therapeutics — one row per annotation, all carrying
the same sequence cell. Deduplicated on the sequence itself.

**Multi-chain entries pack every chain into the sequence cell**, labelled inline. Written
out whole, such a cell is a chimera that matches nothing. There are 62 distinct label
forms across the 161 unique cells — `Heavy chain:`, `A-chain:`, `Light Chain 2`,
`Alpha chain :`, `H-GAMMA-1 Chain:`, `L-KAPPA chain:`, `Beta-Chain(L` — differing in case,
in separator (space, hyphen, `;`, `,`, or nothing at all), in whether the colon is there,
and in whether a product name precedes them (`Ofatumumab Heavy Chain`). So the label is
matched loosely and each resulting segment is then reduced to its longest run of valid
residues, which drops any product name left behind.

**Some cells are not sequences.** Three-letter notation (`MPR-HAR-GLY-...`), prose
("ReoPro-like antibody"), and `NA` all appear; they are skipped.

Every letter in `LIGHTCHAIN` is a valid residue code, so a label that survives into a
sequence passes an amino-acid check silently. :func:`assert_no_label_leak` is what makes
that failure loud instead.

Usage: thpdb_to_fasta.py IN.tsv OUT.fasta
"""

import csv
import re
import sys

ID_COL, NAME_COL, SEQ_COL = 1, 2, 3
MIN_LENGTH = 5
AA = set("ACDEFGHIKLMNPQRSTVWYXBZU")

# The qualifier is matched only when it is one of the known words (or a lone capital),
# because those are the forms that appear glued to the previous chain's last residue
# ("...NRGECHeavy-chain:"). A product-name qualifier is always space-separated and is
# removed later by longest_residue_run().
CHAIN_LABEL = re.compile(
    r"(?:heavy|light|alpha|beta|gamma|delta|kappa|lambda|truncated|[A-Z](?![a-z]))?"
    r"[\s-]*chains?\s*\d*\s*[:(]?",
    re.IGNORECASE,
)


def longest_residue_run(segment):
    """The longest whitespace-delimited all-residue token in `segment`.

    Residue runs never contain whitespace in this file, so anything left beside the
    sequence — a product name, a stray word — is a separate token and loses.
    """
    best = ""
    for token in re.split(r"[^A-Za-z]+", segment):
        token = token.upper()
        if len(token) >= MIN_LENGTH and not set(token) - AA and len(token) > len(best):
            best = token
    return best


def split_chains(cell):
    """Yield (chain_label, segment) for every chain in one sequence cell."""
    cell = cell.strip()
    if not cell:
        return
    hits = list(CHAIN_LABEL.finditer(cell))
    hits = [m for m in hits if m.group(0).strip(" -:(")]  # drop empty matches
    if not hits:
        yield "", cell
        return
    if hits[0].start() > 0:
        yield "", cell[: hits[0].start()]
    for i, m in enumerate(hits):
        end = hits[i + 1].start() if i + 1 < len(hits) else len(cell)
        label = re.sub(r"[^A-Za-z0-9]+", "_", m.group(0).strip(" -:(")).strip("_")
        yield label, cell[m.end() : end]


def assert_no_label_leak(seq, name):
    """A label that survived into a sequence is invisible to an amino-acid check."""
    if "CHAIN" in seq:
        i = seq.find("CHAIN")
        raise SystemExit(
            f"chain label leaked into {name}: ...{seq[max(0, i - 30):i + 20]}... "
            "- CHAIN_LABEL does not cover a label form present in this file"
        )


def main(src, dst):
    csv.field_size_limit(1 << 24)
    seen, rows, skipped = {}, 0, 0
    with open(src, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) <= SEQ_COL:
                continue
            rows += 1
            tid = (row[ID_COL] or "").strip()
            name = re.sub(r"[^A-Za-z0-9]+", "_", (row[NAME_COL] or "").strip()).strip("_")
            for label, segment in split_chains(row[SEQ_COL]):
                seq = longest_residue_run(segment)
                if len(seq) < MIN_LENGTH:
                    skipped += 1
                    continue
                record = f"thpdb_{tid}_{name}" + (f"_{label}" if label else "")
                assert_no_label_leak(seq, record)
                seen.setdefault(seq, record)

    with open(dst, "w") as out:
        for seq, name in seen.items():
            out.write(f">{name}\n{seq}\n")
    ids = {n.split("_")[1] for n in seen.values()}
    print(
        f"thpdb: {len(seen)} unique chains covering {len(ids)} therapeutics, "
        f"from {rows} rows -> {dst} ({skipped} non-sequence segments skipped)"
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
