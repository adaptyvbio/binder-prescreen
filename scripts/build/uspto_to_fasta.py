"""Convert EBI's EMBL-style uspto_prt.dat into FASTA.

createdb has no parser for the EMBL flat format: running it on the .dat directly
ingests the whole corpus as a single entry, which is how the deployed patent
database came to hold one record against an 8.9 GB data file.
"""
import sys, re

src, dst = sys.argv[1], sys.argv[2]
n = 0
skipped_empty = 0
with open(src, "r", errors="replace") as fh, open(dst, "w") as out:
    acc = None
    desc = ""
    seq = []
    in_seq = False
    for line in fh:
        if line.startswith("ID   "):
            acc = line[5:].split()[0]
            desc, seq, in_seq = "", [], False
        elif line.startswith("DE   "):
            desc += (" " if desc else "") + line[5:].strip()
        elif line.startswith("SQ   "):
            in_seq = True
        elif line.startswith("//"):
            s = "".join(seq)
            if acc and s:
                out.write(f">uspto_{acc} {desc}\n")
                for i in range(0, len(s), 80):
                    out.write(s[i:i+80] + "\n")
                n += 1
            elif acc:
                skipped_empty += 1
            acc, desc, seq, in_seq = None, "", [], False
        elif in_seq and line.startswith("     "):
            seq.append(re.sub(r"[^A-Za-z]", "", line))
print(f"records={n} skipped_empty={skipped_empty}")
