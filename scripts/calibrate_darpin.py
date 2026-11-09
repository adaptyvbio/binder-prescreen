"""Re-measure the DARPin segment rule (prescreen/darpin.py) against labelled sets.

Positives: PDB chains whose entity name says DARPin / designed ankyrin repeat, plus a few
generically named DARPins. Negatives: reviewed UniProt entries with an ankyrin-repeat
domain (InterPro IPR002110) and PDB chains named "ankyrin" that are not designed.

    python scripts/calibrate_darpin.py WORKDIR

Downloads pdb_seqres.txt.gz (~70 MB) and the UniProt set into WORKDIR, then prints the
segment-score range per class and how the old single-reference gate compares.
"""

from __future__ import annotations

import gzip
import re
import sys
import urllib.request
from pathlib import Path

from prescreen import darpin
from prescreen.compat import proteintyper

SEQRES = "https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
UNIPROT = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    "query=(xref:interpro-IPR002110)+AND+(reviewed:true)&format=fasta"
)
DESIGNED = re.compile(r"darpin|designed ankyrin|SANK E3|consensus ankyrin|DARP", re.I)
NOT_NATURAL = re.compile(r"darpin|designed|synthetic|SANK|consensus|engineered|DARP", re.I)
#: Name-filter mistakes, checked by hand: entities matching "DARP" that are not DARPins
#: (the DARP14 cage subunits, the deltaDARP LOV2 fusion) ...
NOT_DARPIN = re.compile(r"DARP14|deltaDARP", re.I)
#: ... DARPins with a generic entity name ...
GENERIC_DARPIN = {"1svx_A", "2qyj_A", "2v5q_C", "6ob5_C", "9kbu_A"}
#: ... and designed ankyrins from other libraries, reported but not scored.
OTHER_DESIGNED = {"4hll_A", "4zfh_A", "5gik_A", "1n0r_A", "5kwd_A"}


def fetch(url: str, path: Path) -> Path:
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    return path


def read_fasta(text: str):
    name, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main(workdir: str) -> None:
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    seqres = gzip.open(fetch(SEQRES, work / "pdb_seqres.txt.gz"), "rt").read()
    uniprot = fetch(UNIPROT, work / "uniprot_ankyrin.fasta").read_text()

    labelled: dict[str, tuple[str, str]] = {}
    seen: set[str] = set()
    for name, seq in read_fasta(seqres):
        if "mol:protein" not in name or seq in seen:
            continue
        cid = name.split()[0]
        if cid in GENERIC_DARPIN or (
            DESIGNED.search(name) and not NOT_DARPIN.search(name)
        ):
            labelled[cid] = ("darpin", seq)
        elif cid in OTHER_DESIGNED:
            labelled[cid] = ("other_designed", seq)
        elif re.search("ankyrin", name, re.I) and not NOT_NATURAL.search(name):
            labelled[cid] = ("natural", seq)
        else:
            continue
        seen.add(seq)
    for name, seq in read_fasta(uniprot):
        labelled[name.split()[0]] = ("natural", seq)

    scaffold_cdr = proteintyper("scaffold_cdr")
    rows: dict[str, list] = {}
    for cid, (label, seq) in labelled.items():
        seq = "".join(c for c in seq.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
        if len(seq) < 40:
            continue
        ev = darpin.evidence(seq)
        old = (
            len(seq) >= int(0.6 * len(scaffold_cdr.SCAFFOLD_REFS["darpin"]))
            and scaffold_cdr.framework_identity("darpin", seq) >= 0.85
        )
        rows.setdefault(label, []).append((ev["segment_score"], darpin.is_darpin(ev), old, cid))

    for label, rs in sorted(rows.items()):
        scores = sorted(r[0] for r in rs)
        print(
            f"{label:15} n={len(rs):5}  segment score {scores[0]:.3f}-{scores[-1]:.3f}  "
            f"new rule {sum(r[1] for r in rs):5}  old gate {sum(r[2] for r in rs):5}"
        )
        edge = max(rs) if label != "darpin" else min(rs)
        print(f"{'':15} closest to the boundary: {edge[3]} ({edge[0]:.3f})")


if __name__ == "__main__":
    main(sys.argv[1])
