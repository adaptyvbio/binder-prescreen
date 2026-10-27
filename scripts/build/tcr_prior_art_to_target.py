#!/usr/bin/env python3
"""Build a pMHC target set from a merged TCR prior-art table.

A pMHC target's prior art is mostly TCRs, gathered from places that describe them
differently: IEDB and VDJdb report a receptor by its CDR3 pair, the PDB and patent sequence
listings by its chains (often with a signal peptide and a constant domain). The input is
one merged table with a row per binder::

    set,tier,modality,target_peptide,allele,binder_id,seq_type,sequence,cdr3a,cdr3b,
    interface_residues,pdb,reference,note

where ``sequence`` holds the chains of one receptor joined by ``|`` in any order, and the
CDR3s may carry the conserved C and F/W or not. Rows are kept when ``allele`` names
``--allele`` and ``target_peptide`` names one of the ``--mutation`` s (as ``G12D`` or as a
peptide in ``--peptides``) and the row carries a sequence or a CDR3 - so other peptides on
the same allele, and literature rows with nothing extracted, are left out.

Chains are read by TCR numbering (:mod:`prescreen.tcr`), not by position or by the
table's CDR3 columns: the V domain is cut out, signal peptide and constant domain dropped,
and its CDR1-3 come from the IMGT numbering. A sequence that does not number as a TCR is
reported and skipped. The table's CDR3s are used where a row has no sequence.

Writes three files beside each other, ready for ``targets.json``:

- ``<slug>_binders.fasta`` / ``.csv`` - one reference per unique V domain
  (``tcr_chain_alpha`` / ``tcr_chain_beta``) and per unique CDR3 (``tcr_cdr3``). A CDR3 is a
  whole reference too, so a verbatim CDR3 inside any submission is found by containment;
- ``<slug>_binder_cdrs.csv`` - the CDR index, chain-typed ``A``/``B``, with the receptors
  (table rows) each CDR comes from in ``receptor_ids``, so a query's alpha and beta CDR3
  can be required to match the same receptor.

Binder status: a binder of a ``--confirmed-mutation`` with a structure, or an IEDB record
(a positive T-cell or multimer assay), is ``confirmed_binder``; everything else - patent
claims, VDJdb-only receptors, and binders shown only for another mutation - is
``disclosed_binder_claim``, with the mutation in the note.

Usage::

    tcr_prior_art_to_target.py --table prior_art_merged_core.csv --slug kras-pmhc \\
        --prefix KRP --out-dir prescreen/data/targets/kras-pmhc
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from prescreen.tcr import tcr_domains  # noqa: E402

FIELDS = [
    "id", "class", "binder_class_raw", "binder_status", "named_agent", "length",
    "chain_type", "source_db", "accession", "note", "sequence", "citation", "evidence",
    "target_species", "all_sources", "all_accessions", "corpus_occurrences",
    "receptor_ids",
]
CDR_FIELDS = [
    "cdr_id", "cdr_seq", "cdr_region", "chain_type", "length", "n_parent_sequences",
    "best_status", "named_agents", "parent_classes", "source_dbs", "parent_ids",
    "example_citation", "corpus_occurrences", "receptor_ids",
]
CHAIN_CLASS = {"A": "tcr_chain_alpha", "B": "tcr_chain_beta"}
STATUS_RANK = {"confirmed_binder": 2, "disclosed_binder_claim": 1}
#: KRAS codon-12 peptides presented by HLA-A*11:01 (9- and 10-mers).
DEFAULT_PEPTIDES = {
    "G12D": ("VVVGADGVGK", "VVGADGVGK"),
    "G12V": ("VVVGAVGVGK", "VVGAVGVGK"),
}


def _cdr3(raw: str) -> str:
    """IMGT CDR3 (105-117): drop the conserved C104 and F/W118 when they are present."""
    s = (raw or "").strip().upper()
    if len(s) > 4 and s[0] == "C" and s[-1] in "FW":
        s = s[1:-1]
    return s


def _pmids(raw: str) -> list:
    out: list = []
    for tok in re.findall(r"(?:PMID:)?(\d{7,8})\b", raw or ""):
        if tok not in out:
            out.append(tok)
    return out


def _mutations(row: dict, wanted: dict) -> list:
    text = (row.get("target_peptide") or "").upper()
    return [m for m, peps in wanted.items() if m in text or any(p in text for p in peps)]


def _identity(row: dict, pdb_names: dict) -> dict:
    """Name, sources, accession and citation of one table row."""
    kind = row["set"].split()[0].upper()
    note = row.get("note") or ""
    if kind == "PDB":
        entry = row["binder_id"].split("|")[0].upper()
        m = re.search(r"(?:of|Structure of) (?:the )?(\S+) (?:TCR|Fab)\b", row.get("reference") or "")
        names = pdb_names.get(entry) or ([m.group(1)] if m and m.group(1) != "a" else [])
        return {"names": names, "sources": ["pdb"], "accessions": [entry],
                "citation": "", "evidence": "structure_complex", "kind": "pdb"}
    if kind == "PATENT":
        patent = row["binder_id"].split("_")[0]
        seq_ids = re.search(r"SEQ IDs? ([\d,\- ]+)", note)
        acc = f"{patent} SEQ ID {seq_ids.group(1).strip()}" if seq_ids else patent
        return {"names": [], "sources": ["patent"], "accessions": [acc],
                "citation": patent, "evidence": "patent_sequence_listing", "kind": "patent"}
    m = re.search(r"sources ([A-Za-z, ]+)", note)
    sources = [s.strip().lower() for s in (m.group(1) if m else "").split(",") if s.strip()]
    names = [n.strip() for n in row["binder_id"].split(",")
             if n.strip() and not n.startswith("PMID")]
    pmids = _pmids(row.get("reference"))
    return {"names": names, "sources": sources,
            "accessions": [f"{s}:{names[0]}" for s in sources] if names else sources,
            "citation": f"PMID:{pmids[0]}" if pmids else "",
            "evidence": "iedb_tcell_assay" if "iedb" in sources else "vdjdb_record",
            "kind": "db"}


class Ref:
    """One unique reference sequence, accumulating every binder that carries it."""

    def __init__(self, seq: str, cls: str, chain: str):
        self.seq, self.cls, self.chain = seq, cls, chain
        self.names, self.sources, self.accessions = [], [], []
        self.citations, self.notes, self.evidence = [], [], []
        self.status = "disclosed_binder_claim"
        self.cdrs: dict = {}
        #: The table rows (receptors) carrying this sequence. The CDR index records them
        #: so the region arm can ask whether both chains of a query match ONE receptor.
        self.receptors: list = []

    def add(self, ident: dict, status: str, note: str, receptor: str) -> None:
        if receptor not in self.receptors:
            self.receptors.append(receptor)
        for bucket, vals in (
            (self.names, ident["names"]), (self.sources, ident["sources"]),
            (self.accessions, ident["accessions"]), (self.citations, [ident["citation"]]),
            (self.notes, [note]), (self.evidence, [ident["evidence"]]),
        ):
            bucket.extend(v for v in vals if v and v not in bucket)
        if STATUS_RANK[status] > STATUS_RANK[self.status]:
            self.status = status


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--table", required=True, help="merged prior-art CSV")
    ap.add_argument("--slug", required=True)
    ap.add_argument("--prefix", required=True, help="id prefix, e.g. KRP")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--allele", default="11:01")
    ap.add_argument("--mutation", nargs="+", default=["G12D", "G12V"])
    ap.add_argument("--confirmed-mutation", default="G12D",
                    help="the target's own mutation; binders shown only for another "
                    "are kept as claims")
    args = ap.parse_args()
    wanted = {m: DEFAULT_PEPTIDES.get(m, ()) for m in args.mutation}

    rows = list(csv.DictReader(open(args.table, newline="", encoding="utf-8-sig")))
    # A database row naming a PDB entry lends the structure its name (7PB2 is JDI).
    pdb_names: dict = {}
    for row in rows:
        for p in re.split(r"[,; ]+", row.get("pdb") or ""):
            if p and row["set"].split()[0].upper() not in ("PDB", "PATENT"):
                pdb_names.setdefault(p.upper(), []).extend(
                    n.strip() for n in row["binder_id"].split(",") if n.strip())

    refs: dict = {}

    def ref(seq, cls, chain) -> Ref:
        return refs.setdefault((seq, cls), Ref(seq, cls, chain))

    kept, skipped = 0, []
    for row in rows:
        if args.allele not in (row.get("allele") or ""):
            continue
        muts = _mutations(row, wanted)
        seqs = [s.strip().upper() for s in (row.get("sequence") or "").split("|") if s.strip()]
        cdr3s = {"A": _cdr3(row.get("cdr3a")), "B": _cdr3(row.get("cdr3b"))}
        if not muts or not (seqs or any(cdr3s.values())):
            continue
        ident = _identity(row, pdb_names)
        confirmed = args.confirmed_mutation in muts and ident["kind"] == "pdb" or (
            args.confirmed_mutation in muts and "iedb" in ident["sources"]
        )
        status = "confirmed_binder" if confirmed else "disclosed_binder_claim"
        note = f"{'/'.join(muts)} on HLA-A*{args.allele}" + (
            "" if args.confirmed_mutation in muts
            else f"; {args.confirmed_mutation} not shown"
        ) + (f"; {row['note']}" if ident["kind"] == "db" and row.get("note") else "")

        numbered = {}
        for s in seqs:
            for d in tcr_domains(s):
                if d["chain"] in CHAIN_CLASS and d["chain"] not in numbered:
                    numbered[d["chain"]] = (s[d["start"]:d["end"]], d)
        if seqs and not numbered:
            skipped.append(f"{row['binder_id']}: no TCR V domain in its sequence")
            continue
        kept += 1
        for chain in ("A", "B"):
            if chain in numbered:
                v, d = numbered[chain]
                r = ref(v, CHAIN_CLASS[chain], chain)
                r.add(ident, status, note, row["binder_id"])
                r.cdrs = {k: d[k] for k in ("cdr1", "cdr2", "cdr3") if d[k]}
                if cdr3s[chain] and cdr3s[chain] != d["cdr3"]:
                    print(f"  {row['binder_id']} {chain}: table CDR3 {cdr3s[chain]} != "
                          f"numbered {d['cdr3']}; using the numbering", file=sys.stderr)
                cdr3 = d["cdr3"]
            else:
                cdr3 = cdr3s[chain]
            if cdr3:
                ref(cdr3, "tcr_cdr3", chain).add(ident, status, note, row["binder_id"])

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    order = sorted(refs.values(), key=lambda r: (
        r.cls != "tcr_chain_alpha", r.cls != "tcr_chain_beta", r.chain, r.seq))
    ids = {}
    with (out / f"{args.slug}_binders.csv").open("w", newline="") as fh, \
            (out / f"{args.slug}_binders.fasta").open("w") as fa:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for i, r in enumerate(order, start=1):
            rid = ids[id(r)] = f"{args.prefix}B{i:05d}"
            source = next((s for s in ("pdb", "iedb", "vdjdb", "patent") if s in r.sources), "")
            accession = r.accessions[0] if r.accessions else ""
            w.writerow({
                "id": rid, "class": r.cls, "binder_class_raw": r.cls,
                "binder_status": r.status, "named_agent": "; ".join(r.names),
                "length": len(r.seq), "chain_type": r.chain, "source_db": source,
                "accession": accession, "note": " | ".join(r.notes), "sequence": r.seq,
                "citation": "; ".join(r.citations), "evidence": "; ".join(r.evidence),
                "target_species": "human", "all_sources": "; ".join(r.sources),
                "all_accessions": "; ".join(r.accessions), "corpus_occurrences": "",
                "receptor_ids": ";".join(r.receptors),
            })
            # No whitespace: MMseqs2 keys a target on the header's first token, and the
            # agent travels in the header.
            fa.write(f">{rid}|{r.cls}|{r.status}|{source}|{accession.replace(' ', '_')}|"
                     f"{';'.join(n.replace(' ', '_') for n in r.names)}\n{r.seq}\n")

    # The CDR index: every CDR3 reference, plus CDR1/2/3 of each V domain, merged by
    # (sequence, region, chain) with the strongest status among their parents.
    cdrs: dict = {}
    for r in order:
        items = [("cdr3", r.seq)] if r.cls == "tcr_cdr3" else list(r.cdrs.items())
        for region, seq in items:
            parents = cdrs.setdefault((seq, region, r.chain), [])
            if r not in parents:
                parents.append(r)
    with (out / f"{args.slug}_binder_cdrs.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CDR_FIELDS)
        w.writeheader()
        for i, ((seq, region, chain), parents) in enumerate(
            sorted(cdrs.items(), key=lambda kv: (kv[0][1], kv[0][2], kv[0][0])), start=1
        ):
            names = []
            for p in parents:
                names.extend(n for n in p.names if n not in names)
            cites = [c for p in parents for c in p.citations]
            w.writerow({
                "cdr_id": f"{args.prefix}CDR{i:05d}", "cdr_seq": seq,
                "cdr_region": region.upper(), "chain_type": chain, "length": len(seq),
                "n_parent_sequences": len(parents),
                "best_status": max((p.status for p in parents), key=STATUS_RANK.get),
                "named_agents": "; ".join(names),
                "parent_classes": ";".join(sorted({p.cls for p in parents})),
                "source_dbs": "; ".join(sorted({s for p in parents for s in p.sources})),
                "parent_ids": "; ".join(ids[id(p)] for p in parents),
                "example_citation": cites[0] if cites else "",
                "corpus_occurrences": "",
                "receptor_ids": ";".join(
                    sorted({x for p in parents for x in p.receptors})
                ),
            })

    counts: dict = {}
    for r in order:
        counts[r.cls] = counts.get(r.cls, 0) + 1
    for s in skipped:
        print(f"  skipped {s}", file=sys.stderr)
    print(f"{kept} binders -> {len(order)} references {counts}, "
          f"{len(cdrs)} CDR index rows -> {out}")


if __name__ == "__main__":
    main()
