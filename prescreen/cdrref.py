"""The paratope reference: known anti-target CDRs, framework-independent.

The curation track's measurement is the reason this file exists. The known-binder corpus
cross-hits itself at the framework level — a single probe VHH returns 1,843 of 4,964
entries at >= 0.8 query coverage, nearly all of it Ig-framework background. So a
whole-sequence search against that corpus cannot answer "is this a known binder": it
would flag almost every antibody-format submission. The binding region has to be
compared on its own, against regions rather than against whole sequences.

The reference is a deduplicated CDR table (``tnfa_binder_cdrs.csv``): 2,506 unique IMGT
CDRs extracted from the antibody-class entries, each carrying its chain type, how many
parent sequences it came from, the strongest binder status among those parents, and a
named agent where one is known.

Two details that change results:

* **Ambiguous residues.** Patent sequence listings use Markush-style claims, so some
  reference CDRs contain ``X``. Length-normalised identity against an ``X``-containing
  reference is meaningless, so those are excluded from identity arithmetic and matched
  only by a wildcard-aware exact comparison.
* **Chain type.** A heavy CDR3 is compared against heavy CDR3s. antpack types some
  single domains as light-like, so ``K`` and ``L`` are one class ("L").
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .compat import proteintyper

CANONICAL_CHAIN = {"H": "H", "K": "L", "L": "L"}
REGION_LABELS = ("cdr1", "cdr2", "cdr3")


@dataclass
class CdrEntry:
    cdr_id: str
    seq: str
    region: str
    chain: str
    status: str
    agents: str
    n_parents: int
    citation: str
    ambiguous: bool


@dataclass
class CdrMatch:
    """Best match of one query region against the reference CDRs."""

    region: str
    query: str
    identity: float
    exact: bool
    reference_id: str | None = None
    reference_seq: str | None = None
    status: str | None = None
    agents: str | None = None
    citation: str | None = None

    def as_dict(self) -> dict:
        return {
            "region": self.region,
            "query_region": self.query,
            "identity": round(self.identity, 4),
            "exact": self.exact,
            "reference_id": self.reference_id,
            "reference_region": self.reference_seq,
            "binder_status": self.status,
            "named_agent": self.agents,
            "citation": self.citation,
        }


def _wildcard_equal(query: str, reference: str) -> bool:
    """Exact match treating ``X`` in the reference as a wildcard."""
    if len(query) != len(reference):
        return False
    return all(r == "X" or q == r for q, r in zip(query, reference))


class CdrReference:
    """Known anti-target CDRs, indexed by (region, canonical chain)."""

    def __init__(self, csv_path: str | Path):
        self.path = Path(csv_path)
        self.entries: list = []
        self.index: dict = {}
        #: Curated binder classes that contributed at least one CDR. A class absent here
        #: cannot be corroborated by this index however good the whole-sequence match is,
        #: and a caller must not read that silence as "the paratope is new".
        self.covered_classes: set = set()
        with self.path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                seq = (row.get("cdr_seq") or "").strip().upper()
                if not seq:
                    continue
                entry = CdrEntry(
                    cdr_id=row.get("cdr_id") or "",
                    seq=seq,
                    region=(row.get("cdr_region") or "").lower(),
                    chain=CANONICAL_CHAIN.get((row.get("chain_type") or "").upper(), "H"),
                    status=row.get("best_status") or "",
                    agents=row.get("named_agents") or "",
                    n_parents=int(float(row.get("n_parent_sequences") or 1)),
                    citation=row.get("example_citation") or "",
                    ambiguous="X" in seq,
                )
                self.entries.append(entry)
                self.index.setdefault((entry.region, entry.chain), []).append(entry)
                for parent in (row.get("parent_classes") or "").split(";"):
                    if parent.strip():
                        self.covered_classes.add(parent.strip())

    def __len__(self) -> int:
        return len(self.entries)

    def counts(self) -> dict:
        return {f"{r}_{c}": len(v) for (r, c), v in sorted(self.index.items())}

    def match_region(
        self, query_seq: str, region: str, chain: str, min_length: int = 5
    ) -> CdrMatch:
        """Best reference match for one query region.

        Exact matches win outright (identity 1.0, ``exact=True``), including a
        wildcard-aware exact match against an ``X``-containing patent reference. Otherwise
        the best length-normalised Levenshtein identity is returned, computed only over
        unambiguous references.
        """
        cdr_novelty = proteintyper("cdr_novelty")
        query_seq = (query_seq or "").strip().upper()
        best = CdrMatch(region=region, query=query_seq, identity=0.0, exact=False)
        if len(query_seq) < min_length:
            return best
        chain = CANONICAL_CHAIN.get(chain or "H", "H")
        for entry in self.index.get((region, chain), ()):
            if entry.seq == query_seq or (
                entry.ambiguous and _wildcard_equal(query_seq, entry.seq)
            ):
                return CdrMatch(
                    region=region,
                    query=query_seq,
                    identity=1.0,
                    exact=True,
                    reference_id=entry.cdr_id,
                    reference_seq=entry.seq,
                    status=entry.status,
                    agents=entry.agents or None,
                    citation=entry.citation or None,
                )
            if entry.ambiguous:
                continue
            value = float(cdr_novelty.cdr_identity(query_seq, entry.seq))
            if value > best.identity:
                best = CdrMatch(
                    region=region,
                    query=query_seq,
                    identity=value,
                    exact=False,
                    reference_id=entry.cdr_id,
                    reference_seq=entry.seq,
                    status=entry.status,
                    agents=entry.agents or None,
                    citation=entry.citation or None,
                )
        return best

    def match_chains(self, regions: dict, min_length: int = 5) -> dict:
        """Match every region of every chain. Returns ``{chain: {region: CdrMatch}}``."""
        out: dict = {}
        for chain, labels in regions.items():
            canonical = CANONICAL_CHAIN.get(chain, chain)
            out[chain] = {}
            for label in REGION_LABELS:
                seq = labels.get(label)
                if not seq:
                    continue
                out[chain][label] = self.match_region(
                    seq, label, canonical, min_length=min_length
                )
        return out
