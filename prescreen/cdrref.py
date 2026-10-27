"""The paratope reference: known anti-target CDRs, framework-independent.

The curation track's measurement is the reason this file exists. The known-binder corpus
cross-hits itself at the framework level — a single probe VHH returns 1,843 of 4,964
entries at >= 0.8 query coverage, nearly all of it Ig-framework background. So a
whole-sequence search against that corpus cannot answer "is this a known binder": it
would flag almost every antibody-format submission. The binding region has to be
compared on its own, against regions rather than against whole sequences.

The reference is the target set's deduplicated CDR table (``<slug>_binder_cdrs.csv``;
4,075 unique IMGT CDRs for TNF-alpha): regions extracted from the antibody-class entries
plus the bare claimed loops, each carrying its chain type, how many
parent sequences it came from, the strongest binder status among those parents, and a
named agent where one is known.

Two details that change results:

* **Ambiguous residues.** Patent sequence listings use Markush-style claims, so some
  reference CDRs contain ``X``. Length-normalised identity against an ``X``-containing
  reference is meaningless, so those are excluded from identity arithmetic and matched
  only by a wildcard-aware exact comparison.
* **Chain type.** A heavy CDR3 is compared against heavy CDR3s. antpack types some
  single domains as light-like, so ``K`` and ``L`` are one class ("L"). TCR chains keep
  their own letters (``A``, ``B``, ``D``, ``G``), and a TCR's CDR1/CDR2 are compared
  against the same chain only. Every CDR3 - antibody or TCR, either chain - is compared
  against every TCR CDR3: a TCR CDR3 lifted into an antibody (a TCR-mimic) or onto the
  other TCR chain is still the known paratope. Measured on VDJdb against the KRAS-G12D/V
  set, this costs almost nothing: of 179k CDR3s of other specificities, 8 come within 0.80
  of a reference on the opposite chain - most of them beta-like ``ASS...`` CDR3s filed
  under alpha - against 4,348 on their own chain.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .compat import proteintyper

CANONICAL_CHAIN = {"H": "H", "K": "L", "L": "L", "A": "A", "B": "B", "D": "D", "G": "G"}
TCR_CHAINS = ("A", "B", "D", "G")
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
    #: How many sequences in a general antibody corpus carry this exact region. A CDR
    #: seen in hundreds of unrelated antibodies is germline-encoded, so matching it says
    #: the query uses a common framework, not that it carries a known paratope. 0 means
    #: never seen - the distinctive case - and is also the default when the index has not
    #: been annotated, which keeps an un-annotated index behaving as it always did.
    occurrences: int = 0
    #: The receptors (rows of the source table) this region comes from - set for TCR
    #: sets, where specificity belongs to an alpha/beta pair, not to either CDR3.
    receptors: frozenset = frozenset()


@dataclass
class CdrMatch:
    """Best match of one query region against the reference CDRs."""

    region: str
    query: str
    identity: float
    exact: bool
    reference_id: str | None = None
    reference_seq: str | None = None
    occurrences: int | None = None
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
            "corpus_occurrences": self.occurrences,
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
                raw_region = (row.get("cdr_region") or "").lower()
                raw_chain = (row.get("chain_type") or "").upper()
                entry = CdrEntry(
                    cdr_id=row.get("cdr_id") or "",
                    seq=seq,
                    region=raw_region,
                    chain="*" if raw_chain == "*" else CANONICAL_CHAIN.get(raw_chain, "H"),
                    status=row.get("best_status") or "",
                    agents=row.get("named_agents") or "",
                    n_parents=int(float(row.get("n_parent_sequences") or 1)),
                    citation=row.get("example_citation") or "",
                    ambiguous="X" in seq,
                    occurrences=int(float(row.get("corpus_occurrences") or 0)),
                    receptors=frozenset(
                        r for r in (row.get("receptor_ids") or "").split(";") if r
                    ),
                )
                self.entries.append(entry)
                self.index.setdefault((entry.region, entry.chain), []).append(entry)
                for parent in (row.get("parent_classes") or "").split(";"):
                    if parent.strip():
                        self.covered_classes.add(parent.strip())
        # Rarest-first, but a named or confirmed binder's region always ahead of an
        # anonymous one: when several references match exactly, the informative one should
        # be the hit that gets reported and judged.
        for bucket in self.index.values():
            bucket.sort(
                key=lambda e: (
                    not (e.status == "confirmed_binder" or e.agents),
                    e.occurrences,
                    e.cdr_id,
                )
            )

        self.by_id = {e.cdr_id: e for e in self.entries}
        # receptor -> {TCR chain: [CDR3 entries]}, for match_pair
        self.receptor_cdr3: dict = {}
        for e in self.entries:
            if e.region == "cdr3" and e.chain in TCR_CHAINS and not e.ambiguous:
                for r in e.receptors:
                    self.receptor_cdr3.setdefault(r, {}).setdefault(e.chain, []).append(e)

    def __len__(self) -> int:
        return len(self.entries)

    def match_pair(self, cdr3s: dict) -> dict | None:
        """Best match of a TCR query's two CDR3s against ONE reference receptor.

        ``cdr3s`` is ``{chain: cdr3}`` with two TCR chains (``A`` and ``B``, or ``G`` and
        ``D``). Each receptor holding a CDR3 on both chains is scored by the *weaker* of the
        two identities, in the query's orientation or with the chains swapped (a CDR3 moved
        to the other chain is still that receptor's); the best receptor is returned as
        ``{"identity", "receptor", "references", "exact"}``, or ``None`` when no receptor
        has both chains.
        """
        cdr_novelty = proteintyper("cdr_novelty")
        (c1, q1), (c2, q2) = sorted(cdr3s.items())

        def best(q, entries):
            top = (-1.0, None)
            for e in entries:
                v = 1.0 if e.seq == q else float(cdr_novelty.cdr_identity(q, e.seq))
                if v > top[0]:
                    top = (v, e)
            return top

        out = None
        for receptor, chains in sorted(self.receptor_cdr3.items()):
            if len(chains) < 2:
                continue
            for x, y in ((c1, c2), (c2, c1)):
                if x not in chains or y not in chains:
                    continue
                (i1, e1), (i2, e2) = best(q1, chains[x]), best(q2, chains[y])
                score = min(i1, i2)
                if out is None or score > out["identity"]:
                    out = {
                        "identity": round(score, 4),
                        "receptor": receptor,
                        "references": [e1.cdr_id, e2.cdr_id],
                        "exact": i1 == 1.0 and i2 == 1.0,
                    }
        return out

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
        # Wildcard entries are bare claimed loops with no parent chain, so they are
        # consulted for any region on any chain: a patent-claimed binding region is prior
        # art wherever it reappears, not only at the position it was first described in.
        # Sorted rarest-first at index build, so the exact-match early return below picks
        # the most discriminative reference that matches rather than whichever happened to
        # be curated first: an exact hit on a CDR unique to one patent is real evidence,
        # the same hit on a germline-encoded one is not.
        candidates = [
            *self.index.get((region, chain), ()),
            *self.index.get((region, "*"), ()),
            *self.index.get(("*", "*"), ()),
        ]
        if region == "cdr3":
            for tcr_chain in TCR_CHAINS:
                if tcr_chain != chain:
                    candidates.extend(self.index.get((region, tcr_chain), ()))
        for entry in candidates:
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
                    occurrences=entry.occurrences,
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
                    occurrences=entry.occurrences,
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
