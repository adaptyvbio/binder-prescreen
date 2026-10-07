"""Splitting a multi-chain submission, and merging its per-chain search results back.

The Proteinbase submission template joins the chains of a multi-chain entry with ``:``
(a Fab is ``heavy:light``). Those chains are handled two different ways on purpose:

**Classification reads the concatenation.** ``classify()`` recognises a paired H+L inside
one string and reads the subformat off the constant domains present, which is how a Fab
types as ``fab`` rather than as two loose chains. Splitting first would lose that.

**The searches run per chain.** Every reference in every arm is a single chain — the
packaged known-binder set has a median length of 108 aa and only 1.3% of entries over
400 aa. Searching the concatenation dilutes query coverage by exactly the fraction of the
molecule the matching chain is not: a verbatim adalimumab Fab queried as 432 aa against its
own 219-aa heavy chain scores ``fident = 1.00`` but ``qcov = 0.51``, so the composite is
0.51 — under ``near_known_sequence`` (0.80), ``known_sequence`` (0.95) and
``tnf_known_coverage`` (0.80) alike. Those flags are unreachable by construction for any
Fab or IgG submitted in template form. Per chain the same query scores 1.00.

Query ids for the search arms are synthetic (``q0c1``, ``q0c2``, ...) rather than derived
from the submitter's name. MMseqs2 keys a FASTA accession on the first whitespace token, as
do both readers in this package, so a submission named ``my fab v2`` comes back from the
aligner as query ``my``, matches no key in the results dict, and has its hits dropped in
silence. An opaque id cannot collide with anything the submitter typed. Rows are relabelled
with the real submission id before they leave :func:`prescreen.screen`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .compat import proteintyper

#: A V domain is 110-130 residues. A leftover segment shorter than this cannot hold one, so
#: the search for further domains stops rather than numbering noise.
MIN_DOMAIN_SEGMENT = 80


def _germline_identity(seq: str) -> float:
    """Best germline percent identity over the chains antpack numbers in ``seq``."""
    cdr_novelty = proteintyper("cdr_novelty")
    try:
        chains = cdr_novelty.extract_cdrs(seq, molecule_type="auto")
    except Exception:
        return 0.0
    pids = [float(c.get("percent_identity") or 0.0) for c in chains.values()]
    return max(pids) if pids else 0.0


def _spans_in(seq: str, lo: int, hi: int, out: list) -> None:
    """Collect V-domain spans in ``seq[lo:hi]``, then recurse into what they did not cover.

    ``variable_domain_spans`` returns a dict keyed by canonical chain, so it holds at most
    one span per chain type and cannot report more than two domains however many are
    present: on a trivalent VHH it returns two and drops the third. Re-running it on the
    uncovered segments recovers the rest without adding any numbering code of our own.
    """
    if hi - lo < MIN_DOMAIN_SEGMENT:
        return
    cdr_novelty = proteintyper("cdr_novelty")
    try:
        spans = cdr_novelty.variable_domain_spans(seq[lo:hi])
    except Exception:
        return
    ordered = sorted(spans.values())
    if not ordered:
        return
    # Same guard the library applies: overlapping spans mean the numbering placed two
    # domains on the same residues, so this segment cannot be trusted.
    for (_, prev_end), (next_start, _) in zip(ordered, ordered[1:]):
        if next_start <= prev_end:
            return
    bounds = []
    for start, end in ordered:
        a, b = lo + start - 1, lo + end
        out.append((a, b))
        bounds.append((a, b))
    # The gaps before, between and after what was just found.
    edges = [lo] + [x for ab in bounds for x in ab] + [hi]
    for g_lo, g_hi in zip(edges[0::2], edges[1::2]):
        if g_hi - g_lo >= MIN_DOMAIN_SEGMENT:
            _spans_in(seq, g_lo, g_hi, out)


def variable_domains(seq: str, floor: float | None = None) -> list:
    """The variable domains of ``seq``, N- to C-terminal, or ``[]`` if fewer than two.

    Returning ``[]`` for a single-domain input is the contract the callers rely on: a
    nanobody, a lone VH, an IgG heavy chain whose only real V domain is its VH, a peptide
    and a non-antibody are all left to be handled whole, exactly as before.

    Candidate domains are filtered by germline identity. antpack numbers a constant domain
    as readily as a variable one — 41 of 149 ``igg_chain_heavy`` references in the curated
    set split into a real VH plus a CH domain at ~0.35 identity — so a domain counts only
    if it clears the same floor the classifier uses to decide an antibody is an antibody.
    Genuine dual-variable heavy chains and bivalent VHHs sit at 0.92-0.98 and survive.
    """
    from .classify import ANTIBODY_GERMLINE_FLOOR

    floor = ANTIBODY_GERMLINE_FLOOR if floor is None else floor
    found: list = []
    _spans_in(seq, 0, len(seq), found)
    if len(found) < 2:
        return []
    kept = [seq[a:b] for a, b in sorted(found) if _germline_identity(seq[a:b]) >= floor]
    return kept if len(kept) >= 2 else []


@dataclass
class ChainIndex:
    """Both views of a batch, and the mapping between them."""

    #: ``{submission_id: "HHHH...LLLL"}`` — what classification and clustering see.
    concat: dict = field(default_factory=dict)
    #: ``{query_id: one_chain}`` — what the MMseqs2 arms see.
    queries: dict = field(default_factory=dict)
    #: ``{query_id: (submission_id, chain_number, domain_number)}``, both 1-based. The
    #: domain number is 1 for a chain that was searched whole.
    owner: dict = field(default_factory=dict)
    #: ``{submission_id: [query_id, ...]}`` in chain order.
    by_submission: dict = field(default_factory=dict)
    #: ``{submission_id: [chain_length, ...]}`` in chain order.
    lengths: dict = field(default_factory=dict)
    #: Submission ids whose sequence was empty once cleaned; dropped from the batch.
    dropped: list = field(default_factory=list)

    def chain_of(self, query_id: str) -> int:
        return self.owner[query_id][1]

    def domain_of(self, query_id: str) -> int:
        return self.owner[query_id][2]

    def submission_of(self, query_id: str) -> str:
        return self.owner[query_id][0]

    def map_rows(self) -> list:
        """Rows for ``query_map.tsv``: the synthetic id, and what it stands for."""
        return [
            (qid, sid, chain, domain, len(self.queries[qid]))
            for qid, (sid, chain, domain) in self.owner.items()
        ]


def split_submissions(records: dict) -> ChainIndex:
    """Build a :class:`ChainIndex` from ``{submission_id: raw_sequence}``.

    A single-chain submission comes out byte-identical to the previous normalisation
    (``strip`` → upper-case → drop ``*``), so nothing downstream of the concatenation
    changes for it.
    """
    index = ChainIndex()
    for i, (sid, raw) in enumerate(records.items()):
        if not raw:
            index.dropped.append(sid)
            continue
        parts = [c.strip().upper().replace("*", "") for c in str(raw).split(":")]
        parts = [c for c in parts if c]
        if not parts:
            # ":", "   " and "***" all clean to nothing. Letting one through reaches
            # classify(""), which raises and takes the whole batch down with it.
            index.dropped.append(sid)
            continue
        index.concat[sid] = "".join(parts)
        index.by_submission[sid] = []
        index.lengths[sid] = [len(c) for c in parts]
        for k, chain in enumerate(parts, start=1):
            # A chain carrying two or more variable domains is searched and compared as
            # those domains; one carrying fewer is left exactly as it is. The linker and
            # any tag between domains belong to neither and are dropped, which is what
            # stops them diluting whole-sequence coverage against single-chain references.
            units = variable_domains(chain) or [chain]
            for d, unit in enumerate(units, start=1):
                qid = f"q{i}c{k}d{d}"
                index.queries[qid] = unit
                index.owner[qid] = (sid, k, d)
                index.by_submission[sid].append(qid)
    return index


def _relabel(row: dict, sid: str, chain: int, domain: int = 1) -> dict:
    """Stamp a hit row with the submission, chain and domain that produced it."""
    row["query"] = sid
    row["chain"] = chain
    row["domain"] = domain
    return row


def merge_whole(per_chain: dict, index: ChainIndex) -> dict:
    """Collapse ``search_whole``'s ``{query_id: row}`` to one row per submission.

    The best row by composite wins outright; its ``fident`` and ``qcov`` are the ones
    ``flags.evaluate`` will read. They are never recombined across chains — taking the best
    identity from one chain and the best coverage from another would manufacture a
    known-binder verdict out of evidence that exists in no single alignment.
    """
    best: dict = {}
    for sid, qids in index.by_submission.items():
        chosen = None
        for qid in qids:
            row = per_chain.get(qid)
            if row is None:
                continue
            row = _relabel(dict(row), sid, index.chain_of(qid), index.domain_of(qid))
            # Strict ``>`` while walking chains in order, so a tie goes to the lowest
            # chain number — the heavy chain under the template's convention.
            if chosen is None or row["similarity_check"] > chosen["similarity_check"]:
                chosen = row
        if chosen is not None:
            best[sid] = chosen
    return best
