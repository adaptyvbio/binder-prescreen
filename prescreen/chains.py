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


@dataclass
class ChainIndex:
    """Both views of a batch, and the mapping between them."""

    #: ``{submission_id: "HHHH...LLLL"}`` — what classification and clustering see.
    concat: dict = field(default_factory=dict)
    #: ``{query_id: one_chain}`` — what the MMseqs2 arms see.
    queries: dict = field(default_factory=dict)
    #: ``{query_id: (submission_id, chain_number)}``, chain numbers 1-based.
    owner: dict = field(default_factory=dict)
    #: ``{submission_id: [query_id, ...]}`` in chain order.
    by_submission: dict = field(default_factory=dict)
    #: ``{submission_id: [chain_length, ...]}`` in chain order.
    lengths: dict = field(default_factory=dict)
    #: Submission ids whose sequence was empty once cleaned; dropped from the batch.
    dropped: list = field(default_factory=list)

    def chain_of(self, query_id: str) -> int:
        return self.owner[query_id][1]

    def submission_of(self, query_id: str) -> str:
        return self.owner[query_id][0]

    def map_rows(self) -> list:
        """Rows for ``query_map.tsv``: the synthetic id, and what it stands for."""
        return [
            (qid, sid, chain, len(self.queries[qid]))
            for qid, (sid, chain) in self.owner.items()
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
            qid = f"q{i}c{k}"
            index.queries[qid] = chain
            index.owner[qid] = (sid, k)
            index.by_submission[sid].append(qid)
    return index


def _relabel(row: dict, sid: str, chain: int) -> dict:
    """Stamp a hit row with the submission it belongs to and the chain that found it."""
    row["query"] = sid
    row["chain"] = chain
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
            row = _relabel(dict(row), sid, index.chain_of(qid))
            # Strict ``>`` while walking chains in order, so a tie goes to the lowest
            # chain number — the heavy chain under the template's convention.
            if chosen is None or row["similarity_check"] > chosen["similarity_check"]:
                chosen = row
        if chosen is not None:
            best[sid] = chosen
    return best
