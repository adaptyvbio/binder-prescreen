"""T-cell receptor variable domains: recognising them, and reading their CDRs.

A TCR V domain is an immunoglobulin-superfamily V domain with the same IMGT CDR layout as
an antibody, and antpack numbers it as one: a TCR beta chain with its constant domain
types as a paired ``fv`` at 0.55-0.66 germline identity, an alpha chain as a
``single_domain_antibody``, well above ``ANTIBODY_GERMLINE_FLOOR``. Its CDRs then come
back keyed ``H``/``L`` and are compared against antibody CDRs, which is the wrong index
for a pMHC-binding receptor. antpack carries TCR germlines too (chains ``A``, ``B``,
``D``, ``G``), but refuses to mix them with antibody chains in one annotator, so a
sequence is numbered both ways and is a TCR only when the TCR germlines describe it
better *and* well.

Measured on 32 PDB TCR chains (KRAS-G12D, NY-ESO-1, MART-1, Tax, flu and SARS-CoV-2
TCRs): TCR germline identity 0.71-1.00. On 1,500 antibody and scaffold references from the
TNF-alpha, TfR1 and GLP-1R sets none reached 0.5 against TCR germlines, and none scored
higher against them than against antibody germlines.

CDR3 here is IMGT 105-117, without the conserved C104 and F/W118 - the convention of
IEDB, VDJdb and the antibody CDR3s this package already extracts (``CASSLAPGATNEKLFF`` is
reported as ``ASSLAPGATNEKLF``).
"""

from __future__ import annotations

from functools import lru_cache

from .compat import antpack_annotators, proteintyper

TCR_CHAINS = ("A", "B", "D", "G")

#: Germline identity a numbering needs before the domain counts as a TCR. Real TCR chains
#: score >= 0.71, antibody and scaffold references < 0.5.
TCR_GERMLINE_FLOOR = 0.60

#: Shortest segment worth numbering for a further domain, as for antibodies
#: (``chains.MIN_DOMAIN_SEGMENT``).
MIN_DOMAIN_SEGMENT = 80


@lru_cache(maxsize=1)
def _annotators():
    SingleChainAnnotator, _ = antpack_annotators()
    return (
        SingleChainAnnotator(chains=list(TCR_CHAINS), scheme="imgt"),
        SingleChainAnnotator(chains=["H", "K", "L"], scheme="imgt"),
    )


def _number(seq: str) -> dict | None:
    """Number ``seq`` as one TCR V domain, or ``None`` if it is not one."""
    cdr_novelty = proteintyper("cdr_novelty")
    tcr, antibody = _annotators()
    numbering, pid, chain, err = tcr.analyze_seq(seq)
    if err or pid < TCR_GERMLINE_FLOOR:
        return None
    labels = tcr.assign_cdr_labels(numbering, chain)
    cdrs: dict = {"cdr1": [], "cdr2": [], "cdr3": []}
    for i, label in enumerate(labels):
        if label in cdrs:
            cdrs[label].append(i)
    valid, _reason = cdr_novelty.numbering_is_valid(
        seq, numbering, {k: len(v) for k, v in cdrs.items()}, pid
    )
    if not valid:
        return None
    span = [i for i, n in enumerate(numbering) if n != "-"]
    domain = seq[span[0] : span[-1] + 1]
    # The same domain against antibody germlines. A TCR is not an antibody it resembles;
    # this keeps any antibody that TCR germlines happen to fit from being retyped.
    _n, ab_pid, _c, _e = antibody.analyze_seq(domain)
    if ab_pid >= pid:
        return None
    return {
        "chain": chain,
        "start": span[0],
        "end": span[-1] + 1,
        "percent_identity": round(float(pid), 3),
        "antibody_identity": round(float(ab_pid), 3),
        **{k: "".join(seq[i] for i in v) for k, v in cdrs.items()},
        "positions": {k: [i + 1 for i in v] for k, v in cdrs.items()},
    }


def _scan(seq: str, lo: int, hi: int, out: list) -> None:
    """Collect TCR V domains in ``seq[lo:hi]``, then in what the first one left over."""
    if hi - lo < MIN_DOMAIN_SEGMENT:
        return
    found = _number(seq[lo:hi])
    if found is None:
        return
    found["start"] += lo
    found["end"] += lo
    found["positions"] = {
        k: [p + lo for p in v] for k, v in found["positions"].items()
    }
    out.append(found)
    _scan(seq, lo, found["start"], out)
    _scan(seq, found["end"], hi, out)


@lru_cache(maxsize=4096)
def _domains(seq: str) -> tuple:
    out: list = []
    _scan(seq, 0, len(seq), out)
    return tuple(sorted(out, key=lambda d: d["start"]))


def tcr_domains(seq: str) -> list:
    """The TCR V domains of ``seq``, N- to C-terminal; ``[]`` if it carries none.

    Each is ``{"chain", "start", "end", "percent_identity", "antibody_identity", "cdr1",
    "cdr2", "cdr3", "positions"}``, with ``start``/``end`` a 0-based slice and
    ``positions`` 1-based residue numbers into ``seq``. A constant domain does not number
    as a V domain, so an alpha chain with its C-alpha is one domain, and a single-chain
    Valpha-linker-Vbeta two.
    """
    seq = seq.strip().upper().replace("*", "")
    return [dict(d) for d in _domains(seq)]


def tcr_cdrs(seq: str) -> dict:
    """``{chain: {"cdr1", "cdr2", "cdr3"}}`` for the TCR domains of ``seq``.

    Keyed by antpack's chain letter (``A``/``B``/``D``/``G``); where two domains share a
    chain type the N-terminal one is kept, as the antibody extractor keeps one H and one L.
    """
    out: dict = {}
    for d in tcr_domains(seq):
        if d["chain"] not in out:
            out[d["chain"]] = {k: d[k] for k in ("cdr1", "cdr2", "cdr3") if d[k]}
    return out


def tcr_cdr_positions(seq: str) -> dict:
    """``{chain: {"cdr1": [...], ...}}``, 1-based positions, keyed like :func:`tcr_cdrs`."""
    out: dict = {}
    for d in tcr_domains(seq):
        out.setdefault(d["chain"], d["positions"])
    return out
