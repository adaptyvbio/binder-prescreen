"""Scaffold-CDR region extraction for non-antibody binder scaffolds.

Antibody CDR novelty rests on two pieces that have no off-the-shelf
equivalent for alternative scaffolds:

  1. a residue numbering scheme (antpack/IMGT) that makes CDR boundaries
     reproducible across sequences of different length, and
  2. a curated paratope reference database (PLAbDab) to measure identity against.

This module supplies (1) for affibody, DARPin and monobody scaffolds using
MMseqs2-/alignment-based projection of a canonical randomization mask onto the
query, and emits regions in the same ``{chain: {cdr1, cdr2, cdr3}}`` shape that
``cdr_novelty.analyze_cdr_novelty`` already consumes. (2) is built by running the
same projection over family members retrieved from PDB seqres / SwissProt.

The masks are literature-canonical and were each re-derived empirically from
known members (see the feasibility report):
  affibody  13 surface positions on helices 1 and 2 of the 58-aa protein A Z domain
  DARPin     6 randomized paratope positions + 1 framework position per 33-aa repeat
  monobody   BC, DE and FG loops of the 94-aa fibronectin 10th FN3 domain
"""
# ----------------------------------------------------------------------------------
# Vendored verbatim from the upstream library (commit a02fae8), except where marked
# `vendored:`. Do not edit to fix a bug here - fix it upstream and re-copy, so both
# copies keep giving the same numbering and the same identities.
# ----------------------------------------------------------------------------------

from __future__ import annotations

import re

from Bio.Align import PairwiseAligner, substitution_matrices

# --------------------------------------------------------------------------------------
# Framework references and canonical randomization masks
# --------------------------------------------------------------------------------------

# Protein A Z domain (PDB 1Q2N chain A), 58 aa.
Z_DOMAIN = "VDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPK"
# Nord K, Gunneriusson E, Ringdahl J, Stahl S, Uhlen M, Nygren PA (1997)
# Nat Biotechnol 15(8):772-777 — 13 randomized surface positions in helices 1 and 2,
# the original Fc-binding face of the Z domain.
AFFIBODY_MASK = (9, 10, 11, 13, 14, 17, 18, 24, 25, 27, 28, 32, 35)
AFFIBODY_HELIX1 = (9, 10, 11, 13, 14, 17, 18)
AFFIBODY_HELIX2 = (24, 25, 27, 28, 32, 35)

# Human fibronectin 10th type III domain, FNfn10 (PDB 1TTG chain A), 94 aa.
FN3_DOMAIN = (
    "VSDVPRDLEVVAATPTSLLISWDAPAVTVRYYRITYGETGGNSPVQEFTVPGSKS"
    "TATISGLKPGVDYTITVYAVTGRGDSPASSKPISINYRT"
)
# Loop spans in FNfn10 numbering (Koide et al. 1998 J Mol Biol 284:1141-1151).
# Note these cover the paratope loops only: 'side-and-loop' libraries also diversify
# beta-sheet surface positions outside these spans, which this mask does not capture.
MONOBODY_LOOPS = {"BC": (21, 31), "DE": (51, 56), "FG": (75, 88)}

# Full-consensus DARPin NI3C (PDB 2XEE chain A) as the repeat frame reference.
DARPIN_REF = (
    "DLGKKLLEAARAGQDDEVRILMANGADVNAKDKDGYTPLHLAAREGHLEIVEVLLKAGADVNA"
    "KDKDGYTPLHLAAREGHLEIVEVLLKAGADVNAKDKDGYTPLHLAAREGHLEIVEVLLKAGADVNA"
    "QDKFGKTPFDLAIDNGNEDIAEVLQKAA"
)
# Binz HK et al. (2003) J Mol Biol 332:489-503; Wetzel SK et al. (2008) J Mol Biol
# 376:241-257. Within the 33-residue internal repeat, positions 2, 3, 5, 13, 14, 33 are
# the randomized paratope residues and position 26 is the randomized framework residue
# (z, restricted to N/H/Y). DARPIN_MASK is the union of the two, i.e. 7 positions;
# DARPIN_PARATOPE is the published 6-position paratope set on its own.
DARPIN_REPEAT_LEN = 33
DARPIN_MASK = (2, 3, 5, 13, 14, 26, 33)
DARPIN_PARATOPE = (2, 3, 5, 13, 14, 33)
# Conserved motif closing each internal repeat; the next repeat frame starts one
# residue after it.
DARPIN_ANCHOR = re.compile(r"G[AS]DVN[AV]")

SCAFFOLD_REFS = {
    "affibody": Z_DOMAIN,
    "monobody": FN3_DOMAIN,
    "darpin": DARPIN_REF,
}


def _aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    return a


_ALIGNER = _aligner()


def _position_map(ref: str, query: str) -> dict:
    """Map 1-based reference positions onto query residues via global alignment."""
    aln = _ALIGNER.align(ref, query)[0]
    pmap = {}
    for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(e1 - s1):
            pmap[s1 + k + 1] = query[s2 + k]
    return pmap


def _span_from_alignment(ref: str, query: str, start: int, end: int) -> str:
    """Return the query substring aligned to reference span [start, end] (1-based).

    Unlike a position map this preserves insertions, so a loop that is longer or
    shorter in the query than in the reference is still returned in full. This is
    what makes the monobody FG loop (length-variable between designs) usable.
    """
    aln = _ALIGNER.align(ref, query)[0]
    blocks = list(zip(aln.aligned[0], aln.aligned[1]))
    q_lo = q_hi = None
    for (s1, e1), (s2, e2) in blocks:
        if e1 <= start - 1 or s1 >= end:
            continue
        lo = max(s1, start - 1)
        hi = min(e1, end)
        q_lo = s2 + (lo - s1) if q_lo is None else q_lo
        q_hi = s2 + (hi - s1)
    if q_lo is None or q_hi is None:
        return ""
    return query[q_lo:q_hi]


def framework_identity(family: str, query: str) -> float:
    """Identity of the query to the family framework reference OUTSIDE the mask.

    This is the discriminator between a designed scaffold variant (framework
    essentially intact, mask region rewritten) and a natural family homologue
    (framework itself diverged).
    """
    ref = SCAFFOLD_REFS[family]
    pmap = _position_map(ref, query)
    if family == "affibody":
        mask = set(AFFIBODY_MASK)
    elif family == "monobody":
        mask = {p for s, e in MONOBODY_LOOPS.values() for p in range(s, e + 1)}
    else:
        mask = set()
        for off in darpin_repeat_offsets(ref):
            mask |= {off + p for p in DARPIN_MASK}
    fw = [(p, c) for p, c in pmap.items() if p not in mask]
    if not fw:
        return 0.0
    return sum(1 for p, c in fw if c == ref[p - 1]) / len(fw)


def darpin_repeat_offsets(seq: str) -> list:
    """Return 0-based offsets of each complete 33-residue internal repeat frame."""
    anchors = [m.start() for m in DARPIN_ANCHOR.finditer(seq)]
    offsets = []
    for a in anchors:
        start = a + 6 + 1  # past the anchor motif, past the residue closing the repeat
        if start + DARPIN_REPEAT_LEN <= len(seq):
            offsets.append(start)
    return offsets


def project_affibody(seq: str) -> dict:
    """Return affibody paratope regions keyed cdr1 (helix 1), cdr2 (helix 2), cdr3 (all 13)."""
    pmap = _position_map(Z_DOMAIN, seq)
    pick = lambda ps: "".join(pmap.get(p, "-") for p in ps)
    return {
        "cdr1": pick(AFFIBODY_HELIX1),
        "cdr2": pick(AFFIBODY_HELIX2),
        "cdr3": pick(AFFIBODY_MASK),
    }


def project_monobody(seq: str) -> dict:
    """Return monobody paratope regions keyed cdr1 (BC), cdr2 (DE), cdr3 (BC+DE+FG).

    ``cdr3`` carries the whole paratope, not just FG, because ``cdr3`` is the region the
    production scale scores — the same convention as the affibody (all 13 mask positions)
    and the DARPin (all repeats concatenated). Keying FG alone there left BC and DE
    scored by nothing: ``framework_identity`` already excludes all three loops from the
    0.85 gate, so a design with BC and DE fully rewritten and FG copied read as a known
    paratope. The 0.61 level-3/4 cut is derived from the full 31-of-94 paratope share
    (see ``novelty_scale.md``), so it matches this definition and not an FG-only one.
    """
    spans = {
        name: _span_from_alignment(FN3_DOMAIN, seq, *MONOBODY_LOOPS[name])
        for name in ("BC", "DE", "FG")
    }
    # A loop that could not be anchored comes back "". Concatenating whatever survived
    # would silently shorten the paratope, and identity normalizes by the longer string:
    # a query missing BC caps at 20/31 = 0.645 against a full-length member, which the
    # 0.70 cut reads as de novo even when every residue present is byte-identical. Report
    # no paratope at all instead, so the emptiness checks the callers already apply
    # (``paratope._has_regions``, ``scaffold_search``'s focus guard) route this to the
    # structural fallback or to "undetermined" rather than to a fabricated verdict.
    if not all(spans.values()):
        return {"cdr1": spans["BC"], "cdr2": spans["DE"], "cdr3": ""}
    return {
        "cdr1": spans["BC"],
        "cdr2": spans["DE"],
        "cdr3": spans["BC"] + spans["DE"] + spans["FG"],
    }


def project_darpin(seq: str) -> dict:
    """Return DARPin paratope regions.

    Repeat count varies between designs, so the mask is projected per detected
    internal repeat and the per-repeat residues are concatenated N- to C-terminal.
    cdr1/cdr2 carry the first two repeats individually; cdr3 carries the full
    concatenated paratope, which is the CDR-H3 analogue the production scale reads.
    """
    offsets = darpin_repeat_offsets(seq)
    per_repeat = []
    for off in offsets:
        rep = seq[off : off + DARPIN_REPEAT_LEN]
        if len(rep) < DARPIN_REPEAT_LEN:
            continue
        per_repeat.append("".join(rep[p - 1] for p in DARPIN_MASK))
    return {
        "cdr1": per_repeat[0] if len(per_repeat) > 0 else "",
        "cdr2": per_repeat[1] if len(per_repeat) > 1 else "",
        "cdr3": "".join(per_repeat),
        "n_repeats": len(per_repeat),
    }


_PROJECTORS = {
    "affibody": project_affibody,
    "monobody": project_monobody,
    "darpin": project_darpin,
}


def extract_scaffold_regions(seq: str, family: str) -> dict:
    """Scaffold analogue of ``cdr_novelty.extract_cdrs``.

    Returns ``{"S": {"cdr1", "cdr2", "cdr3"}}`` — a single pseudo-chain "S", since
    these scaffolds are single-domain and have no heavy/light pairing. The shape
    matches what ``analyze_cdr_novelty`` consumes so the comparison and scoring
    stages are reused unchanged.
    """
    if family not in _PROJECTORS:
        raise ValueError(f"no projector for scaffold family {family!r}")
    regions = _PROJECTORS[family](seq)
    regions.pop("n_repeats", None)
    return {"S": regions}
