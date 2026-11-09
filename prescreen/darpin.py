"""DARPin recognition and paratope projection, segment by segment.

A DARPin is an N-cap, two or more 33-residue internal repeats and a C-cap. The library
check in ``scaffold_cdr.framework_identity`` aligns the whole query to one reference, the
three-repeat NI3C consensus with the original C-cap, and gates at 0.85. That penalises
two things that say nothing about whether a sequence is a DARPin:

* the repeat count - a two-repeat design is scored against a three-repeat reference;
* the cap variant - the Mii C-cap (``...KTAFDISIDNGNEDLAEILQKLN``, Interlandi et al. 2008)
  that H10-2-G3 (PDB 2JAB) and many later DARPins carry costs ~7 framework mismatches.

Generative redesigns of a published DARPin also touch framework positions outside the
seven randomised ones. Together these put real DARPins under the gate: two-repeat
redesigns of H10-2-G3 (PDB 2JAB) scored 0.79-0.83 and were screened as plain
miniproteins, so the whole-sequence match to their own template decided the verdict.

Here each segment is scored on its own - the N-cap and C-cap by local identity to the
known cap variants, each internal repeat by identity to the consensus repeat at its 26
non-randomised positions, found by a sliding scan rather than the ``G[AS]DVN[AV]`` anchor
regex - and combined into one length-weighted score.

Calibration (``scripts/calibrate_darpin.py``, PDB seqres and UniProt, 2026-11-09): 158
unique DARPin chains named as such in the PDB, including fusions, scored >= 0.827, while
1,376 natural ankyrin-repeat proteins (1,266 reviewed UniProt entries with InterPro
IPR002110, plus 110 PDB chains) scored <= 0.782. Only GABP-beta1 comes close (0.78): the
DARPin N-cap was derived from GABP-beta, so its caps score 0.87/0.82, but its repeats
stay at 0.62. The second-repeat floor exists for that case. The single-reference gate
accepted 156 of the 158 DARPins and this rule accepts all 158. It accepts no natural
ankyrin and none of 9,950 affibody or 4,115 monobody sequences. With random substitutions
spread over the whole sequence, NI3C stays recognised up to ~20 % (the old gate held to
~10 %) and 2JAB up to ~15 % (old: < 5 %).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

#: Consensus internal repeat (NI3C), framed as in ``scaffold_cdr``: starting one residue
#: after the ``GADVNA`` motif that closes the previous repeat.
REPEAT = "DKDGYTPLHLAAREGHLEIVEVLLKAGADVNAK"
REPEAT_LEN = len(REPEAT)
#: Randomised positions in a repeat, 1-based - ``scaffold_cdr.DARPIN_MASK``, restated so
#: this module does not depend on the vendored copy for its numbering.
MASK = (2, 3, 5, 13, 14, 26, 33)
_FRAMEWORK = np.array([p - 1 for p in range(1, REPEAT_LEN + 1) if p not in MASK])
_REPEAT_FW = np.frombuffer(REPEAT.encode(), dtype=np.uint8)[_FRAMEWORK]

N_CAP = "DLGKKLLEAARAGQDDEVRILMANGADVNA"
#: The original consensus C-cap and the stabilised Mii C-cap.
C_CAPS = ("QDKFGKTPFDLAIDNGNEDIAEVLQKAA", "QDKFGKTAFDISIDNGNEDLAEILQKLN")

#: Minimum framework identity for a window to count as an internal repeat when framing
#: the paratope. The caps share the ``GADVNA``/``LLKxG`` motifs and score 0.46-0.50 on
#: this scan; internal repeats of designed DARPins scored >= 0.69.
REPEAT_FRAME_MIN = 0.60
#: Minimum length-weighted segment score for a DARPin call (see the module docstring).
SEGMENT_SCORE_MIN = 0.80
#: Minimum framework identity of the second-best repeat: a DARPin has at least two
#: internal repeats. The lowest designed DARPin sits at 0.69 and GABP-beta1 at 0.62.
SECOND_REPEAT_MIN = 0.65
#: Repeat windows closer than this are the same repeat seen at a shifted frame.
_MIN_REPEAT_SPACING = 30
#: Scan floor for reporting a window at all; below it a window is background.
_SCAN_FLOOR = 0.45


@lru_cache(maxsize=1)
def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner(mode="local")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    return aligner


def _cap_identity(cap: str, seq: str) -> float:
    """Identical aligned residues of the best local alignment, over the cap length."""
    aln = _aligner().align(cap, seq)[0]
    matches = sum(
        cap[s1 + k] == seq[s2 + k]
        for (s1, e1), (s2, _e2) in zip(*aln.aligned)
        for k in range(e1 - s1)
    )
    return matches / len(cap)


def repeat_windows(seq: str, floor: float = _SCAN_FLOOR) -> list[tuple[int, float]]:
    """Non-overlapping 33-residue windows scored against the consensus repeat framework.

    Returns ``(offset, identity)`` pairs, 0-based, N- to C-terminal. Windows are taken
    best-first and a window within ``_MIN_REPEAT_SPACING`` of a better one is dropped, so
    each repeat is reported once, at its best frame.
    """
    b = np.frombuffer(seq.encode(), dtype=np.uint8)
    n = len(b) - REPEAT_LEN + 1
    if n <= 0:
        return []
    scores = (b[np.arange(n)[:, None] + _FRAMEWORK[None, :]] == _REPEAT_FW).mean(axis=1)
    picked: list[tuple[int, float]] = []
    for i in np.argsort(-scores, kind="stable"):
        if scores[i] < floor:
            break
        if all(abs(int(i) - j) >= _MIN_REPEAT_SPACING for j, _ in picked):
            picked.append((int(i), float(scores[i])))
    return sorted(picked)


def evidence(seq: str) -> dict:
    """Per-segment DARPin evidence and the combined score."""
    seq = seq.strip().upper().replace("*", "")
    windows = repeat_windows(seq)
    top = sorted((s for _, s in windows), reverse=True)[:2]
    top += [0.0] * (2 - len(top))
    n_cap = _cap_identity(N_CAP, seq) if seq else 0.0
    c_cap = max(_cap_identity(c, seq) for c in C_CAPS) if seq else 0.0
    weight = len(N_CAP) + len(C_CAPS[0]) + 2 * len(_FRAMEWORK)
    score = (
        len(N_CAP) * n_cap + len(C_CAPS[0]) * c_cap + len(_FRAMEWORK) * sum(top)
    ) / weight
    return {
        "n_cap_identity": round(n_cap, 4),
        "c_cap_identity": round(c_cap, 4),
        "repeat_identities": [round(s, 4) for _, s in windows if s >= REPEAT_FRAME_MIN],
        "second_repeat_identity": round(top[1], 4),
        "segment_score": round(score, 4),
        "segment_score_min": SEGMENT_SCORE_MIN,
        "second_repeat_min": SECOND_REPEAT_MIN,
    }


def is_darpin(ev: dict) -> bool:
    """Whether :func:`evidence` describes a DARPin."""
    return (
        ev["segment_score"] >= SEGMENT_SCORE_MIN
        and ev["second_repeat_identity"] >= SECOND_REPEAT_MIN
    )


def repeat_offsets(seq: str) -> list[int]:
    """0-based offsets of the complete internal repeats, N- to C-terminal.

    Replaces ``scaffold_cdr.darpin_repeat_offsets`` for the prescreen. That function
    frames a repeat after every ``G[AS]DVN[AV]`` anchor, so it drops a repeat whose
    closing motif was redesigned (``NGMDVNA``, ``GADVDA``) and frames the C-cap as a
    fourth repeat whenever a tag or linker follows it.
    """
    return [
        off
        for off, s in repeat_windows(seq, floor=REPEAT_FRAME_MIN)
        if off + REPEAT_LEN <= len(seq)
    ]


def project(seq: str) -> dict:
    """DARPin paratope regions, the shape ``regions.extract`` returns for a scaffold.

    The seven masked positions of each repeat, concatenated N- to C-terminal, are the
    focus region (``cdr3``). ``cdr1``/``cdr2`` carry the first two repeats. This is the
    convention of ``scaffold_cdr.project_darpin``, so
    ``scaffold_search.darpin_repeat_identity`` still applies.
    """
    seq = seq.strip().upper().replace("*", "")
    per_repeat = [
        "".join(seq[off + p - 1] for p in MASK) for off in repeat_offsets(seq)
    ]
    if not per_repeat:
        return {}
    return {
        "S": {
            "cdr1": per_repeat[0],
            "cdr2": per_repeat[1] if len(per_repeat) > 1 else "",
            "cdr3": "".join(per_repeat),
        }
    }
