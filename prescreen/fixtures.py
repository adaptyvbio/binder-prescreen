"""Construction of labelled test sequences for calibrating the cut points.

Three kinds of synthetic submission are needed to see whether the filter separates what
it should:

``graft``
    A known binder's focus region (CDR3 or projected paratope) spliced into a *different*
    framework of the same category. Whole-sequence identity to the known binder stays
    low, region identity is 1.0. This is the case the region arm exists for, and the case
    a whole-sequence prescreen passes.
``framework_reuse``
    The reverse chimera: a known binder's framework carrying an unrelated region. Whole-
    sequence identity is high, region identity is low. These must *not* be flagged —
    reusing a published scaffold is allowed.
``mutant``
    A known binder with ``k`` random substitutions, giving a graded identity series for
    locating the cut points.

All construction is positional splicing of real sequences; nothing here invents sequence.
"""

from __future__ import annotations

import random

from .compat import proteintyper
from .regions import FOCUS, extract, focus_region

AA = "ACDEFGHIKLMNPQRSTVWY"


def focus_span(seq: str, category: str, region_kind: str) -> tuple[int, int] | None:
    """1-based inclusive span of the focus region in ``seq``, or ``None``.

    Antibody spans come from antpack numbering; scaffold and whole-sequence spans are
    located by substring search on the extracted region, which is exact because the
    extractors slice the query itself.
    """
    if region_kind == "antibody_cdrs":
        cdr_novelty = proteintyper("cdr_novelty")
        positions = cdr_novelty.cdr_residue_positions(seq, molecule_type="auto")
        for chain in ("H", "L"):
            pos = (positions.get(chain) or {}).get(FOCUS) or []
            if pos:
                return min(pos), max(pos)
        return None
    regions = extract(seq, category, region_kind)
    _chain, focus = focus_region(regions)
    if not focus:
        return None
    idx = seq.find(focus)
    if idx < 0:
        return None
    return idx + 1, idx + len(focus)


def graft(host_seq: str, donor_seq: str, category: str, region_kind: str) -> str | None:
    """Splice the donor's focus region into the host's focus span."""
    host_span = focus_span(host_seq, category, region_kind)
    if host_span is None:
        return None
    _chain, donor_focus = focus_region(extract(donor_seq, category, region_kind))
    if not donor_focus:
        return None
    start, end = host_span
    return host_seq[: start - 1] + donor_focus + host_seq[end:]


def framework_reuse(
    known_seq: str, donor_region: str, category: str, region_kind: str
) -> str | None:
    """Keep the known binder's framework, replace its focus region with ``donor_region``."""
    span = focus_span(known_seq, category, region_kind)
    if span is None:
        return None
    start, end = span
    return known_seq[: start - 1] + donor_region + known_seq[end:]


def mutant(seq: str, n_subs: int, seed: int = 0, protect: tuple | None = None) -> str:
    """Return ``seq`` with ``n_subs`` random substitutions.

    Args:
        seq: starting sequence.
        n_subs: number of positions to substitute.
        seed: RNG seed, so a calibration series is reproducible.
        protect: 1-based positions to leave untouched (e.g. a focus span).
    """
    rng = random.Random(seed)
    chars = list(seq)
    blocked = set()
    if protect:
        start, end = protect
        blocked = set(range(start - 1, end))
    choices = [i for i in range(len(chars)) if i not in blocked]
    rng.shuffle(choices)
    for idx in choices[: max(0, n_subs)]:
        current = chars[idx]
        chars[idx] = rng.choice([a for a in AA if a != current])
    return "".join(chars)


def identity_series(
    seq: str,
    fractions=(0.02, 0.05, 0.10, 0.20, 0.30, 0.50),
    seed: int = 0,
    protect: tuple | None = None,
) -> dict:
    """Build ``{label: sequence}`` mutants at the requested substitution fractions."""
    out = {}
    for frac in fractions:
        n = max(1, int(round(frac * len(seq))))
        out[f"mut{int(frac * 100):02d}"] = mutant(seq, n, seed=seed, protect=protect)
    return out
