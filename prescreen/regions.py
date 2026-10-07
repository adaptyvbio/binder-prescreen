"""Binding-region extraction, one function per category family.

Every extractor returns the same shape that ``cdr_novelty.compare_regions`` consumes —
``{chain: {"cdr1": str, "cdr2": str, "cdr3": str}}`` — so a single comparison core works
for antibodies, alternative scaffolds and plain sequences alike. The ``cdr3`` slot is the
*focus* region in all three cases:

==================  =========================================  ======================
region_kind         cdr1 / cdr2 / cdr3                         focus (``cdr3``)
==================  =========================================  ======================
antibody_cdrs       CDR1 / CDR2 / CDR3 per chain (IMGT)        CDR3
scaffold_paratope   per-element paratope spans, cdr3 = all     whole paratope
whole               cdr3 only                                  the entire sequence
==================  =========================================  ======================

The whole-sequence case deliberately reuses the ``cdr3`` slot rather than inventing a
fourth code path: for a peptide or an unannotatable miniprotein the binding region
*is* the sequence, and routing it through the same slot means the flag logic reads one
number regardless of category.
"""

from __future__ import annotations

from .compat import proteintyper

FOCUS = "cdr3"


def extract(seq: str, category: str, region_kind: str) -> dict:
    """Extract the comparison regions for one sequence.

    Args:
        seq: amino-acid sequence.
        category: category slug from :func:`prescreen.classify.classify`.
        region_kind: ``antibody_cdrs``, ``scaffold_paratope`` or ``whole``.

    Returns:
        ``{chain: {region_label: sequence}}``; empty dict if extraction failed, which the
        caller must treat as "region arm not applicable" rather than as novelty.
    """
    seq = seq.strip().upper().replace("*", "")
    if region_kind == "antibody_cdrs":
        cdr_novelty = proteintyper("cdr_novelty")
        try:
            chains = cdr_novelty.extract_cdrs(seq, molecule_type="auto")
        except Exception:
            return {}
        return {
            chain: {k: v for k, v in regions.items() if k.startswith("cdr") and v}
            for chain, regions in chains.items()
        }
    if region_kind == "scaffold_paratope":
        scaffold_cdr = proteintyper("scaffold_cdr")
        try:
            return scaffold_cdr.extract_scaffold_regions(seq, category) or {}
        except Exception:
            return {}
    return {"W": {FOCUS: seq}}


def identity_fn_for(category: str):
    """Region-identity function for a category.

    DARPin paratopes are per-repeat concatenations, so length-normalised identity would
    cap two designs with different repeat counts well below 1; the library ships a
    repeat-aware identity for exactly this case. Everything else uses the standard
    length-normalised Levenshtein identity.
    """
    cdr_novelty = proteintyper("cdr_novelty")
    if category == "darpin":
        scaffold_search = proteintyper("scaffold_search")
        return scaffold_search.darpin_repeat_identity
    return cdr_novelty.cdr_identity


def focus_region(regions: dict) -> tuple[str | None, str]:
    """Return ``(chain, sequence)`` of the focus region, preferring the heavy chain."""
    for chain in ("H", "S", "W"):
        if chain in regions and regions[chain].get(FOCUS):
            return chain, regions[chain][FOCUS]
    for chain, rs in regions.items():
        if rs.get(FOCUS):
            return chain, rs[FOCUS]
    return None, ""
