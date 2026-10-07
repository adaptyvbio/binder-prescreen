"""CDR extraction and CDR-novelty analysis for antibodies / nanobodies / scFvs.

This module assumes the input sequence has *already* been identified upstream as an
antibody, nanobody (VHH), or scFv (e.g. via a structural search) and that whole-sequence
similarity has already been computed (mmseqs2 vs SwissProt/PDB/USPTO/THPDB/PLAbDab).

What this module adds is the *CDR-level* answer: given a query and the sequences of its
whole-sequence hits, it extracts the CDRs (primarily CDRH3) and decides whether each CDR
is already present in a known antibody (near-identity match) or is novel.

It relies on antpack (https://github.com/jlparkI/AntPack) for IMGT/Kabat/etc. numbering and
"""
# ----------------------------------------------------------------------------------
# Vendored verbatim from the upstream library (commit a02fae8), except where marked
# `vendored:`. Do not edit to fix a bug here - fix it upstream and re-copy, so both
# copies keep giving the same numbering and the same identities.
# ----------------------------------------------------------------------------------

from __future__ import annotations

from typing import Optional

from antpack import PairedChainAnnotator, SingleChainAnnotator
from loguru import logger

from .novelty_scales import CDR_DE_NOVO_CUTOFF

# Region labels antpack's assign_cdr_labels emits for CDR positions.
_CDR_LABELS = ("cdr1", "cdr2", "cdr3")

# antpack's percent_identity measures similarity to a GERMLINE consensus, so an
# engineered binder scores low while still being a perfectly numbered antibody: the
# reference VHH-like design 01M45FZYSPBWZX7R0R3Q4H1Z4K scores 0.670 and PDB 2YWZ (a
# published VNAR) scores 0.473. This floor therefore only rejects nonsense; whether the
# numbering is trustworthy is decided by the grammar below. The value is still reported.
_MIN_IDENTITY = 0.25

# Conserved IMGT positions that define the immunoglobulin V fold: the two cysteines of
# the canonical disulfide, the tryptophan after CDR1, and the J-segment W/F. A real V
# domain keeps nearly all of them however heavily its loops are engineered, while an
# unrelated fold that antpack numbers anyway keeps at most one or two.
IMGT_ANCHORS = {"23": "C", "41": "W", "104": "C", "118": "WF"}
MIN_ANCHORS = 3

# Plausible CDR lengths. Non-antibodies that clear the anchor test tend to fail here,
# usually with an empty or single-residue CDR3.
CDR_LENGTH_BANDS = {"cdr1": (3, 20), "cdr2": (3, 20), "cdr3": (3, 35)}

_PAIRED_TYPES = {"scfv", "antibody", "mab", "fab"}
_SINGLE_TYPES = {"nanobody", "vhh", "sdab", "vh", "single"}


def numbering_is_valid(seq, numbering, cdr_lengths, percent_identity):
    """Decide whether an antpack numbering describes a real immunoglobulin V domain.

    Two complementary tests, because neither alone separates the classes this service
    types. Germline identity cannot: a published VNAR scores 0.473, *below* a monobody
    at 0.420. Conserved anchors alone cannot either: a de novo miniprotein can match 3
    of 4 by chance, but has no CDR3 to speak of.

    Parameters
    ----------
    seq : str
        The sequence that was numbered.
    numbering : list of str
        antpack numbering, same length as ``seq``, "-" outside the domain.
    cdr_lengths : dict
        ``{"cdr1": int, "cdr2": int, "cdr3": int}`` from the assigned labels.
    percent_identity : float
        antpack's germline identity, used only as a nonsense floor.

    Returns
    -------
    tuple
        ``(is_valid, reason)``. ``reason`` is None when valid, else a short
        explanation so callers can tell "not an antibody" from "numbering failed".
    """
    if not numbering or not any(n != "-" for n in numbering):
        return False, "antpack produced no numbering"
    if percent_identity < _MIN_IDENTITY:
        return False, f"germline identity {percent_identity:.3f} below {_MIN_IDENTITY}"

    position_residue = {pos: aa for pos, aa in zip(numbering, seq) if pos != "-"}
    found = sum(
        1
        for pos, accepted in IMGT_ANCHORS.items()
        if position_residue.get(pos, "") in accepted
    )
    if found < MIN_ANCHORS:
        return False, f"only {found}/{len(IMGT_ANCHORS)} conserved IMGT anchors present"

    for label, (low, high) in CDR_LENGTH_BANDS.items():
        length = cdr_lengths.get(label, 0)
        if not low <= length <= high:
            return False, f"{label} length {length} outside {low}-{high}"

    return True, None


def _canonical_chain(chain_name: str) -> str:
    """Map antpack chain names to a canonical key: 'H' for heavy, 'L' for light (K or L)."""
    return "H" if chain_name == "H" else "L"


def _cdrs_from_chain(seq, numbering, chain_name, percent_identity, scheme, annotator):
    """Extract CDR1/2/3 sequences for one numbered chain.

    Returns a dict ``{"cdr1":..,"cdr2":..,"cdr3":.., "percent_identity":.., "chain_name":..}``
    or ``None`` if the chain is absent or the numbering does not describe a V domain.
    """
    # numbering is the same length as the full input sequence, with "-" outside this domain.
    if not numbering or not any(n != "-" for n in numbering):
        return None
    labels = annotator.assign_cdr_labels(numbering, chain_name, scheme=scheme)

    cdrs = {"cdr1": "", "cdr2": "", "cdr3": ""}
    for residue, label in zip(seq, labels):
        if label in cdrs:
            cdrs[label] += residue

    valid, reason = numbering_is_valid(
        seq, numbering, {k: len(v) for k, v in cdrs.items()}, percent_identity
    )
    if not valid:
        logger.debug(f"Rejected {chain_name} numbering: {reason}")
        return None

    cdrs["percent_identity"] = percent_identity
    cdrs["chain_name"] = chain_name
    return cdrs


def extract_cdrs(seq: str, molecule_type: str = "auto", scheme: str = "imgt") -> dict:
    """Number ``seq`` with antpack and extract its CDRs per chain.

    Args:
        seq: amino-acid sequence (scFv = VH + linker + VL; nanobody = single VHH; etc.).
        molecule_type: one of ``"scfv"``, ``"antibody"``, ``"nanobody"``/``"vhh"``, or
            ``"auto"``. ``"auto"`` tries a paired (VH+VL) interpretation first and falls
            back to a single heavy domain if no light chain is found.
        scheme: numbering scheme ("imgt", "kabat", "martin", "aho").

    Returns
    -------
        ``{"H": {...}, "L": {...}}`` keyed by canonical chain ("H"/"L"). A nanobody / single
        VH yields only ``"H"``. Each value holds ``cdr1``/``cdr2``/``cdr3`` strings plus
        ``percent_identity`` and the antpack ``chain_name``.
    """
    mol = molecule_type.lower()
    result: dict = {}

    def run_paired() -> dict:
        annotator = PairedChainAnnotator(scheme=scheme, receptor_type="mab")
        heavy, light = annotator.analyze_seq(seq)
        out = {}
        for numbering, pid, chain_name, _err in (heavy, light):
            cdrs = _cdrs_from_chain(seq, numbering, chain_name, pid, scheme, annotator)
            if cdrs is not None:
                out[_canonical_chain(chain_name)] = cdrs
        return out

    def run_single() -> dict:
        annotator = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
        numbering, pid, chain_name, _err = annotator.analyze_seq(seq)
        cdrs = _cdrs_from_chain(seq, numbering, chain_name, pid, scheme, annotator)
        if cdrs is None:
            return {}
        # A nanobody/VHH/VNAR is a heavy-type single domain whatever antpack's closest
        # germline says: PDB 2YWZ (a VNAR) is typed "K". Keying it "L" would hide its
        # CDR3 from every consumer that reads the heavy chain, so an explicitly
        # single-domain request always reports "H". "auto" keeps antpack's own call,
        # where a lone light chain is a real possibility.
        key = "H" if mol in _SINGLE_TYPES else _canonical_chain(chain_name)
        return {key: cdrs}

    if mol in _SINGLE_TYPES:
        result = run_single()
    elif mol in _PAIRED_TYPES:
        result = run_paired()
    elif mol == "auto":
        result = run_paired()
        # Fall back to single-domain if the paired view found no/poor heavy chain.
        if "H" not in result:
            single = run_single()
            if single:
                result = single
    else:
        raise ValueError(f"Unknown molecule_type: {molecule_type!r}")

    return result


def _cdr_positions_from_chain(
    seq, numbering, chain_name, percent_identity, scheme, annotator
):
    """Extract 1-based CDR residue positions for one numbered chain.

    Returns ``{"cdr1": [pos, ...], "cdr2": [...], "cdr3": [...]}`` of 1-based indices into
    the full input sequence, or ``None`` if the chain is absent or the numbering does not
    describe a V domain. Mirrors :func:`_cdrs_from_chain` but keeps positions instead of
    residue characters, and applies the same validity grammar.
    """
    if not numbering or not any(n != "-" for n in numbering):
        return None
    labels = annotator.assign_cdr_labels(numbering, chain_name, scheme=scheme)
    positions: dict = {"cdr1": [], "cdr2": [], "cdr3": []}
    for i, label in enumerate(labels):
        if label in positions:
            positions[label].append(i + 1)

    valid, reason = numbering_is_valid(
        seq, numbering, {k: len(v) for k, v in positions.items()}, percent_identity
    )
    if not valid:
        logger.debug(f"Rejected {chain_name} numbering: {reason}")
        return None
    return positions


def cdr_residue_positions(
    seq: str, molecule_type: str = "auto", scheme: str = "imgt"
) -> dict:
    """Number ``seq`` with antpack and return its CDR residue positions per chain.

    Like :func:`extract_cdrs` but returns 1-based residue *positions* into ``seq`` instead of
    CDR strings, for intersecting with structure-derived residue numbers (e.g. interface
    residues).

    Parameters
    ----------
    seq : str
        Amino-acid sequence to number.
    molecule_type : str, optional
        One of "antibody", "fab", "scfv", "nanobody", "vhh" or "auto", by default "auto".
        "auto" tries the paired annotator first and falls back to the single-chain one.
    scheme : str, optional
        Numbering scheme passed to antpack, by default "imgt".

    Returns
    -------
        ``{"H": {"cdr1": [...], "cdr2": [...], "cdr3": [...]}, "L": {...}}`` keyed by canonical
        chain ("H"/"L"), positions 1-based into ``seq``. Empty dict if ``seq`` cannot be
        numbered as an antibody / nanobody / scFv.
    """
    mol = molecule_type.lower()

    def run_paired() -> dict:
        annotator = PairedChainAnnotator(scheme=scheme, receptor_type="mab")
        heavy, light = annotator.analyze_seq(seq)
        out = {}
        for numbering, pid, chain_name, _err in (heavy, light):
            pos = _cdr_positions_from_chain(
                seq, numbering, chain_name, pid, scheme, annotator
            )
            if pos is not None:
                out[_canonical_chain(chain_name)] = pos
        return out

    def run_single() -> dict:
        annotator = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
        numbering, pid, chain_name, _err = annotator.analyze_seq(seq)
        pos = _cdr_positions_from_chain(
            seq, numbering, chain_name, pid, scheme, annotator
        )
        return {_canonical_chain(chain_name): pos} if pos is not None else {}

    if mol in _SINGLE_TYPES:
        return run_single()
    if mol in _PAIRED_TYPES:
        return run_paired()
    if mol == "auto":
        result = run_paired()
        if "H" not in result:
            single = run_single()
            if single:
                result = single
        return result
    raise ValueError(f"Unknown molecule_type: {molecule_type!r}")


def _domain_span_from_chain(
    seq, numbering, chain_name, percent_identity, scheme, annotator
):
    """Return the 1-based ``(start, end)`` extent of one numbered V domain, or ``None``.

    Applies the same validity grammar as :func:`_cdrs_from_chain`, so a chain antpack
    numbers but that is not a real V domain yields ``None`` rather than a span.
    """
    if not numbering or not any(n != "-" for n in numbering):
        return None
    labels = annotator.assign_cdr_labels(numbering, chain_name, scheme=scheme)
    cdr_lengths = dict.fromkeys(_CDR_LABELS, 0)
    for label in labels:
        if label in cdr_lengths:
            cdr_lengths[label] += 1

    valid, reason = numbering_is_valid(seq, numbering, cdr_lengths, percent_identity)
    if not valid:
        logger.debug(f"Rejected {chain_name} numbering: {reason}")
        return None

    numbered = [i for i, n in enumerate(numbering) if n != "-"]
    return numbered[0] + 1, numbered[-1] + 1


def variable_domain_spans(
    seq: str, molecule_type: str = "auto", scheme: str = "imgt"
) -> dict:
    """Number ``seq`` with antpack and return the extent of each V domain it contains.

    Like :func:`cdr_residue_positions` but returns the whole domain rather than its CDRs.
    An scFv carries two V domains in one chain, and the databases every whole-sequence
    search runs against hold single chains, so the two have to be searched separately —
    see :func:`split_variable_domains`.

    Parameters
    ----------
    seq : str
        Amino-acid sequence to number, chain separators already stripped.
    molecule_type : str, optional
        See :func:`extract_cdrs`, by default "auto".
    scheme : str, optional
        Numbering scheme passed to antpack, by default "imgt".

    Returns
    -------
        ``{"H": (start, end), "L": (start, end)}`` keyed by canonical chain, 1-based and
        inclusive. Empty dict if ``seq`` cannot be numbered as a V domain at all.
    """
    mol = molecule_type.lower()

    def run_paired() -> dict:
        annotator = PairedChainAnnotator(scheme=scheme, receptor_type="mab")
        heavy, light = annotator.analyze_seq(seq)
        out = {}
        for numbering, pid, chain_name, _err in (heavy, light):
            span = _domain_span_from_chain(
                seq, numbering, chain_name, pid, scheme, annotator
            )
            if span is not None:
                out[_canonical_chain(chain_name)] = span
        return out

    def run_single() -> dict:
        annotator = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
        numbering, pid, chain_name, _err = annotator.analyze_seq(seq)
        span = _domain_span_from_chain(
            seq, numbering, chain_name, pid, scheme, annotator
        )
        if span is None:
            return {}
        key = "H" if mol in _SINGLE_TYPES else _canonical_chain(chain_name)
        return {key: span}

    if mol in _SINGLE_TYPES:
        return run_single()
    if mol in _PAIRED_TYPES:
        return run_paired()
    if mol == "auto":
        result = run_paired()
        if "H" not in result:
            single = run_single()
            if single:
                result = single
        return result
    raise ValueError(f"Unknown molecule_type: {molecule_type!r}")


def split_variable_domains(
    seq: str, molecule_type: str = "auto", scheme: str = "imgt"
) -> list | None:
    """Split a single-chain construct into its V domains, N- to C-terminal.

    Returns ``None`` unless at least two disjoint V domains are found, so a nanobody, a
    lone VH and anything that is not an antibody are all left for the caller to handle as
    one query. Residues between the domains — the (GGGGS)n linker of an scFv, and any
    tag — are *dropped*: they belong to neither domain, and carrying them into a
    single-chain database search is what drags whole-sequence identity down.

    Parameters
    ----------
    seq : str
        Amino-acid sequence, chain separators already stripped.
    molecule_type : str, optional
        See :func:`extract_cdrs`, by default "auto".
    scheme : str, optional
        Numbering scheme passed to antpack, by default "imgt".

    Returns
    -------
    list or None
        V-domain subsequences in N-to-C order, or ``None``.
    """
    try:
        spans = variable_domain_spans(seq, molecule_type, scheme)
    except Exception:
        return None

    ordered = sorted(spans.values())
    if len(ordered) < 2:
        return None
    # Consecutive domains along one chain cannot overlap. If they do, the numbering
    # placed two domains on the same residues and the split cannot be trusted.
    for (_, prev_end), (next_start, _) in zip(ordered, ordered[1:]):
        if next_start <= prev_end:
            logger.debug(f"Overlapping V domains in {len(seq)}-aa query; not splitting")
            return None
    return [seq[start - 1 : end] for start, end in ordered]


def cdr_interface_fractions(
    binder_seq: str,
    interface_resnums,
    molecule_type: str = "auto",
    scheme: str = "imgt",
) -> dict | None:
    """Fraction of a binder's CDR residues at the interface, and vice versa.

    Args:
        binder_seq: the designed binder amino-acid sequence.
        interface_resnums: iterable of the binder chain's interface residue *numbers*
            (1-based; matching antpack's 1-based sequence positions, i.e. the binder chain is
            numbered sequentially from 1 in the predicted complex).
        molecule_type: see :func:`extract_cdrs`.
        scheme: numbering scheme.

    Returns
    -------
        ``None`` if ``binder_seq`` cannot be numbered as an antibody (both metrics skipped).
        Otherwise ``{"cdr_fraction_interface", "interface_fraction_cdr", "n_cdr",
        "n_interface", "n_overlap"}`` with both fractions on a 0-100 scale:

        * ``cdr_fraction_interface`` = 100 * |CDR ∩ interface| / |CDR|
        * ``interface_fraction_cdr`` = 100 * |CDR ∩ interface| / |interface| (0 if no interface)
    """
    try:
        positions = cdr_residue_positions(binder_seq, molecule_type, scheme)
    except Exception:
        return None

    cdr_positions = {
        pos for cdrs in positions.values() for lst in cdrs.values() for pos in lst
    }
    if not cdr_positions:
        return None

    interface = {int(r) for r in interface_resnums}
    n_cdr = len(cdr_positions)
    n_interface = len(interface)
    n_overlap = len(cdr_positions & interface)

    return {
        "cdr_fraction_interface": 100.0 * n_overlap / n_cdr,
        "interface_fraction_cdr": (100.0 * n_overlap / n_interface)
        if n_interface
        else 0.0,
        "n_cdr": n_cdr,
        "n_interface": n_interface,
        "n_overlap": n_overlap,
    }


def _levenshtein(a: str, b: str) -> int:
    """Classic Levenshtein edit distance (stdlib only)."""
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def cdr_identity(a: str, b: str) -> float:
    """Length-normalized identity for two CDR strings, in [0, 1]."""
    if not a or not b:
        return 0.0
    return 1.0 - _levenshtein(a, b) / max(len(a), len(b))


def compare_regions(
    query_regions: dict,
    hit_regions: list,
    hit_meta: list | None = None,
    cutoff: float = CDR_DE_NOVO_CUTOFF,
    focus: str = "cdr3",
    identity_fn=None,
) -> tuple:
    """Compare query regions against pre-extracted hit regions; return ``(chains, focus_by_chain)``.

    This is the identity-comparison core of :func:`analyze_cdr_novelty`, factored out so
    it can be driven by a non-antibody region extractor as well. It knows nothing about
    antibodies: ``query_regions`` and each entry of ``hit_regions`` are simply
    ``{chain: {region_label: sequence}}`` mappings. :mod:`scaffold_search` uses it with
    affibody / DARPin / monobody paratope regions.

    ``identity_fn`` overrides :func:`cdr_identity` for families whose regions are not
    comparable as plain strings — the DARPin's paratope is a per-repeat concatenation, so
    length normalization would cap two designs with different repeat counts well below 1.
    """
    identity = identity_fn or cdr_identity
    chains: dict = {}
    focus_verdict_by_chain: dict = {}
    for chain, qc in query_regions.items():
        chains[chain] = {}
        for cdr in _CDR_LABELS:
            q_cdr = qc.get(cdr, "")
            if not q_cdr:
                continue
            best_id = 0.0
            best_seq: str | None = None
            best_idx: int | None = None
            for i, hc in enumerate(hit_regions):
                ref = hc.get(chain, {}).get(cdr, "") if hc.get(chain) else ""
                if not ref:
                    continue
                ident = identity(q_cdr, ref)
                if ident > best_id:
                    best_id, best_seq, best_idx = ident, ref, i
            is_novel = best_id < cutoff
            chains[chain][cdr] = {
                "query": q_cdr,
                "closest_hit_cdr": best_seq,
                "identity": round(best_id, 4),
                "closest_hit_index": best_idx,
                "closest_hit": hit_meta[best_idx]
                if (hit_meta and best_idx is not None)
                else None,
                "is_novel": is_novel,
            }
            if cdr == focus:
                focus_verdict_by_chain[chain] = is_novel
    return chains, focus_verdict_by_chain


def analyze_cdr_novelty(
    query_seq: str,
    hit_seqs: list | None = None,
    molecule_type: str = "auto",
    scheme: str = "imgt",
    cutoff: float = CDR_DE_NOVO_CUTOFF,
    focus: str = "cdr3",
    hit_cdrs: list | None = None,
    hit_meta: list | None = None,
    structure_path: str | None = None,
    tmalign_path: str | None = None,
) -> dict:
    """Decide whether the query's CDRs are novel relative to its whole-sequence hits.

    Args:
        query_seq: the query antibody / nanobody / scFv sequence.
        hit_seqs: sequences of the upstream whole-sequence hits to compare against. Each is
            numbered with antpack to extract its CDRs. Ignored if ``hit_cdrs`` is given.
        molecule_type: see :func:`extract_cdrs`.
        scheme: numbering scheme.
        cutoff: identity cutoff for novelty; a CDR is *not novel* if its best identity to any
            hit's same chain+CDR is >= this value.
        focus: which CDR drives the headline verdict (default ``"cdr3"``).
        hit_cdrs: optional pre-extracted hit CDRs as a list of ``{"H": {...}, "L": {...}}``
            dicts (e.g. from :func:`plabdab_search.plabdab_hit_cdrs`, where heavy and light are
            stored separately). If provided, ``hit_seqs`` is not numbered.
        hit_meta: optional list parallel to the hits with metadata (e.g. ``{"id", "title"}``)
            surfaced as ``closest_hit`` on each CDR result.
        structure_path: predicted structure, used only when the numbering cannot place the
            framework; the CDRs are then recovered by structural projection.
        tmalign_path: TMalign binary for that structural fallback.

    Returns
    -------
        ``{"molecule_type", "scheme", "cutoff", "query_cdrs", "chains", "cdr3_novel"}``.
        ``chains`` maps canonical chain -> CDR -> ``{query, closest_hit_cdr, identity,
        closest_hit_index, closest_hit, is_novel}``. ``cdr3_novel`` is driven by the heavy
        chain's focus CDR (CDRH3), falling back to any chain if there is no heavy chain.
    """
    query = extract_cdrs(query_seq, molecule_type, scheme)
    paratope_method = "numbering"
    if not query and structure_path:
        # Numbering could not place this framework. Recover the CDRs structurally rather
        # than abandoning the design: the fold is the same even when the sequence has
        # drifted past what a germline-based numbering recognises.
        try:
            from proteintyper_lib.paratope import identify_paratope
        except ModuleNotFoundError as exc:
            # Not vendored: it pulls in the structure stack (foldseek, TMalign).
            # The prescreen is sequence-only and never passes `structure_path`, so
            # this branch is unreachable there; install the upstream library to use it.
            raise RuntimeError(
                "structural paratope recovery needs the upstream `paratope` module, "
                "which is not vendored into prescreen"
            ) from exc

        resolved = identify_paratope(
            query_seq,
            "nanobody" if molecule_type in _SINGLE_TYPES else molecule_type,
            structure_path=structure_path,
            tmalign_path=tmalign_path,
        )
        if resolved["method"] == "structure":
            query = {resolved["chain_key"]: resolved["regions"]}
            paratope_method = "structure"
            logger.info("CDRs recovered by structural projection")
    if not query:
        raise ValueError(
            "Could not number the query sequence as an antibody/nanobody/scFv."
        )

    # Use pre-extracted hit CDRs if given; otherwise number every hit once (skip uninterpretable).
    if hit_cdrs is None:
        if hit_seqs is None:
            raise ValueError("Provide either hit_seqs or hit_cdrs.")
        hit_cdrs = []
        for hit in hit_seqs:
            try:
                hc = extract_cdrs(hit, molecule_type, scheme)
            except Exception:
                hc = {}
            hit_cdrs.append(hc)

    chains, focus_verdict_by_chain = compare_regions(
        query, hit_cdrs, hit_meta=hit_meta, cutoff=cutoff, focus=focus
    )

    # CDRH3-centric headline: the heavy-chain focus CDR drives the verdict (that's "CDRH3").
    # Fall back to any available chain if there is no heavy chain (shouldn't happen for a
    # successfully numbered antibody/nanobody/scFv).
    undetermined_reason = None
    if "H" in focus_verdict_by_chain:
        cdr3_novel = focus_verdict_by_chain["H"]
    elif molecule_type.lower() in _PAIRED_TYPES:
        # Tables 1 and 2 say CDR-H3 drives the verdict for a Fab/scFv. extract_cdrs only
        # falls back to the single-chain annotator for "auto", so a heavy chain the
        # grammar rejected leaves {"L": ...} behind - and scoring that would silently
        # report CDR-L3 as the paratope. Undetermined is the honest answer.
        cdr3_novel = None
        undetermined_reason = (
            f"no valid heavy chain for a {molecule_type}; "
            f"numbered chains: {sorted(focus_verdict_by_chain)}"
        )
    elif focus_verdict_by_chain:
        cdr3_novel = next(iter(focus_verdict_by_chain.values()))
    else:
        # No focus CDR was comparable. That is an absence of evidence, not evidence of
        # novelty - returning True here would have the service claim a de novo paratope
        # it never established. Callers must treat None as "undetermined".
        cdr3_novel = None
        undetermined_reason = "no comparable focus CDR in any chain"
    # any(), not truthiness of the list: a hit whose chains could not be numbered is
    # appended as {}, so a list of nothing-but-{} is still truthy while offering nothing to
    # compare against. best_id then stays 0.0 and every CDR reads as novel on no evidence.
    if not any(hit_cdrs) and cdr3_novel:
        cdr3_novel = None
        undetermined_reason = "no hit in the comparison set could be numbered"

    return {
        "undetermined_reason": undetermined_reason,
        "paratope_method": paratope_method,
        "molecule_type": molecule_type,
        "scheme": scheme,
        "cutoff": cutoff,
        "focus": focus,
        "query_cdrs": {
            c: {k: v for k, v in d.items() if k in _CDR_LABELS}
            for c, d in query.items()
        },
        "chains": chains,
        "focus_novel_by_chain": focus_verdict_by_chain,
        "cdr3_novel": cdr3_novel,
    }


# --------------------------------------------------------------------------------------
# Synthetic negative controls: all-alanine CDRs (guaranteed-novel) for full antibody,
# nanobody, and scFv. These replace every CDR residue with 'A' while leaving the framework
# intact, so antpack still numbers them but no real antibody has these CDRs.
# --------------------------------------------------------------------------------------

# scFv scaffold: VH + (GGGGS)x3 linker + VL. (Same as the notebook's shy_shark_cypress.)
SCFV_SCAFFOLD = (
    "QVQLQQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSISRDTSKSQVFFKMNSLQTDDTAIYYCARALTYYDYEFAYWGQGTLVTVSA"
    "GGGGSGGGGSGGGGS"
    "DILLTQSPVILSVSPGERVSFSCRASQSIGTNIHWYQQRTNGSPKLLIRYASESISGIPSRFSGSGSGTDFTLSINSVDPEDIADYYCQQNNNWPTTFGAGTKLELK"
)

# Antibody scaffold: a paired VH+VL domain (the notebook's golden_gecko_stone).
ANTIBODY_SCAFFOLD = (
    "SMEVQLLESGGGLVSPGGSLRLSCAASGFTFSYYYMGWVRQAPGKGLEWVSGISPSSGYTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARYYYGYYYSHMDYWGQGTLVTVSS"
    "GGGGSGGGGSGGGGS"
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSRSGLHTFGQGTKLEIK"
)

# Nanobody scaffold: caplacizumab (ALX-0081) VHH, a single heavy variable domain.
NANOBODY_SCAFFOLD = "EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS"


def make_ala_cdr_variant(
    scaffold_seq: str, molecule_type: str, scheme: str = "imgt", fraction: float = 1.0
) -> str:
    """Return ``scaffold_seq`` with a ``fraction`` of each CDR's residues replaced by alanine.

    Framework positions are left untouched so antpack still numbers the result. ``fraction=1.0``
    (default) gives all-alanine CDRs (guaranteed-novel negative control); ``fraction=0.5``
    mutates every other CDR residue (~50%) for a partial-novelty case. Positions are selected
    deterministically (evenly spaced within each CDR).
    """
    mol = molecule_type.lower()
    seq_chars = list(scaffold_seq)
    step = 1 if fraction >= 1.0 else max(1, round(1.0 / fraction))

    def mask_with(annotator, numbering, chain_name):
        if not numbering or all(n == "-" for n in numbering):
            return
        labels = annotator.assign_cdr_labels(numbering, chain_name, scheme=scheme)
        # Group positions per CDR so the fraction applies within each CDR independently.
        by_cdr: dict = {}
        for i, label in enumerate(labels):
            if label in _CDR_LABELS:
                by_cdr.setdefault(label, []).append(i)
        for idxs in by_cdr.values():
            for i in idxs[::step]:
                seq_chars[i] = "A"

    if mol in _SINGLE_TYPES:
        annotator = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
        numbering, _pid, chain_name, _err = annotator.analyze_seq(scaffold_seq)
        mask_with(annotator, numbering, chain_name)
    else:
        annotator = PairedChainAnnotator(scheme=scheme, receptor_type="mab")
        heavy, light = annotator.analyze_seq(scaffold_seq)
        for numbering, _pid, chain_name, _err in (heavy, light):
            mask_with(annotator, numbering, chain_name)

    return "".join(seq_chars)


def ala_cdr_fixtures(scheme: str = "imgt") -> dict:
    """Return ready-made all-Ala-CDR negative controls for each molecule type.

    Returns ``{name: (variant_seq, molecule_type)}``.
    """
    return {
        "ala_cdr_antibody": (
            make_ala_cdr_variant(ANTIBODY_SCAFFOLD, "antibody", scheme),
            "antibody",
        ),
        "ala_cdr_nanobody": (
            make_ala_cdr_variant(NANOBODY_SCAFFOLD, "nanobody", scheme),
            "nanobody",
        ),
        "ala_cdr_scfv": (make_ala_cdr_variant(SCFV_SCAFFOLD, "scfv", scheme), "scfv"),
    }
