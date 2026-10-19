"""Sequence-level category classification for competition submissions.

The prescreen must decide, from sequence alone and in milliseconds, which *kind* of
binder a submission is, because that decides which binding region is compared against
the known-binder reference set. Nothing here touches a structure: typing a molecule by
TM-align against reference PDBs is both slower and unavailable at submission time, when
no model has been folded yet.

Decision order (first match wins):

1. **Antibody formats** — antpack numbering via ``cdr_novelty.extract_cdrs``. That call
   already enforces the library's numbering-validity test (IMGT anchors, CDR length
   bands, germline percent identity), so a successful extraction *is* the evidence that
   the sequence is an antibody variable domain. Paired (H+L) versus single-domain is read
   off the returned chains; the subformat (Fv / scFv / Fab / IgG) is then read off the
   presence of constant-domain segments.
2. **Alternative scaffolds** — framework identity to the canonical affibody (protein A Z
   domain), monobody (FNfn10) and DARPin (NI3C consensus) references from
   ``scaffold_cdr``, gated at the library's ``SCAFFOLD_FRAMEWORK_IDENTITY`` (0.85).
3. **Length regime** — peptide, miniprotein, or small/large protein, with no claim about
   fold. These get a whole-sequence comparison only.

The category also carries a ``region_kind``, which is the contract with
:mod:`prescreen.regions`: ``antibody_cdrs``, ``scaffold_paratope`` or ``whole``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict

from .compat import antpack_annotators, proteintyper

# Human constant-domain probes (germline reference segments) used only as presence/absence
# landmarks to separate Fv / scFv / Fab / IgG. Truncated to the first 24 residues: enough
# to be diagnostic, short enough that a single window scan is cheap.
CONSTANT_PROBES = {
    "CH1": "ASTKGPSVFPLAPSSKSTSGGTAA",
    "HINGE_CH2": "EPKSCDKTHTCPPCPAPELLGGPS",
    "CH3": "GQPREPQVYTLPPSRDELTKNQVS",
    "CL_KAPPA": "RTVAAPSVFIFPPSDEQLKSGTAS",
    "CL_LAMBDA": "GQPKAAPSVTLFPPSSEELQANKA",
}
PROBE_IDENTITY = 0.70

# Flexible-linker patterns that mark a single-chain paired construct (scFv).
LINKER_MOTIFS = ("GGGGSGGGGS", "GGGGSGGGGSGGGGS", "GSTSGSGKPGSGEGSTKG", "GGSGGSGGS")

# Length regimes (residues). Peptide/miniprotein boundaries follow the structural
# classifier's: isPeptide at <= 30 aa, isMiniprotein at < 15 kDa (~135 aa).
PEPTIDE_MAX = 30
MINIPROTEIN_MAX = 135
SMALL_PROTEIN_MAX = 300

#: Shortest sequence that could hold an antibody variable domain. Below this, numbering
#: producing nothing is the expected answer and not worth explaining.
MIN_DOMAIN_LENGTH = 80

#: Confidence for a category decided on length alone, after the antibody and scaffold
#: tests both declined to describe the sequence. Deliberately not 1.0: nothing was
#: identified, a bucket was chosen.
LENGTH_ONLY_CONFIDENCE = 0.5

# Camelid VHH hallmark residues at IMGT positions 42, 49, 50, 52 (Vincke et al. 2009
# J Biol Chem 284:3273; the "tetrad" that replaces the VH-VL interface of a human VH).
VHH_HALLMARKS = {"42": "FY", "49": "EQ", "50": "RC", "52": "GLF"}
VHH_HUMAN_VH = {"42": "VA", "49": "GA", "50": "LMI", "52": "WY"}
VHH_HALLMARK_MIN = 2

ANTIBODY_CATEGORIES = (
    "fv",
    "scfv",
    "fab",
    "igg",
    "nanobody",
    "vh_domain",
    "single_domain_antibody",
)
SCAFFOLD_CATEGORIES = ("affibody", "monobody", "darpin", "centyrin", "fc_tfr")
#: An Fc-family call only when the Fc is most of the sequence: a binder fused to an Fc must
#: keep its own category, or its paratope would be replaced by the Fc patch.
FC_MAX_LENGTH_FRACTION = 1.15

# ``framework_identity`` normalises by the number of ALIGNED framework positions, not by
# the reference length, so it carries no coverage term: a query that aligns to a handful
# of positions and matches them scores 1.0. A single residue "A" reaches DARPin identity
# 1.0 and clears the 0.85 gate. This is the same trap the prior-art arm documents for
# ``fident`` without ``qcov`` — identity alone is meaningless without the span it covers.
# Require the query to be a plausible length for the family before the gate can fire; a
# sequence well under its reference simply is not that scaffold.
SCAFFOLD_MIN_LENGTH_FRACTION = 0.6

# The Proteinbase submission template carries a ``molecule_class`` the submitter declares
# (https://proteinbase.com/templates/competition-submission-template.csv). It is a
# cross-check, never a substitute: trusting it would let a mis-declared submission
# redirect which region is compared, which is exactly what the germline floor below
# exists to prevent. Each declared class maps to the categories that agree with it —
# deliberately loose, because a VHH whose hallmark tetrad is mutated types as
# ``vh_domain`` and that is not a mis-declaration.
DECLARED_CLASS_MAP = {
    "nanobody": {"nanobody", "vh_domain", "single_domain_antibody"},
    "vhh": {"nanobody", "vh_domain", "single_domain_antibody"},
    "sdab": {"nanobody", "vh_domain", "single_domain_antibody"},
    "vnar": {"single_domain_antibody", "vh_domain", "nanobody"},
    "scfv": {"scfv", "fv"},
    "fv": {"fv", "scfv"},
    "fab": {"fab", "fv"},
    "fab_kappa": {"fab", "fv"},
    "fab_lambda": {"fab", "fv"},
    "igg": {"igg", "fab", "fv"},
    "igg_kappa": {"igg", "fab", "fv"},
    "igg_lambda": {"igg", "fab", "fv"},
    "heavy_chain": {"igg", "fab", "fv", "vh_domain", "nanobody"},
    "light_chain": {"igg", "fab", "fv", "single_domain_antibody"},
    "antibody": {"igg", "fab", "fv", "scfv", "nanobody", "vh_domain"},
    "affibody": {"affibody"},
    "monobody": {"monobody", "centyrin"},
    "centyrin": {"centyrin", "monobody"},
    "fc": {"fc_tfr"},
    "darpin": {"darpin"},
    "peptide": {"peptide", "miniprotein"},
    "miniprotein": {"miniprotein", "peptide", "small_protein"},
    "minibinder": {"miniprotein", "peptide", "small_protein"},
    # "single_chain" in the template means one chain of anything that is not an antibody
    # format, so it agrees with every non-antibody category.
    "single_chain": {
        "peptide", "miniprotein", "small_protein", "unclassified",
        "affibody", "monobody", "darpin", "centyrin", "fc_tfr",
    },
}


def compare_declared(declared: str | None, category: str) -> str | None:
    """Compare a submitter's declared ``molecule_class`` with the assigned category.

    Returns ``"agree"``, ``"mismatch"``, ``"unknown"`` (the declared class is not one we
    know how to compare), or ``None`` when nothing was declared.
    """
    if not declared:
        return None
    allowed = DECLARED_CLASS_MAP.get(declared.strip().lower())
    if allowed is None:
        return "unknown"
    return "agree" if category in allowed else "mismatch"

# antpack returns a numbering for *any* input — it force-fits the sequence to the closest
# antibody scheme and reports a germline percent-identity. A non-antibody (e.g. a TNFR
# cysteine-rich ectodomain, with or without an Fc fusion as in etanercept) numbers at
# <= 0.37 identity; genuine variable domains sit far above this — on the curated anti-TNF
# set every real antibody chain's best domain is >= 0.48 (the lone VNAR, a divergent shark
# single-domain), with VHH >= 0.83 and Fv >= 0.74. A numbered domain is therefore accepted
# as an antibody only if its BEST chain's germline identity clears this floor. A detected
# constant domain alone is NOT sufficient: an Fc-fusion of a non-antibody binder (etanercept
# is TNFR2 + IgG1 Fc) carries a real Fc but no paratope, so it must not be typed as an
# antibody and compared on CDRs it does not have.
ANTIBODY_GERMLINE_FLOOR = 0.45


def _numbering_reason(seq: str) -> str | None:
    """Why ``extract_cdrs`` produced nothing for a sequence long enough to be a domain.

    ``cdr_novelty.numbering_is_valid`` already returns a reason - "only 2/4 conserved IMGT
    anchors present", "cdr3 length 41 outside 3-35" - but every caller inside the vendored
    module discards it. Re-running the check here is cheap (one numbering pass) and turns
    a silent fall-through into something a reviewer can act on: an anchor failure means
    this is not an antibody, while a CDR-length failure means it is one the grammar would
    not accept, which is a very different thing to be told.
    """
    cdr_novelty = proteintyper("cdr_novelty")
    try:
        SingleChainAnnotator, _ = antpack_annotators()
        annotator = SingleChainAnnotator(chains=["H", "K", "L"], scheme="imgt")
        numbering, pid, chain_name, err = annotator.analyze_seq(seq)
        if err:
            return err
        labels = annotator.assign_cdr_labels(numbering, chain_name)
        lengths: dict = {}
        for label, residue in zip(labels, seq):
            if label.startswith("cdr"):
                lengths[label] = lengths.get(label, 0) + 1
        valid, reason = cdr_novelty.numbering_is_valid(seq, numbering, lengths, pid)
        return None if valid else reason
    except Exception as exc:  # never let diagnostics break a classification
        return repr(exc)


@dataclass
class Classification:
    """What kind of binder a submitted sequence is, and how to compare it."""

    category: str
    region_kind: str
    confidence: float
    length: int
    evidence: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return asdict(self)


def _best_window_identity(seq: str, probe: str) -> float:
    """Best fraction of identical residues between ``probe`` and any window of ``seq``."""
    n, m = len(seq), len(probe)
    if m == 0 or n < m:
        return 0.0
    best = 0
    for i in range(n - m + 1):
        hits = 0
        for a, b in zip(seq[i : i + m], probe):
            if a == b:
                hits += 1
        if hits > best:
            best = hits
    return best / m


def constant_domains(seq: str) -> dict:
    """Which antibody constant-domain landmarks are present in ``seq``."""
    return {
        name: round(ident, 3)
        for name, probe in CONSTANT_PROBES.items()
        if (ident := _best_window_identity(seq, probe)) >= PROBE_IDENTITY
    }


def has_linker(seq: str) -> str | None:
    for motif in LINKER_MOTIFS:
        if motif in seq:
            return motif
    return None


def imgt_positions(seq: str, scheme: str = "imgt") -> dict:
    """Map IMGT position label to residue for a single antibody domain."""
    SingleChainAnnotator, _ = antpack_annotators()
    annotator = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
    numbering, _pid, _chain, _err = annotator.analyze_seq(seq)
    return {
        num: aa for num, aa in zip(numbering, seq) if num and num not in ("-", "0")
    }


def vhh_hallmarks(seq: str) -> dict:
    """Score the camelid hallmark tetrad (IMGT 42, 49, 50, 52) of a single heavy domain.

    Returns a dict with the observed residues, how many match VHH-type residues, and how
    many match the human-VH consensus. A lone heavy domain with VHH hallmarks is a
    nanobody; one with VH residues is an isolated VH of a conventional antibody.
    """
    try:
        pos = imgt_positions(seq)
    except Exception as exc:
        return {"observed": {}, "matches": 0, "vh_matches": 0, "error": repr(exc)}
    observed = {p: pos.get(p, "") for p in VHH_HALLMARKS}
    matches = sum(1 for p, aa in observed.items() if aa and aa in VHH_HALLMARKS[p])
    vh_matches = sum(1 for p, aa in observed.items() if aa and aa in VHH_HUMAN_VH[p])
    return {"observed": observed, "matches": matches, "vh_matches": vh_matches}


def _antibody_subformat(seq: str, chains: dict) -> tuple[str, dict]:
    """Name the antibody format given the numbered chains and constant-domain landmarks."""
    const = constant_domains(seq)
    linker = has_linker(seq)
    paired = "H" in chains and "L" in chains
    evidence = {
        "chains": sorted(chains),
        "constant_domains": const,
        "linker": linker,
        "germline_identity": {
            c: round(float(chains[c].get("percent_identity") or 0.0), 3) for c in chains
        },
    }
    if paired:
        if "HINGE_CH2" in const or "CH3" in const:
            return "igg", evidence
        if "CH1" in const or "CL_KAPPA" in const or "CL_LAMBDA" in const:
            return "fab", evidence
        if linker:
            return "scfv", evidence
        return "fv", evidence
    # Single domain. A lone heavy domain is only a VHH if it carries the camelid
    # hallmark residues; an isolated VH from a human IgG numbers identically but is not a
    # nanobody, and the distinction matters because the two are searched against
    # different PLAbDab subsets. antpack types some VNARs as light-like, so a lone
    # light-type domain is reported generically rather than guessed into a subtype.
    if "H" in chains:
        hallmarks = vhh_hallmarks(seq)
        evidence["vhh_hallmarks"] = hallmarks
        if hallmarks["matches"] >= VHH_HALLMARK_MIN:
            return "nanobody", evidence
        return "vh_domain", evidence
    return "single_domain_antibody", evidence


def classify(seq: str, scaffold_gate: float | None = None) -> Classification:
    """Classify one submitted sequence.

    Args:
        seq: amino-acid sequence, uppercase, no gaps.
        scaffold_gate: framework-identity floor for an alternative-scaffold call.
            Defaults to ``scaffold_search.SCAFFOLD_FRAMEWORK_IDENTITY`` (0.85).

    Returns:
        A :class:`Classification`. ``category`` is ``"unclassified"`` only for sequences
        that are neither an antibody domain nor a known scaffold nor short enough for a
        length call, i.e. proteins above ``SMALL_PROTEIN_MAX``.
    """
    seq = seq.strip().upper().replace("*", "")
    n = len(seq)

    cdr_novelty = proteintyper("cdr_novelty")
    scaffold_cdr = proteintyper("scaffold_cdr")
    scaffold_search = proteintyper("scaffold_search")
    gate = (
        scaffold_search.SCAFFOLD_FRAMEWORK_IDENTITY
        if scaffold_gate is None
        else scaffold_gate
    )

    # 1. Antibody variable domains.
    try:
        chains = cdr_novelty.extract_cdrs(seq, molecule_type="auto")
    except Exception as exc:  # antpack rejects some inputs outright
        chains = {}
        antibody_error = repr(exc)
    else:
        antibody_error = None
    antibody_reject = None
    # Why the numbering produced nothing. cdr_novelty computes this reason and then drops
    # it into a debug log that compat._quiet_loguru() disables, so it never reaches a
    # caller - which is how a grammar rule could reject real designed antibodies for a
    # year without anyone seeing it. Re-derive it here, for a sequence long enough to hold
    # a V domain, so "could not check" is distinguishable from "checked, not an antibody".
    numbering_reason = None
    if not chains and not antibody_error and n >= MIN_DOMAIN_LENGTH:
        numbering_reason = _numbering_reason(seq)
    if chains:
        category, evidence = _antibody_subformat(seq, chains)
        pids = [float(c.get("percent_identity") or 0.0) for c in chains.values()]
        max_pid = max(pids) if pids else 0.0
        # Accept the numbering as a real antibody only if the best chain clears the germline
        # floor; otherwise it is a force-fit onto a non-antibody and we fall through to the
        # scaffold / length branches. Confidence is the BEST chain's identity (a spurious
        # second chain at low identity must not drag it down). A constant domain is reported
        # in the evidence and drives the subformat, but does not on its own make an Fc-fusion
        # of a non-antibody binder count as an antibody.
        if max_pid >= ANTIBODY_GERMLINE_FLOOR:
            return Classification(
                category=category,
                region_kind="antibody_cdrs",
                confidence=round(max_pid, 3),
                length=n,
                evidence=evidence,
            )
        antibody_reject = {
            "chains": list(chains),
            "max_germline_identity": round(max_pid, 3),
            "floor": ANTIBODY_GERMLINE_FLOOR,
        }

    # 2. Alternative scaffolds, by framework identity to the canonical reference, but only
    # for a query long enough to be one (see SCAFFOLD_MIN_LENGTH_FRACTION).
    fw = {
        family: round(scaffold_cdr.framework_identity(family, seq), 4)
        for family in SCAFFOLD_CATEGORIES
    }
    min_lengths = {
        family: int(SCAFFOLD_MIN_LENGTH_FRACTION * len(scaffold_cdr.SCAFFOLD_REFS[family]))
        for family in SCAFFOLD_CATEGORIES
    }
    eligible = {f: v for f, v in fw.items() if n >= min_lengths[f]}
    if n > FC_MAX_LENGTH_FRACTION * len(scaffold_cdr.FC_REF):
        eligible.pop("fc_tfr", None)
    best_family = max(eligible, key=eligible.get) if eligible else None
    if best_family is not None and fw[best_family] >= gate:
        return Classification(
            category=best_family,
            region_kind="scaffold_paratope",
            confidence=fw[best_family],
            length=n,
            evidence={
                "framework_identity": fw,
                "gate": gate,
                "min_length": min_lengths,
            },
        )

    # 3. Length regime.
    if n <= PEPTIDE_MAX:
        category = "peptide"
    elif n <= MINIPROTEIN_MAX:
        category = "miniprotein"
    elif n <= SMALL_PROTEIN_MAX:
        category = "small_protein"
    else:
        category = "unclassified"
    return Classification(
        category=category,
        region_kind="whole",
        # A length regime is not an identification. Reporting 1.0 here said "certainly a
        # small protein" about a sequence the antibody and scaffold tests had both
        # declined to describe, which reads as confidence about the wrong thing.
        confidence=LENGTH_ONLY_CONFIDENCE if category != "unclassified" else 0.0,
        length=n,
        evidence={
            "framework_identity": fw,
            "gate": gate,
            "min_length": min_lengths,
            "antibody_numbering": (
                "rejected"
                if antibody_error
                else "below_germline_floor"
                if antibody_reject
                else "no_valid_domain"
            ),
            "antibody_reject": antibody_reject,
            "antibody_numbering_reason": numbering_reason,
        },
    )
