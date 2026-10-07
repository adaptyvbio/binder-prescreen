"""Novelty cut points and the 1-4 scoring rules they drive.

Pure logic, no service or orchestration imports, so the calibration can be
checked directly.

Two parallel scales live here. The antibody scale scores a design by how novel its
CDRs are against PLAbDab; the scaffold scale does the same for the randomized
paratope positions of an affibody, DARPin or monobody. Both follow one rule order:
the de-novo test short-circuits first, which is what makes the four levels
partition without overlap.

The level-3/4 cut is per class. That arm asks whether the *framework* is prior art
while the axis measures the *whole sequence*; the two differ by the framework
fraction of the chain, so no single number serves classes whose paratopes occupy
different shares of the molecule. See ``novelty_scale.md`` for the derivation and
the per-class occupancy tables.
"""
# ----------------------------------------------------------------------------------
# Vendored verbatim from the upstream library (commit a02fae8), except where marked
# `vendored:`. Do not edit to fix a bug here - fix it upstream and re-copy, so both
# copies keep giving the same numbering and the same identities.
# ----------------------------------------------------------------------------------

from loguru import logger

# Classification slugs that have CDRs and are eligible for CDR-novelty scoring,
# mapped to the molecule_type expected by the PLAbDab CDR pipeline.
CDR_MOLECULE_TYPES = {
    "nanobody": "nanobody",
    "vnar": "nanobody",  # single variable domain
    "scfv": "scfv",
    "fab": "fab",
}

# CDR novelty cut points, tuned on non-redundant SAbDab.
# CDRH3 best-identity cutoff defining a "de novo" CDR.
CDR_DE_NOVO_CUTOFF = 0.70
# Above this CDRH3 identity the known CDR is essentially unchanged (novelty 1 vs 2).
CDR_KNOWN_IDENTITY = 0.95
# Global max sequence similarity below which a de novo CDR scores 4 (vs 3).
CDR_GLOBAL_NOVEL_SIM = 0.70
# Global max sequence similarity below which a modified known CDR scores 2 (vs 1).
CDR_GLOBAL_KNOWN_SIM = 0.95

# Scaffold (non-antibody) paratope novelty cut points. The four thresholds are the
# antibody ones: validated on 432 designed/randomized affibody, DARPin and monobody
# sequences, where every fully randomized paratope falls at or below 0.467 and every
# parent control sits at 1.000, so any cut in that band separates them.
SCAFFOLD_DE_NOVO_CUTOFF = 0.70
SCAFFOLD_KNOWN_IDENTITY = 0.95
SCAFFOLD_GLOBAL_KNOWN_SIM = 0.95


def scaffold_known_identity(paratope_len):
    """Level-1/2 identity cut for a paratope of ``paratope_len`` residues.

    Identity on an L-residue paratope only takes the values k/L, so a flat 0.95 is not
    always attainable. The affibody's paratope is its 13 mask positions, where the
    highest value below 1.0 is 12/13 = 0.923 — under a flat cut, "known paratope"
    collapses to "exact match", and a design one substitution away from a catalogued
    affibody is scored level 2 rather than level 1. The same holds for a DARPin scored
    over few repeats (2 repeats -> 13/14 = 0.929).

    So the cut is lowered to one modification, (L-1)/L, whenever 0.95 would be stricter
    than that, and left at 0.95 otherwise. The two agree from L=20 up, which covers every
    monobody in the prior-art set (lengths 23-35) — but only since the paratope became
    the full BC+DE+FG span; keyed on FG alone (L=14) it was unreachable too.

    Parameters
    ----------
    paratope_len : int
        Length of the scored paratope region for this query. 0 or None falls back to the
        flat cut, since nothing can be inferred without a region.

    Returns
    -------
    float
        Identity at or above which the paratope counts as unchanged.
    """
    if not paratope_len:
        return SCAFFOLD_KNOWN_IDENTITY
    return min(SCAFFOLD_KNOWN_IDENTITY, (paratope_len - 1) / paratope_len)


# The level-4 global cut is PER FAMILY, and it is derivable rather than tuned.
#
# The level-3/4 distinction asks whether the FRAMEWORK is prior art: de novo paratope on a
# known framework is level 3, de novo paratope on a de novo framework is level 4. Whole-
# sequence identity is only a proxy for framework identity, and it is a SCALED proxy:
#     I_whole = f * F_framework + (1 - f) * c
# where f is the framework fraction of the chain and c the chance identity of a randomized
# paratope against real prior art. Fitted on 246 of 270 calibration designs spanning framework
# divergence 0-50% (affibody n=80, DARPin n=90, monobody n=76; the rest had no significant hit
# or a best hit that could not be projected onto the family framework):
#     affibody  slope 0.746 (f = 0.776), intercept 0.065, r2 = 0.965
#     DARPin    slope 0.859 (f = 0.866), intercept 0.022, r2 = 0.995
#     monobody  slope 0.538 (f = 0.670), intercept 0.157, r2 = 0.856
# The fitted slope recovers the framework fraction, so the SAME framework-identity threshold
# maps to a DIFFERENT whole-sequence cut per class. Holding the antibody/nanobody cut at its
# calibrated 0.70 (with c = 0.240, the CDR-H3 unrelated-pair median from the 16.27M-pair
# SAbDab null) pins that threshold at F* = 0.841, which then gives:
#     derived: nanobody/VNAR 0.700, scFv 0.689, Fab 0.692, affibody 0.692, DARPin 0.744,
#              monobody 0.609
#     adopted: the first four are held at 0.70 (all within 0.01 of derived, and 0.69 vs 0.70
#              gives identical agreement for the affibody), so only DARPin -> 0.74 and
#              monobody -> 0.61 actually change. CDR_GLOBAL_NOVEL_SIM is untouched.
# Agreement with the framework-identity rule on the calibration set: affibody 0.988 (a flat
# 0.70 does equally well, so it is left unchanged), DARPin 1.000 vs 0.867 at a flat 0.70,
# monobody 0.974 vs 0.553 at a flat 0.70.
#
# NOTE the monobody moves DOWN, not up. An earlier revision of this file set it to 0.74 to
# push fully randomized monobodies to level 4; that was wrong under the definition above,
# because those designs keep an intact FNfn10 framework and are therefore level 3. At 0.74
# agreement with the framework rule is 0.513, the worst of the options tested.
#
# PREFERRED FUTURE FORM: score the level-3/4 arm on framework identity directly, against a
# single cut of 0.841 for every class. That removes the framework-fraction scaling, needs no
# per-class table, and widens the level-3 band from 0.085-0.141 (proxy, class-dependent) to a
# uniform 0.159. These proxy cuts exist for compatibility with the current max_seq_sim input.
# Families absent from this map fall back to the antibody value.
SCAFFOLD_GLOBAL_NOVEL_SIM_DEFAULT = 0.70
# The level-3/4 boundary on the framework axis. Set the same way CDR_GLOBAL_NOVEL_SIM was:
# the 95th percentile of the UNRELATED-pair identity distribution (5% false-positive rate).
# On the project's 16,270,660-pair SAbDab null that rule gives 0.7000 for global identity,
# reproducing CDR_GLOBAL_NOVEL_SIM exactly, and for the framework regions it gives:
#     heavy framework (fwh) 0.857    light framework (fwl) 0.876    VH+VL weighted 0.866
# An earlier revision used 0.841 here, which was back-derived by ASSUMING the 0.70 global cut
# was correct for nanobodies instead of reading the framework null directly. Scaffold families
# still lack their own unrelated-framework null, so the antibody heavy value is used for them;
# re-deriving the scaffold cuts at 0.857 instead of 0.841 shifts them by only 0.010-0.014
# (affibody 0.704, DARPin 0.758, monobody 0.618), within the spread of those fits.
SCAFFOLD_FRAMEWORK_KNOWN_IDENTITY = 0.857
CDR_FRAMEWORK_KNOWN_IDENTITY_HEAVY = 0.857
CDR_FRAMEWORK_KNOWN_IDENTITY_LIGHT = 0.876

# Per-class antibody level-4 global cut. Validated separately per class on 360 designs
# (20 cluster-representative SAbDab antibodies per class x 6 framework-divergence levels,
# all CDRs randomized), searched against PDB seqres + SwissProt with self-exclusion:
#     Fab       slope 0.727, intercept 0.063, r2 0.984 -> derived 0.693; 0.70 enforces F=0.876
#     scFv      slope 0.685, intercept 0.098, r2 0.951 -> derived 0.691; 0.70 enforces F=0.879
#     nanobody  slope 0.665, intercept 0.092, r2 0.928 -> derived 0.662; 0.70 enforces F=0.914
# Fab and scFv land within 0.01 of 0.70 and are LEFT UNCHANGED. The nanobody does not: at 0.70
# the cut demands a 0.914 framework identity against the 0.857 target, and 18 of 120 nanobody
# designs with a KNOWN framework are pushed to level 4 instead of 3 (agreement 0.842; at 0.66
# it is 0.950). The cause is structural - a VHH carries three CDRs in one ~122-aa domain, so
# its paratope is 25.2% of the chain against 21.9% for a Fab, and the whole-sequence axis is
# framework identity scaled by the framework fraction.
#
# CAUTION for multi-domain constructs: a Fab or scFv searched as ONE concatenated query only
# aligns one domain against the single-chain databases, which drops whole-sequence identity to
# 0.62-0.63 even with a fully intact framework - below 0.70, so an intact-framework design is
# mis-called level 4. Split such constructs per domain and recombine length-weighted before
# applying these cuts.
CDR_GLOBAL_NOVEL_SIM_BY_CLASS = {
    "fab": 0.70,  # validated, unchanged
    "scfv": 0.70,  # validated, unchanged
    "nanobody": 0.66,  # derived 0.662; 0.70 over-calls novelty for 15% of known-framework VHHs
    "vnar": 0.66,  # mapped to nanobody upstream by CDR_MOLECULE_TYPES
}
SCAFFOLD_GLOBAL_NOVEL_SIM = {
    "affibody": 0.70,  # derived 0.692; held at 0.70, identical agreement (0.988)
    "darpin": 0.74,
    "monobody": 0.61,
}

# Classification slugs that carry a known randomization mask, i.e. that can be
# routed onto the scaffold paratope scale instead of the general scale.
SCAFFOLD_MOLECULE_TYPES = {
    "affibody": "affibody",
    "monobody": "monobody",
    "darpin": "darpin",
}


def is_prior_art(max_seq_sim, cutoff=CDR_GLOBAL_KNOWN_SIM):
    """True when the whole sequence is already in a public database.

    At or above ``cutoff`` the design is prior art whatever its paratope scores, and both
    scales short-circuit to level 1 on it. The two axes cannot disagree honestly: the
    paratope is a sizeable share of the chain in every class scored here — 25% for a VHH,
    22% for the affibody mask, ~21% for the DARPin mask, 33% for the monobody paratope —
    so rewriting it caps whole-sequence identity near 0.75, and even the narrowest case,
    a CDR-H3-only rewrite at ~11% of a VHH chain, caps it near 0.89. A "de novo paratope"
    reported alongside >= 0.95 global identity is therefore never a novel design; it means
    the paratope database does not hold the parent, which is what
    :func:`cdr_novelty_score` warns about.

    Parameters
    ----------
    max_seq_sim : float or None
        Best whole-sequence similarity in [0, 1]. None (the search produced no usable
        result) is not prior art — nothing was established either way.
    cutoff : float, optional
        Identity at or above which the sequence counts as known.

    Returns
    -------
    bool
    """
    return max_seq_sim is not None and max_seq_sim >= cutoff


def cdr_novelty_score(cdr_result, max_seq_sim, molecule_type=None):
    """Antibody-specific novelty (1-4) from the CDR verdict and global sequence similarity.

    - 4: de novo CDR and high global novelty (max seq similarity < 70%)
    - 3: de novo CDR but global sequence similarity >= 70%
    - 2: known CDR with changes and global sequence similarity < 95%
    - 1: known CDR, or the whole sequence is prior art (global similarity >= 95%)
    - None: the CDR verdict is undetermined, so this scale does not apply. The caller
      falls back to the general sequence+structure scale rather than guessing.
    """
    verdict = cdr_result.get("cdr3_novel")
    if verdict is None:
        return None
    de_novo = bool(verdict)

    # Prior art outranks the paratope verdict; see is_prior_art for why the two cannot
    # honestly disagree. The warning is the signal that a prior-art database is incomplete:
    # it is how the PLAbDab paired-only gap surfaced, where every nanobody read as de novo
    # because no single-domain antibody was in the comparison set at all.
    if is_prior_art(max_seq_sim):
        if de_novo:
            logger.warning(
                f"Paratope called de novo but the whole sequence matches a known one at "
                f"{max_seq_sim:.3f}; scoring as prior art (level 1). The prior-art "
                f"database for {molecule_type!r} is missing this design's parent."
            )
        return 1

    # Best CDRH3 identity to a known antibody: prefer heavy chain, else any chain.
    chains = cdr_result.get("chains", {})
    cdr3 = (chains.get("H") or {}).get("cdr3")
    if cdr3 is None:
        for chain in chains.values():
            if chain.get("cdr3"):
                cdr3 = chain["cdr3"]
                break
    cdr_id = (cdr3 or {}).get("identity", 0.0)

    # Per-class level-4 cut. Fab and scFv are validated at CDR_GLOBAL_NOVEL_SIM; the nanobody
    # needs a lower value because its paratope is a larger share of the chain. See the
    # CDR_GLOBAL_NOVEL_SIM_BY_CLASS comment for the calibration.
    if molecule_type in CDR_GLOBAL_NOVEL_SIM_BY_CLASS:
        global_novel_cut = CDR_GLOBAL_NOVEL_SIM_BY_CLASS[molecule_type]
    else:
        # Falling back here silently reverts a calibrated per-class cut. For a nanobody
        # that matters: 0.70 pushes 18 of 120 known-framework VHH designs to level 4,
        # which is the whole reason Table 3 lowers it to 0.66.
        global_novel_cut = CDR_GLOBAL_NOVEL_SIM
        logger.warning(
            f"No per-class novelty cut for molecule_type {molecule_type!r}; "
            f"falling back to {CDR_GLOBAL_NOVEL_SIM}. Known classes: "
            f"{sorted(CDR_GLOBAL_NOVEL_SIM_BY_CLASS)}"
        )

    # max_seq_sim None means the global axis is unavailable (a failed or error-shaped
    # similarity search). The paratope axis is still valid, so the design is still scored —
    # but only on what the global axis could have ESTABLISHED. It is the sole evidence for
    # level 4, and it can only promote a known paratope to level 1, so without it the
    # answer is 3 on the de novo arm and the paratope clause alone on the known arm.
    # Discarding the paratope verdict instead would send a design with a demonstrably
    # known paratope to the general scale, which reads a missing search as "no similarity"
    # and scores it 4 — the most novel verdict there is, on no evidence at all.
    if de_novo:
        return 3 if max_seq_sim is None or max_seq_sim >= global_novel_cut else 4
    if cdr_id < CDR_KNOWN_IDENTITY and (
        max_seq_sim is None or max_seq_sim < CDR_GLOBAL_KNOWN_SIM
    ):
        return 2
    return 1


def scaffold_novelty_score(scaffold_result, max_seq_sim, family=None):
    """Paratope novelty (1-4) for a scaffold design. Same rule order as :func:`cdr_novelty_score`.

    The level-4 arm is reachable only when whole-sequence similarity is below this family's
    entry in ``SCAFFOLD_GLOBAL_NOVEL_SIM`` (see the note on that constant for the measured
    per-family calibration). With the adopted cuts, a fully randomized monobody reaches
    level 4 while a fully randomized affibody or DARPin does not and stops at level 3 — for
    those two no cut on this axis separates them from known designs.

    Returns None when the paratope verdict is undetermined, so an unprojectable mask is
    never scored as de novo.
    """
    verdict = scaffold_result.get("cdr3_novel")
    if verdict is None:
        return None
    de_novo = bool(verdict)

    # Resolved before the prior-art branch below, which names the family in its warning.
    if family is None:
        family = scaffold_result.get("molecule_type")

    # Prior art outranks the paratope verdict; see is_prior_art and the equivalent branch
    # in cdr_novelty_score.
    if is_prior_art(max_seq_sim, SCAFFOLD_GLOBAL_KNOWN_SIM):
        if de_novo:
            logger.warning(
                f"Paratope called de novo but the whole sequence matches a known one at "
                f"{max_seq_sim:.3f}; scoring as prior art (level 1). The {family!r} "
                f"prior-art database is missing this design's parent."
            )
        return 1

    chains = scaffold_result.get("chains", {})
    paratope = (chains.get("S") or {}).get("cdr3")
    if paratope is None:
        for chain in chains.values():
            if chain.get("cdr3"):
                paratope = chain["cdr3"]
                break
    region_id = (paratope or {}).get("identity", 0.0)

    global_novel_cut = SCAFFOLD_GLOBAL_NOVEL_SIM.get(
        family, SCAFFOLD_GLOBAL_NOVEL_SIM_DEFAULT
    )

    # None max_seq_sim: see the equivalent branch in cdr_novelty_score.
    if de_novo:
        return 3 if max_seq_sim is None or max_seq_sim >= global_novel_cut else 4
    # Length-aware so the cut is actually attainable; see scaffold_known_identity. The
    # length that matters is the one cdr_identity normalized by -- max(query, hit) -- not
    # the query's alone, or a short query against a longer member gets a cut looser than
    # the one modification the rule is meant to allow.
    known_cut = scaffold_known_identity(
        max(
            len((paratope or {}).get("query") or ""),
            len((paratope or {}).get("closest_hit_cdr") or ""),
        )
    )
    if region_id < known_cut and (
        max_seq_sim is None or max_seq_sim < SCAFFOLD_GLOBAL_KNOWN_SIM
    ):
        return 2
    return 1
