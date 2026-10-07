"""Search a query scaffold design against a scaffold-family reference set and judge paratope novelty.

Scaffold analogue of :mod:`plabdab_search`. Where PLAbDab supplies "every known
antibody" for CDR novelty, this module supplies "every known member of this
scaffold family" for paratope novelty, built from PDB seqres and SwissProt.

Pipeline (mirrors ``plabdab_search`` stage for stage):
  1. ``mmseqs easy-search`` the query against the scaffold seed db -> family call.
  2. Confirm with framework identity: a designed scaffold variant keeps the family
     framework intact (>= ``SCAFFOLD_FRAMEWORK_IDENTITY``) while rewriting the mask
     positions; a natural family homologue has a diverged framework and is sent to
     the general novelty scale instead.
  3. ``mmseqs easy-search`` the query against the family reference db -> candidate members.
  4. Project the randomization mask onto each member to get its paratope regions.
  5. Compare query vs member regions with ``cdr_novelty.compare_regions`` and report.

Only step 4 is scaffold-specific. Steps 3 and 5 reuse the antibody machinery, and the
1-4 novelty scoring is the antibody rule set with renamed constants.
"""
# ----------------------------------------------------------------------------------
# Vendored verbatim from the upstream library (commit a02fae8), except where marked
# `vendored:`. Do not edit to fix a bug here - fix it upstream and re-copy, so both
# copies keep giving the same numbering and the same identities.
# ----------------------------------------------------------------------------------

from __future__ import annotations

import os
import subprocess
import tempfile

from .cdr_novelty import cdr_identity, compare_regions
from .scaffold_cdr import (
    DARPIN_MASK,
    SCAFFOLD_REFS,
    extract_scaffold_regions,
    framework_identity,
)

# Prebuilt mmseqs databases: the three framework seeds, and one reference db per family.
# vendored: the upstream default is a fixed deployment path; here it follows
# $PRESCREEN_DB_ROOT so the scaffold databases sit beside the prior-art arms. Note that
# the prescreen itself only reads this module's constants and darpin_repeat_identity -
# none of the searching functions below are on its path, so these databases are optional.
DEFAULT_DB_DIR = os.path.join(
    os.environ.get("PRESCREEN_DB_ROOT", "dbs").split(os.pathsep)[0], "scaffolds"
)
SEED_DB_NAME = "scaffold_seeds"
DEFAULT_MMSEQS = "mmseqs"

# A designed scaffold variant keeps its framework; a natural homologue does not.
# Calibrated on 432 designed/randomized scaffold sequences vs 194 natural family
# homologues from PDB: the band around 0.85 is empty in all three families
# (natural max 0.841-0.850, designed min 0.852-0.867).
SCAFFOLD_FRAMEWORK_IDENTITY = 0.85

_M8_COLUMNS = ("query", "target", "fident", "alnlen", "qlen", "evalue", "qcov", "tcov")


def darpin_repeat_identity(query, hit):
    """Repeat-count-independent identity between two DARPin paratope concatenations.

    A DARPin's paratope is ``len(DARPIN_MASK)`` residues per repeat, concatenated, and the
    repeat count varies between designs. Comparing the concatenations as strings normalizes
    by ``max(len(a), len(b))``, so two designs with an *identical* per-repeat paratope but
    different repeat counts cap out well below 1 — 3-vs-4 repeats at 0.75, and 2-vs-3 at
    0.667, which is below the 0.70 de novo cut. A design that copies a known paratope and
    adds one repeat would therefore read as de novo.

    Instead each query repeat is matched against its best counterpart in the hit and the
    per-repeat identities are averaged, which returns 1.0 whenever every query repeat is
    reproduced somewhere in the hit, whatever the counts.
    """
    size = len(DARPIN_MASK)
    q = [query[i : i + size] for i in range(0, len(query), size)]
    h = [hit[i : i + size] for i in range(0, len(hit), size)]
    q = [r for r in q if r]
    h = [r for r in h if r]
    if not q or not h:
        return 0.0
    return sum(max(cdr_identity(r, s) for s in h) for r in q) / len(q)


class ScaffoldGateError(ValueError):
    """The query is not a scaffold variant this scale applies to.

    Distinct from a genuine failure: a natural family homologue, or a design whose
    class and framework disagree, has a *determinate* answer — the scaffold scale does
    not apply and the general sequence+structure scale does. Callers catch this
    separately from ``Exception`` so the two are not logged as the same thing.
    """


def _easy_search(
    query_seq, db_path, mmseqs_bin=DEFAULT_MMSEQS, sensitivity=7.5, max_seqs=20000
):
    """Run ``mmseqs easy-search`` for one sequence and return parsed rows."""
    with tempfile.TemporaryDirectory() as tmp:
        q_fasta = os.path.join(tmp, "query.fasta")
        out_m8 = os.path.join(tmp, "out.m8")
        with open(q_fasta, "w") as fh:
            fh.write(f">query\n{query_seq}\n")
        cmd = [
            mmseqs_bin,
            "easy-search",
            q_fasta,
            db_path,
            out_m8,
            os.path.join(tmp, "mmtmp"),
            "--format-output",
            ",".join(_M8_COLUMNS),
            "-s",
            str(sensitivity),
            "--max-seqs",
            str(max_seqs),
            "-e",
            "1e-3",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"mmseqs easy-search failed:\n{proc.stderr.strip()}")
        rows = []
        with open(out_m8) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) == len(_M8_COLUMNS):
                    rows.append(dict(zip(_M8_COLUMNS, cols)))
    return rows


def detect_scaffold_family(
    query_seq, db_dir=DEFAULT_DB_DIR, mmseqs_bin=DEFAULT_MMSEQS, **kw
):
    """Return ``(family, framework_identity)`` for the best seed hit, gate NOT applied.

    ``family`` is whichever family's seed the query matches, whatever its framework
    identity, so a caller can tell "no hit at all" (``(None, None)``) apart from
    "matched, but the framework has diverged". The framework gate itself lives in
    :func:`analyze_scaffold_novelty`, so the family a design is measured against is
    always the family its mask is projected with.
    """
    rows = _easy_search(query_seq, os.path.join(db_dir, SEED_DB_NAME), mmseqs_bin, **kw)
    if not rows:
        return None, None
    best = min(rows, key=lambda r: float(r["evalue"]))
    family = best["target"].split("|")[0]
    if family not in SCAFFOLD_REFS:
        return None, None
    return family, framework_identity(family, query_seq)


def search_scaffold_family(
    query_seq,
    family,
    db_dir=DEFAULT_DB_DIR,
    mmseqs_bin=DEFAULT_MMSEQS,
    max_hits=20000,
    **kw,
):
    """Return family members similar to the query, best e-value first, with their regions.

    ``max_hits`` must stay above every family's prior-art set, which the patent corpus grew
    to 9,997 affibody, 4,114 monobody and 1,028 DARPin members. The cap was 2000 while the
    comment here still claimed it cleared sets of a few hundred, so the two larger families
    were being silently truncated. The ordering it truncates on is e-value — whole-sequence
    similarity, which among members sharing an identical framework ranks essentially
    arbitrarily with respect to the *paratope* similarity the verdict depends on.
    """
    rows = _easy_search(query_seq, os.path.join(db_dir, family), mmseqs_bin, **kw)
    best = {}
    for r in rows:
        t = r["target"]
        if t not in best or float(r["evalue"]) < float(best[t]["evalue"]):
            best[t] = r
    hits = sorted(best.values(), key=lambda r: float(r["evalue"]))[:max_hits]
    out = []
    for h in hits:
        seq = _member_sequence(h["target"], family, db_dir)
        if not seq:
            continue
        try:
            regions = extract_scaffold_regions(seq, family)
        except Exception:
            continue
        out.append(
            {
                "id": h["target"],
                "fident": float(h["fident"]),
                "evalue": float(h["evalue"]),
                "regions": regions,
            }
        )
    return out


def _member_sequence(target, family, db_dir):
    """Look up a reference member's sequence from the family FASTA shipped with the db."""
    cache = _member_sequence._cache.setdefault(family, None)
    if cache is None:
        cache = {}
        fasta = os.path.join(db_dir, f"{family}.fasta")
        if os.path.exists(fasta):
            cur, buf = None, []
            with open(fasta) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if cur:
                            cache[cur] = "".join(buf)
                        cur, buf = line[1:].split()[0], []
                    else:
                        buf.append(line.strip())
            if cur:
                cache[cur] = "".join(buf)
        _member_sequence._cache[family] = cache
    return cache.get(target)


_member_sequence._cache = {}


def analyze_scaffold_novelty(
    query_seq,
    family=None,
    cutoff=0.70,
    focus="cdr3",
    exclude_ids=(),
    db_dir=DEFAULT_DB_DIR,
    structure_path=None,
    tmalign_path=None,
    **kw,
):
    """End-to-end paratope-novelty call for a non-antibody scaffold design.

    Owns the scaffold gate: the query must match a family seed *and* keep that family's
    framework, and when ``family`` is given (from the design's class) the detected family
    must agree with it. Anything else raises :class:`ScaffoldGateError`, which means
    "use the general novelty scale", not "something broke".

    ``structure_path`` and ``tmalign_path`` drive the structural fallback used when the
    sequence mask will not project. They are explicit parameters rather than part of
    ``**kw``, which stays reserved for :func:`_easy_search` tuning (``sensitivity``,
    ``max_seqs``).

    Returns the same dict shape as ``cdr_novelty.analyze_cdr_novelty`` with ``db`` set
    to the scaffold family and an added ``framework_identity``, so downstream scoring
    and the metric payload need no new shape.
    """
    detected, fw = detect_scaffold_family(query_seq, db_dir=db_dir, **kw)
    if family is None:
        family = detected
    if family is None:
        raise ScaffoldGateError("query does not match a supported scaffold family")
    if detected is None:
        raise ScaffoldGateError(
            f"no significant hit to the {family} framework seed set"
        )
    if detected != family:
        raise ScaffoldGateError(
            f"classified as {family} but the framework matches {detected}"
        )
    if fw < SCAFFOLD_FRAMEWORK_IDENTITY:
        raise ScaffoldGateError(
            f"{family} framework identity {fw:.3f} below "
            f"{SCAFFOLD_FRAMEWORK_IDENTITY}: a natural homologue, not a scaffold variant"
        )

    query_regions = extract_scaffold_regions(query_seq, family)
    paratope_method = "mask"
    if not (query_regions.get("S", {}).get("cdr3") or "").strip("-") and structure_path:
        # The mask would not project. For a DARPin that happens whenever the repeat
        # consensus is perturbed, since detection leans on a regex; recover the paratope
        # from the fold instead of reporting an empty one.
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
            query_seq, family, structure_path=structure_path, tmalign_path=tmalign_path
        )
        if resolved["method"] == "structure":
            query_regions = {"S": resolved["regions"]}
            paratope_method = "structure"
    members = [
        m
        for m in search_scaffold_family(query_seq, family, db_dir=db_dir, **kw)
        if m["id"] not in exclude_ids
    ]
    hit_regions = [m["regions"] for m in members]
    hit_meta = [
        {"id": m["id"], "fident": m["fident"], "evalue": m["evalue"]} for m in members
    ]

    chains, focus_by_chain = compare_regions(
        query_regions,
        hit_regions,
        hit_meta=hit_meta,
        cutoff=cutoff,
        focus=focus,
        identity_fn=darpin_repeat_identity if family == "darpin" else None,
    )
    # None, not True: an unprojectable mask leaves no region to compare, and that is an
    # absence of evidence rather than evidence of a de novo paratope.
    novel = focus_by_chain.get("S")
    undetermined_reason = None
    if novel is None:
        undetermined_reason = "paratope could not be projected onto the query"
    elif not members:
        novel, undetermined_reason = None, "comparison set was empty"
    elif not any((hr.get("S") or {}).get(focus) for hr in hit_regions):
        # extract_scaffold_regions SUCCEEDS with empty regions for a member whose mask will
        # not project (a DARPin db entry the repeat regex misses, say), so members can be
        # non-empty while nothing in it is actually comparable. Without this the query is
        # called de novo against a set that offered no paratope to compare with.
        novel, undetermined_reason = (
            None,
            "no member of the comparison set had a usable paratope",
        )
    return {
        "undetermined_reason": undetermined_reason,
        "paratope_method": paratope_method,
        "molecule_type": family,
        "scheme": "scaffold_mask",
        "cutoff": cutoff,
        "focus": focus,
        "query_cdrs": {c: dict(d) for c, d in query_regions.items()},
        "chains": chains,
        "focus_novel_by_chain": focus_by_chain,
        "cdr3_novel": novel,
        "framework_identity": fw,
        "db": f"scaffold:{family}",
        "n_hits": len(members),
    }
