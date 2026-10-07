"""Prior-art reference databases for the submission prescreen.

Question this module answers: *is a submitted sequence already public prior art?*
It resolves the MMseqs2 target databases, runs the searches with settings that do
not silently lose short queries, and returns one best hit per (query, arm) scored
with the same composite the production ``db_check`` service uses.

It deliberately says nothing about whether a scaffold is de novo: reusing a known
scaffold is allowed in the competition, so only the *submitted sequence itself*
is compared against prior art.

Database root
-------------
Databases are looked up in a list of roots, first match wins:

1. the ``roots`` argument,
2. ``$PRESCREEN_DB_ROOT`` (``:``-separated, like ``PATH``),
3. :data:`FALLBACK_ROOTS`.

so the same code runs against the workstation mount and against a deployed image
with the databases on a volume.

Search settings
---------------
Two length regimes, because MMseqs2's defaults silently drop short queries. On a
120-query PDB panel with known answers, default ``easy-search`` recovered 0/6 exact
matches of 10-14 aa and 69% of all queries under 50 aa; :data:`PROFILE_SHORT`
recovers 100%. The three settings that matter are a relaxed E-value (a 12-residue
exact match scores far above the default 1e-3 cut-off in a million-sequence
database), ``--mask 0`` (low-complexity masking of the target removes the only
k-mers a short query has), and contiguous rather than spaced k-mers at ``-k 6``.

Queries at or above :data:`LENGTH_SPLIT` are already recovered at defaults, so they
are searched with :data:`PROFILE_STANDARD` and cost no extra time.
"""

# ------------------------------------------------------------------------------------
# Vendored into the prescreen package from the reference-database track's deliverable
# (artifact prescreen_refdb.py, version da296348-86bb-431c-8e41-d983a1c26558). One
# deliberate deviation from that original: ``tseq`` is appended to M8_COLUMNS so each
# hit carries the hit's own sequence. The prescreen re-numbers those hit sequences and
# compares their binding regions with the query's, which is how a shared framework is
# told apart from a known molecule. The first eight columns and the
# ``similarity = fident * qcov`` composite are unchanged, so results stay comparable
# with the deployed ``db_check`` service and with that track's benchmarks.
# ------------------------------------------------------------------------------------


from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import pathlib
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "Arm", "ARMS", "PROFILE_SHORT", "PROFILE_STANDARD", "LENGTH_SPLIT",
    "ARM_EXTRA", "ARM_MAX_THREADS", "ReferenceDbMissing", "db_roots", "resolve", "available", "read_fasta",
    "profile_for_length", "search_fasta", "best_prior_art", "manifest", "M8_COLUMNS",
]

# Tried in order when neither the `roots` argument nor $PRESCREEN_DB_ROOT is set.
# Build them with scripts/build_reference_dbs.sh, which writes exactly this layout.
FALLBACK_ROOTS = (
    "dbs",                                  # a ./dbs next to where you are working
    str(pathlib.Path.home() / "prescreen-dbs"),
    "/public-sequence-dbs",                 # a shared mount, if your deployment has one
)
ROOT_ENV = "PRESCREEN_DB_ROOT"

# Same eight fields the production db_check service requests, so its m8 files and
# this module's are interchangeable and the composite is comparable.
M8_COLUMNS = ("query", "target", "fident", "alnlen", "qlen", "mismatch", "qcov", "tcov", "tseq")

LENGTH_SPLIT = 50  # queries shorter than this go through PROFILE_SHORT

PROFILE_STANDARD = (
    "-s", "7.5",
    "-e", "1e-3",
    "--max-seqs", "300",
    "--split-memory-limit", "2G",
)

PROFILE_SHORT = (
    "-s", "7.5",
    "-k", "6",
    "--spaced-kmer-mode", "0",
    "--mask", "0",
    "--min-ungapped-score", "5",
    "--comp-bias-corr", "0",
    "-e", "1e4",
    "--max-seqs", "300",
    "--split-memory-limit", "2G",
)

# The production baseline, kept so a run can be reproduced for comparison.
PROFILE_BASELINE: tuple[str, ...] = ()

# Per-arm additions appended after the profile. The patent arm is 10.2M sequences
# and MMseqs2 refuses to run it in the default footprint: on a 15 GB box it needs
# an explicit target split and must not preload the database, or the prefilter is
# OOM-killed. A smaller --split-memory-limit does NOT help - below about 1 GB
# MMseqs2 aborts with "Cannot fit databases into ...".
ARM_EXTRA: dict[str, tuple[str, ...]] = {
    "uspto": ("--split-memory-limit", "3G", "--db-load-mode", "3"),
}
# Arms whose thread count is capped regardless of the caller's request, because
# per-thread alignment buffers are what pushes the patent arm over the limit.
ARM_MAX_THREADS: dict[str, int] = {"uspto": 4}


@dataclass(frozen=True)
class Arm:
    """One prior-art reference database."""

    name: str
    arm: str           # general | antibody | designs | patent
    visibility: str    # public | internal
    subdir: str = ""   # directory under the root; defaults to `name`
    db: str = ""       # db prefix inside that directory; defaults to `name`
    source: str = ""
    note: str = ""

    @property
    def rel(self) -> str:
        return f"{self.subdir or self.name}/{self.db or self.name}"


ARMS: dict[str, Arm] = {a.name: a for a in [
    Arm("pdb", "general", "public",
        source="https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz",
        note="all PDB seqres chains; ~5.7% are mol:na and never match a protein query"),
    Arm("swissprot", "general", "public",
        source="https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
               "knowledgebase/complete/uniprot_sprot.fasta.gz"),
    Arm("plabdab", "antibody", "public",
        source="https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/paired_sequences.csv.gz",
        note="paired dump only: every entry has a heavy AND a light chain, so it "
             "contains no single-domain antibodies - search plabdab_nano alongside it"),
    Arm("plabdab_nano", "antibody", "public",
        source="https://opig.stats.ox.ac.uk/webapps/plabdab-nano/static/downloads/all_sequences.csv.gz",
        note="VHH / VNAR / sdAb prior art, one record per unique sequence"),
    Arm("therasabdab", "antibody", "public",
        source="https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/"
               "TheraSAbDab_SeqStruc_OnlineDownload.csv",
        note="therapeutic antibody chains incl. the approved anti-TNF biologics"),
    Arm("thpdb", "general", "public",
        source="https://webs.iiitd.edu.in/raghava/thpdb/sequences/allseq",
        note="FDA-approved therapeutic proteins and peptides; upstream FASTA dated 2017"),
    Arm("proteinbase_public", "designs", "public",
        source="proteinbase.com public API (/api/proteins); scripts/build_reference_dbs.sh builds this arm",
        note="published designs from earlier rounds - the designed-binder prior art"),
    Arm("uspto", "patent", "public",
        source="https://ftp.ebi.ac.uk/pub/databases/patentdata/uspto_prt.dat.gz",
        note="patent protein sequences; createdb must be run on a FASTA conversion, "
             "not on the EMBL-style .dat"),
]}

PUBLIC_ARMS = tuple(a.name for a in ARMS.values() if a.visibility == "public")


class ReferenceDbMissing(FileNotFoundError):
    """A reference database an arm is searched against is not present.

    Raised rather than quietly searching whatever is available: a missing arm
    reads exactly like a sequence with no prior art, which is the answer the
    prescreen must never give by accident.
    """


def db_roots(roots=None) -> list[Path]:
    if roots is None:
        env = os.environ.get(ROOT_ENV, "")
        roots = [p for p in env.split(os.pathsep) if p] or list(FALLBACK_ROOTS)
    elif isinstance(roots, (str, Path)):
        roots = [roots]
    return [Path(p).expanduser() for p in roots]


def resolve(name: str, roots=None) -> Path:
    """Return the MMseqs2 db prefix for `name`, or raise ReferenceDbMissing."""
    if name not in ARMS:
        raise KeyError(f"unknown arm {name!r}; known: {sorted(ARMS)}")
    rel = ARMS[name].rel
    tried = []
    for root in db_roots(roots):
        prefix = root / rel
        tried.append(str(prefix))
        if prefix.with_suffix(".dbtype").exists() and prefix.with_suffix(".index").exists():
            return prefix
    raise ReferenceDbMissing(f"arm {name!r} not found; looked for " + ", ".join(tried))


def available(roots=None) -> dict[str, str | None]:
    out = {}
    for name in ARMS:
        try:
            out[name] = str(resolve(name, roots))
        except ReferenceDbMissing:
            out[name] = None
    return out


def read_fasta(path) -> dict[str, str]:
    seqs, key = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                key = line[1:].split()[0]
                seqs[key] = ""
            elif key is not None:
                seqs[key] += line.strip()
    return seqs


def profile_for_length(n: int) -> str:
    return "short" if n < LENGTH_SPLIT else "standard"


_PROFILES = {"short": PROFILE_SHORT, "standard": PROFILE_STANDARD, "baseline": PROFILE_BASELINE}


def _merge_flags(profile, extra):
    """Overlay `extra` on `profile`, replacing rather than repeating a flag.

    MMseqs2 rejects a command line that names the same option twice, so a
    per-arm override of e.g. --split-memory-limit has to displace the profile's
    value instead of being appended after it.
    """
    def pairs(flags):
        out, i = [], 0
        while i < len(flags):
            if i + 1 < len(flags) and not flags[i + 1].startswith("--") and not (
                    flags[i + 1].startswith("-") and flags[i + 1][1:2].isalpha()
                    and not flags[i + 1][1:2].isdigit()):
                out.append((flags[i], flags[i + 1])); i += 2
            else:
                out.append((flags[i], None)); i += 1
        return out

    merged = dict(pairs(profile))
    merged.update(dict(pairs(extra)))
    flat = []
    for k, v in merged.items():
        flat.append(k)
        if v is not None:
            flat.append(v)
    return tuple(flat)


def _easy_search(query_fasta, db_prefix, out_m8, extra, threads, tmp_root):
    tmp = tempfile.mkdtemp(dir=tmp_root)
    cmd = [
        "mmseqs", "easy-search", str(query_fasta), str(db_prefix), str(out_m8), tmp,
        "--format-output", ",".join(M8_COLUMNS), "--threads", str(threads), *extra,
    ]
    env = dict(os.environ, KMP_AFFINITY="disabled")  # else MMseqs2 aborts with OMP #179
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"mmseqs easy-search failed for {db_prefix}:\n{proc.stderr.strip()[-2000:]}"
        )
    return cmd


def search_fasta(fasta, arms=None, roots=None, out_dir=None, threads=8,
                 profile="auto", tmp_root=None, keep_m8=True):
    """Search every sequence in `fasta` against each arm. Returns a DataFrame.

    One row per (query, arm): the best hit by the production composite
    ``similarity = fident * qcov``. Queries with no hit in an arm are absent from
    that arm's rows, so a caller that needs an explicit "no prior art" row should
    reindex against the query set.

    `profile` is "auto" (route each query by length), or one of "short",
    "standard", "baseline" to force a single setting set for every query.
    """
    import pandas as pd

    arms = list(arms or PUBLIC_ARMS)
    prefixes = {a: resolve(a, roots) for a in arms}      # raises before any search
    seqs = read_fasta(fasta)
    if not seqs:
        raise ValueError(f"no sequences read from {fasta}")

    out_dir = Path(out_dir or tempfile.mkdtemp(prefix="prescreen_"))
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_root = str(tmp_root or out_dir)

    if profile == "auto":
        groups = {}
        for k, s in seqs.items():
            groups.setdefault(profile_for_length(len(s)), {})[k] = s
    else:
        groups = {profile: seqs}

    frames, timings = [], []
    for gname, gseqs in groups.items():
        qf = out_dir / f"query_{gname}.fasta"
        with open(qf, "w") as fh:
            for k, s in gseqs.items():
                fh.write(f">{k}\n{s}\n")
        for arm in arms:
            m8 = out_dir / f"{arm}.{gname}.m8"
            import time
            t0 = time.perf_counter()
            extra = _merge_flags(_PROFILES[gname], ARM_EXTRA.get(arm, ()))
            nthreads = min(threads, ARM_MAX_THREADS.get(arm, threads))
            _easy_search(qf, prefixes[arm], m8, extra, nthreads, tmp_root)
            dt = time.perf_counter() - t0
            timings.append({"arm": arm, "profile": gname, "n_queries": len(gseqs),
                            "seconds": round(dt, 2)})
            if m8.exists() and m8.stat().st_size:
                df = pd.read_csv(m8, sep="\t", header=None, names=list(M8_COLUMNS))
                df["arm"] = arm
                df["profile"] = gname
                frames.append(df)
            if not keep_m8:
                m8.unlink(missing_ok=True)

    cols = [*M8_COLUMNS, "arm", "profile", "similarity"]
    if not frames:
        res = pd.DataFrame(columns=cols)
    else:
        res = pd.concat(frames, ignore_index=True)
        res["similarity"] = res["fident"] * res["qcov"]
        res = (res.sort_values("similarity", ascending=False)
                  .drop_duplicates(["query", "arm"])
                  .reset_index(drop=True))
    res.attrs["timings"] = timings
    res.attrs["out_dir"] = str(out_dir)
    return res


def best_prior_art(hits, queries=None):
    """Collapse per-arm hits to one best public-prior-art row per query."""
    import pandas as pd

    pub = hits[hits["arm"].map(lambda a: ARMS[a].visibility == "public")]
    best = (pub.sort_values("similarity", ascending=False)
               .drop_duplicates("query")
               .set_index("query"))
    if queries is None:
        return best.reset_index()
    rows = []
    for q in queries:
        if q in best.index:
            rows.append(best.loc[q].to_dict() | {"query": q})
        else:
            rows.append({"query": q, "arm": "no_match", "similarity": 0.0,
                         "fident": 0.0, "qcov": 0.0, "target": None})
    return pd.DataFrame(rows)


def manifest(path=None):
    p = Path(path or Path(__file__).with_name("reference_manifest.json"))
    return json.loads(p.read_text())


def _cli(argv=None):
    import argparse
    import pandas as pd

    ap = argparse.ArgumentParser(prog="prescreen_refdb", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("search", help="search a FASTA against the prior-art arms")
    s.add_argument("-i", "--input", required=True)
    s.add_argument("-o", "--out-dir", default=None)
    s.add_argument("--arms", nargs="*", default=None)
    s.add_argument("--roots", nargs="*", default=None)
    s.add_argument("--threads", type=int, default=8)
    s.add_argument("--profile", default="auto", choices=["auto", "short", "standard", "baseline"])

    sub.add_parser("arms", help="list arms and whether each database resolves")

    a = ap.parse_args(argv)
    if a.cmd == "arms":
        for name, p in available().items():
            arm = ARMS[name]
            print(f"{name:22s} {arm.arm:9s} {arm.visibility:9s} {p or 'MISSING'}")
        return 0

    hits = search_fasta(a.input, arms=a.arms, roots=a.roots, out_dir=a.out_dir,
                        threads=a.threads, profile=a.profile)
    out_dir = Path(hits.attrs["out_dir"])
    hits.to_csv(out_dir / "prior_art_hits.csv", index=False)
    queries = list(read_fasta(a.input))
    best_prior_art(hits, queries).to_csv(out_dir / "prior_art_best.csv", index=False)
    pd.DataFrame(hits.attrs["timings"]).to_csv(out_dir / "timings.csv", index=False)
    print(f"{len(queries)} queries, {len(hits)} (query, arm) hits -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
