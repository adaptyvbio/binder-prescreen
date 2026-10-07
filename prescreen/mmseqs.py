"""Thin MMseqs2 wrappers: one batched search, one clustering call, nothing clever.

Two things here are load-bearing for a prescreen and are the reason this is not just a
copy of the deployed ``dbc`` command:

* **Batching.** One ``easy-search`` call per database for the whole submission batch,
  never one call per sequence. Database load dominates runtime, so N sequences cost
  roughly the same as one.
* **Length-aware settings.** Default MMseqs2 settings silently return nothing for short
  queries (k-mer prefilter), which is how a 25-aa peptide copied verbatim from a paper
  can pass a prescreen as "no match". Short queries are therefore searched with a
  separate, sensitive parameter set, and both passes are merged.

``KMP_AFFINITY=disabled`` is exported for every call: without it MMseqs2 aborts with
``OMP: Error #179: pthread_setaffinity_np() failed`` in this sandbox.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

# ``tseq`` is requested so the hit's own sequence comes back with the alignment. That is
# what makes a paratope-aware prior-art check possible without a second database lookup:
# the hit can be re-numbered and its binding region compared with the query's.
FORMAT_FIELDS = (
    "query,target,fident,alnlen,qlen,tlen,mismatch,qcov,tcov,evalue,bits,tseq"
)
COLUMNS = [
    "query",
    "target",
    "fident",
    "alnlen",
    "qlen",
    "tlen",
    "mismatch",
    "qcov",
    "tcov",
    "evalue",
    "bits",
    "tseq",
]

# The two length-routed parameter sets come from :mod:`prescreen.refdb`, so the curated
# set and the public databases are searched with identical settings and one source of
# truth. ``LENGTH_SPLIT`` is the boundary (50 aa): below it, short-query settings.
def _as_dict(flags) -> dict:
    out, i = {}, 0
    while i < len(flags):
        key = flags[i]
        value = flags[i + 1] if i + 1 < len(flags) and not flags[i + 1].startswith("-") else None
        out[key] = value
        i += 2 if value is not None else 1
    return out


def settings_for_length(n: int) -> dict:
    """Parameter set for a query of length ``n``."""
    from .refdb import LENGTH_SPLIT, PROFILE_SHORT, PROFILE_STANDARD

    return _as_dict(PROFILE_SHORT if n < LENGTH_SPLIT else PROFILE_STANDARD)


def _env() -> dict:
    env = dict(os.environ)
    env["KMP_AFFINITY"] = "disabled"
    return env


def write_fasta(records: dict, path: str | Path) -> Path:
    """Write ``{id: sequence}`` to ``path`` and return it."""
    path = Path(path)
    with path.open("w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")
    return path


def read_fasta(path: str | Path) -> dict:
    """Read a FASTA file into ``{id: sequence}``, keying on the first whitespace token."""
    out: dict = {}
    key = None
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                key = line[1:].split()[0]
                out[key] = ""
            elif key is not None:
                out[key] += line
    return out


def _run(cmd: list, label: str) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        cmd, capture_output=True, text=True, env=_env(), check=False
    )
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "")[-600:]
        raise RuntimeError(f"mmseqs {label} failed ({proc.returncode}): {tail}")
    return proc


def easy_search(
    query_fasta: str | Path,
    target_db: str | Path,
    out_m8: str | Path,
    settings: dict | None = None,
    mmseqs_bin: str = "mmseqs",
    threads: int = 0,
    tmp_dir: str | Path | None = None,
) -> Path:
    """Run ``mmseqs easy-search`` and return the path to the tabular output."""
    out_m8 = Path(out_m8)
    tmp = Path(tmp_dir) if tmp_dir else Path(tempfile.mkdtemp(prefix="prescreen-"))
    tmp.mkdir(parents=True, exist_ok=True)
    cmd = [
        mmseqs_bin,
        "easy-search",
        str(query_fasta),
        str(target_db),
        str(out_m8),
        str(tmp),
        "--format-output",
        FORMAT_FIELDS,
    ]
    if settings is None:
        settings = settings_for_length(10**6)
    for key, value in settings.items():
        cmd += [key] if value is None else [key, str(value)]
    if threads:
        cmd += ["--threads", str(threads)]
    _run(cmd, f"easy-search -> {out_m8.name}")
    return out_m8


def easy_cluster(
    fasta: str | Path,
    out_prefix: str | Path,
    min_seq_id: float = 0.95,
    coverage: float = 0.80,
    mmseqs_bin: str = "mmseqs",
    threads: int = 0,
    tmp_dir: str | Path | None = None,
) -> Path:
    """Cluster a FASTA against itself; returns the ``*_cluster.tsv`` path."""
    out_prefix = Path(out_prefix)
    tmp = Path(tmp_dir) if tmp_dir else Path(tempfile.mkdtemp(prefix="prescreen-clu-"))
    tmp.mkdir(parents=True, exist_ok=True)
    cmd = [
        mmseqs_bin,
        "easy-cluster",
        str(fasta),
        str(out_prefix),
        str(tmp),
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        "1",
    ]
    if threads:
        cmd += ["--threads", str(threads)]
    _run(cmd, "easy-cluster")
    return Path(f"{out_prefix}_cluster.tsv")


def parse_m8(path: str | Path) -> list:
    """Parse an MMseqs2 tabular file written with :data:`FORMAT_FIELDS`."""
    rows = []
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return rows
    numeric = {"fident", "qcov", "tcov", "evalue", "bits"}
    integer = {"alnlen", "qlen", "tlen", "mismatch"}
    with p.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(COLUMNS):
                continue
            row = {}
            for col, val in zip(COLUMNS, parts):
                if col in numeric:
                    row[col] = float(val)
                elif col in integer:
                    row[col] = int(val)
                else:
                    row[col] = val
            rows.append(row)
    return rows
