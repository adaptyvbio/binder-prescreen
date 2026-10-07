"""Within-batch duplicate detection.

Submissions arrive in waves and the same sequence (or a one-residue variant of it) is
often submitted repeatedly, sometimes by different entrants. Clustering the batch against
itself is one MMseqs2 call and makes those groups explicit before a reviewer reads the
report.
"""

from __future__ import annotations

from pathlib import Path

from . import mmseqs
from .config import Config


def cluster_batch(records: dict, cfg: Config, workdir: str | Path) -> dict:
    """Cluster a submission batch against itself.

    Returns ``{query_id: {"cluster_representative": id, "cluster_size": int}}``.
    Singletons are reported with size 1 and themselves as representative.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out = {
        qid: {"cluster_representative": qid, "cluster_size": 1} for qid in records
    }
    if len(records) < 2:
        return out
    fasta = mmseqs.write_fasta(records, workdir / "batch.fasta")
    tsv = mmseqs.easy_cluster(
        fasta,
        workdir / "batch_clu",
        min_seq_id=cfg.batch_cluster_identity,
        coverage=cfg.batch_cluster_coverage,
        mmseqs_bin=cfg.mmseqs_bin,
        threads=cfg.threads,
        tmp_dir=workdir / "tmp_clu",
    )
    members: dict = {}
    with Path(tsv).open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            rep, member = parts
            members.setdefault(rep, []).append(member)
    for rep, group in members.items():
        for member in group:
            if member in out:
                out[member] = {
                    "cluster_representative": rep,
                    "cluster_size": len(group),
                }
    return out
