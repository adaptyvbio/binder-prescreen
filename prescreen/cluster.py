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
    # Synthetic ids, for the reason :mod:`prescreen.chains` gives for the search arms:
    # MMseqs2 keys a FASTA accession on the first whitespace token, so a submission named
    # "my fab v2" comes back from easy-cluster as "my", matches no key below, and is
    # silently reported as a singleton — while "my fab v3" collides with it on the same
    # truncated header. An opaque id cannot collide with anything the submitter typed.
    safe = {f"b{i}": seq for i, seq in enumerate(records.values())}
    owner = dict(zip(safe, records))
    fasta = mmseqs.write_fasta(safe, workdir / "batch.fasta")
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
        # Back to submission ids before anything leaves this function.
        rep_id = owner.get(rep, rep)
        for member in group:
            member_id = owner.get(member, member)
            if member_id in out:
                out[member_id] = {
                    "cluster_representative": rep_id,
                    "cluster_size": len(group),
                }
    return out
