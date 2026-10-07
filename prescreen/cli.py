"""Batch CLI: ``prescreen <submissions> --target <binders.fasta> -o <outdir>``.

Reads FASTA or CSV, screens the whole batch in one pass, and writes a flat CSV plus the
full JSON record (every hit, every threshold, every piece of evidence).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import click

from . import Config, mmseqs, screen, to_rows
from .flags import FLAG_ORDER


def _read_submissions(
    path: Path, id_column: str, sequence_column: str, class_column: str
) -> tuple[dict, dict]:
    """Read submissions, returning ``({id: sequence}, {id: declared_class})``.

    CSV follows the Proteinbase competition template
    (https://proteinbase.com/templates/competition-submission-template.csv):
    ``name,sequence,molecule_class``, with the chains of a multi-chain entry joined by
    ``:``. A legacy ``id`` column is still accepted. ``molecule_class`` is optional and
    is only cross-checked against the classifier, never trusted in its place.
    """
    if path.suffix.lower() not in (".csv", ".tsv"):
        return mmseqs.read_fasta(path), {}

    delim = "\t" if path.suffix.lower() == ".tsv" else ","
    records: dict = {}
    declared: dict = {}
    with path.open(newline="") as fh:
        for i, row in enumerate(csv.DictReader(fh, delimiter=delim)):
            seq = (row.get(sequence_column) or "").strip()
            if not seq:
                continue
            # The template names the column "name"; accept "id" so a CSV written against
            # the older convention still runs.
            key = ""
            for col in (id_column, "name", "id"):
                key = (row.get(col) or "").strip()
                if key:
                    break
            key = key or f"row{i + 1}"
            if key in records:
                # Silently overwriting would drop a submission from the report without
                # anyone noticing it was ever screened.
                raise click.ClickException(
                    f"duplicate submission name {key!r} in {path} (row {i + 2})"
                )
            records[key] = seq
            cls = (row.get(class_column) or "").strip()
            if cls:
                declared[key] = cls
    return records, declared


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("submissions", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--target",
    "target_fasta",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Curated known-binder FASTA for the target (default: packaged TNF-alpha set).",
)
@click.option(
    "--target-metadata",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="CSV of per-reference annotations (id, class, source, citation...).",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("prescreen_out"),
    help="Output directory for report.csv and report.json.",
)
@click.option("--db-root", default=None, help="Prior-art database root directory.")
@click.option(
    "--skip-prior-art",
    is_flag=True,
    help="Run only the target-binder arm (no public database search).",
)
@click.option(
    "--no-patent",
    "skip_patent_arm",
    is_flag=True,
    help="Skip the USPTO patent arm (10.2M sequences; ~6-7 min per batch, amortised over "
    "all queries). Searched by default — without it, a binder claimed only in a patent "
    "is invisible.",
)
@click.option("--threads", default=0, type=int, help="MMseqs2 threads (0 = auto).")
@click.option("--id-column", default="name", help="CSV id column.")
@click.option("--sequence-column", default="sequence", help="CSV sequence column.")
@click.option(
    "--class-column",
    default="molecule_class",
    help="CSV column holding the submitter's declared molecule class; cross-checked "
    "against the classifier, never used in its place.",
)
@click.option(
    "--flagged-only", is_flag=True, help="Write only submissions whose verdict is not 'pass'."
)
def main(
    submissions: Path,
    target_fasta: Path | None,
    target_metadata: Path | None,
    output: Path,
    db_root: str | None,
    skip_prior_art: bool,
    skip_patent_arm: bool,
    threads: int,
    id_column: str,
    sequence_column: str,
    class_column: str,
    flagged_only: bool,
) -> None:
    """Screen SUBMISSIONS against prior art and known target binders.

    SUBMISSIONS is a FASTA file, or a CSV in the Proteinbase competition template format
    (name,sequence,molecule_class; chains of a multi-chain entry joined by ":").
    """
    records, declared = _read_submissions(
        submissions, id_column, sequence_column, class_column
    )
    if not records:
        click.echo("no sequences read", err=True)
        sys.exit(1)

    cfg = Config.from_env()
    if db_root:
        cfg.db_root = db_root
    cfg.threads = threads
    # Only override the config when the flag was actually given, so the env var survives.
    if skip_patent_arm:
        cfg.include_patent_arm = False

    output.mkdir(parents=True, exist_ok=True)
    screened = screen(
        records,
        cfg=cfg,
        target_fasta=target_fasta,
        target_metadata=target_metadata,
        workdir=output / "work",
        skip_prior_art=skip_prior_art,
        declared=declared,
    )

    rows = to_rows(screened)
    if flagged_only:
        rows = [r for r in rows if r["verdict"] != "pass"]
    csv_path = output / "report.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["id"])
        writer.writeheader()
        writer.writerows(rows)
    (output / "report.json").write_text(json.dumps(screened, indent=1))

    counts = {}
    for rec in screened["results"].values():
        counts[rec["verdict"]] = counts.get(rec["verdict"], 0) + 1
    click.echo(f"screened {len(records)} sequences -> {csv_path}")
    for flag in FLAG_ORDER:
        if flag in counts:
            click.echo(f"  {flag:28s} {counts[flag]}")
    click.echo("timing (s): " + json.dumps(screened["timing"]))


if __name__ == "__main__":
    main()
