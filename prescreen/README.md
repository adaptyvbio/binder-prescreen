# `prescreen` — sequence-level prior-art & known-binder filter

A submission-time prescreen for protein-design binder competitions. For each submitted
sequence it answers two questions, both from sequence alone and both fast enough to run
when a submission arrives:

1. **Is this sequence already public prior art?** Batched MMseqs2 search against PDB,
   SwissProt, PLAbDab (+ PLAbDab-nano for single-domain formats), Thera-SAbDab, THPdb,
   the Proteinbase design corpus, and 10.2 M USPTO patent sequences.
2. **Is it a known binder of the target, or does it carry a known binding region?**
   Whole-sequence search against a curated known-binder set, plus a framework-independent
   comparison of the binding region (antibody CDR3, projected scaffold paratope) against
   the known paratopes.

It **does not** judge whether the scaffold is de novo. Reusing a published framework is
allowed in these competitions, so the molecule category (nanobody, DARPin, affibody,
monobody, miniprotein, peptide, …) is used only to choose which binding region to compare.
A novel scaffold carrying a genuinely new paratope against the same epitope **passes by
design** — see the blind-spots section of `PRESCREEN_SPEC.md`.

## Why not just the existing whole-sequence service

The deployed `dbc` / `sequence_similarity` service answers question 1 only, at default
MMseqs2 sensitivity. This package differs where it matters for a prescreen:

- **Short sequences are not silently missed.** Default MMseqs2 settings return nothing for
  a 10–30-aa peptide (too few k-mers survive the prefilter). Queries below 50 aa are
  searched with a sensitive parameter set; measured recovery under 20 aa went from 31 % to
  100 % with no normal-length case made worse.
- **The binding region is compared, not just the whole sequence.** The known-binder corpus
  cross-hits itself at the framework level — one probe VHH returned 1,843 hits at ≥0.80
  coverage, 647 of them at ≥0.70 identity, all framework background (only 15 reached ≥0.90
  identity) — so a whole-sequence cut cannot separate "known binder" from "known
  framework". The CDR / paratope arm is what makes that call.
- **A known framework with a new paratope passes; a known paratope on a new framework is
  flagged.** This is the core requirement and the thing the whole-sequence service gets
  wrong in both directions.

## Install & data

```
pip install -e .    # antpack==0.3.8.6.2 (pinned, pre-license-gate), biopython, click,
                    # loguru, pandas
```

The package is self-contained. The four numbering and paratope modules it shares with its
upstream library — `cdr_novelty` (antpack IMGT numbering, CDR identity), `scaffold_cdr`
(affibody / monobody / DARPin masks), `scaffold_search` and `novelty_scales` — are vendored
verbatim in `prescreen/_proteintyper/`, so the same sequence is numbered the same way either
side without the prescreen needing that library on disk. Point `$PROTEINTYPER_SRC` at an
upstream checkout's `src` directory to run against it instead and check for drift.

The curated TNF-α reference set ships in `prescreen/data/` (`tnfa_binders.fasta`,
`tnfa_binder_cdrs.csv`). The prior-art databases are resolved from a `:`-separated root
list — `PRESCREEN_DB_ROOT`, then built-in fallbacks (the workstation mount, then a
deployed volume). See `reference_manifest.json` for per-arm counts, provenance and
snapshot dates, and `build_reference_dbs.sh` to rebuild any arm.

## Use

```python
from prescreen import Config, screen, to_rows
out = screen({"sub_001": "QVQLVESGG..."}, cfg=Config.from_env())
rows = to_rows(out)          # one flat dict per submission
```

```
prescreen submissions.fasta -o report/              # FASTA or CSV in; report.csv + report.json out
prescreen submissions.csv -o report/ --flagged-only
```

The USPTO patent arm (10.2 M sequences) is searched on every run: a prior-art screen
that skips patents is not a prior-art screen. It costs seconds per query, amortised over
the batch; `--no-patent` (or `Config(include_patent_arm=False)`) turns it off.

## Verdicts

The verdict is the highest-ranked flag that fired (`FLAG_ORDER` in `flags.py`); every flag
that fired is listed, and framework-reuse observations are reported separately as
`context` and never disqualify. Full definitions, cut points and the calibration behind
them are in `PRESCREEN_SPEC.md`.
