---
name: binder-prescreen
description: Run the sequence-level binder prescreen in this repository — screen submitted protein designs for prior art and for known binders of the target. Use when asked to screen, triage or prescreen a batch of designed sequences (FASTA or CSV), to check whether a design copies a known binder or carries a known paratope, to run a prior-art or patent search over candidate binders, or to build the MMseqs2 reference databases the screen needs. Covers installing the package, building the reference arms, running the CLI or the Python API, and checking which arms were actually searched.
---

# Running the binder prescreen

This skill is the mechanics of driving the tool in this repository. For what the verdicts
mean and the cut points behind them, read `README.md` and `PRESCREEN_SPEC.md`.

## 1. Install

From the repository root:

```bash
uv venv && uv pip install -e .
```

MMseqs2 must be on `$PATH` — it is not a pip dependency:

```bash
micromamba install -c conda-forge -c bioconda mmseqs2
# or the static binary:
curl -sSL https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz | tar xz
export PATH="$PWD/mmseqs/bin:$PATH"
```

Do not "upgrade" `antpack`. It is pinned to `0.3.8.6.2`; 0.4 put the numbering tools behind
a license key and imports a Qt GUI that fails on a headless machine. The package sets
`KMP_AFFINITY=disabled` itself — without it MMseqs2 aborts with an OpenMP affinity error.

Check the install with `uv run prescreen --help`.

## 2. Pick a mode before you run anything

| | needs | cost |
|---|---|---|
| `--skip-prior-art` | nothing downloaded | ~25 s for 4 sequences |
| full screen | `$PRESCREEN_DB_ROOT` with the reference arms built | ~70 s for 4 sequences |

`--skip-prior-art` runs the whole **target** arm — whole sequence against the curated
binder set that ships in `prescreen/data/`, the paratope comparison, and the verbatim
fragment check. Reach for it when the question is "is this a known binder", or when the
databases are not built yet. It cannot answer "is this sequence already public".

The full screen adds the public corpora and answers both. Use it when prior art matters.

## 3. Build the reference databases (only for the full screen)

```bash
bash scripts/build_reference_dbs.sh ~/prescreen-dbs                  # all eight arms
SKIP_USPTO=1 bash scripts/build_reference_dbs.sh ~/prescreen-dbs     # skip the patent arm
bash scripts/build_reference_dbs.sh ~/prescreen-dbs thpdb pdb        # or name the arms
export PRESCREEN_DB_ROOT=~/prescreen-dbs
```

No credentials for any arm. Everything except `uspto` is ~1.7 GB and a couple of minutes on
a fast connection; `uspto` is the expensive one. An arm that is already built is skipped, so
an interrupted run resumes — `FORCE=1` rebuilds. Each arm is verified by a real search
before the script calls it done.

**If an arm fails, the script keeps going and exits non-zero.** Each arm comes from a
different third-party host and any of them can be down. Read the summary it prints: it names
the failed arms and gives the command to retry just those. Do not report a screen as complete
while an arm is missing — rerun the failed arm first, or say which arms were skipped.

`$PRESCREEN_DB_ROOT` takes a `:`-separated list, first match wins, so a shared mount and a
few locally built arms can be mixed.

## 4. Run it

```bash
uv run prescreen submissions.fasta -o report/
uv run prescreen submissions.fasta --skip-prior-art -o report/
uv run prescreen submissions.csv --id-column id --sequence-column sequence -o report/
uv run prescreen submissions.fasta --patent -o report/        # + the USPTO arm, seconds per query
uv run prescreen submissions.fasta -o report/ --flagged-only  # drop the passes from report.csv
```

**Pass the whole batch in one call. Never loop over sequences.** The search is batched — one
MMseqs2 call per database for the entire input — so database load dominates and N sequences
cost about what one does. A loop multiplies the cost by N and, because within-batch
clustering only sees one sequence at a time, silently loses the `batch_duplicate` flag.

Output lands in the directory given to `-o`: `report.csv` (one flat row per submission) and
`report.json` (every hit, every threshold, every intermediate).

From Python, when you need the records rather than the files:

```python
from prescreen import Config, screen, to_rows

out  = screen({"sub_001": "EVQLVESGGG..."}, cfg=Config.from_env())
rows = to_rows(out)        # the report.csv schema, one dict per submission
```

`--target` (or `target_fasta=`) points the known-binder arm at a different curated set, so
the same filter works for any target, not just the one that ships.

## 5. Before reporting a clean result, check what was searched

An arm whose database is not on disk is **skipped silently**. Nothing is printed, and a
skipped arm reads exactly like a sequence with no prior art — so a "nothing found" result
from an incomplete database root is not evidence of anything.

```bash
uv run python -c "from prescreen import refdb; print(refdb.available())"   # before
```

```python
pa = json.load(open("report/report.json"))["results"]["<id>"]["prior_art"]
pa["arms_searched"], pa["arms_missing"]                             # after
```

If `arms_missing` is non-empty, say so when reporting the result, and name the arms that
were skipped.
