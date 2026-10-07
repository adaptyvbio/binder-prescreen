# TNF-α binder prescreen

A sequence-level prior-art and known-binder filter for protein-design binder
competitions. Given a submitted sequence it answers two questions, from sequence alone:

1. **Is this already public?** MMseqs2 against PDB, SwissProt, PLAbDab (+ PLAbDab-nano
   for single-domain formats), Thera-SAbDab, THPdb, the published Proteinbase designs,
   and — off the fast path — 10.2 M USPTO patent sequences.
2. **Is it a known TNF-α binder?** Whole-sequence *and* **binding-region** comparison
   against a curated set of 4,964 known or claimed anti-TNF binders that ships with the
   package.

It deliberately does **not** judge whether the scaffold is de novo. Reusing a published
framework is allowed in these competitions, so the molecule category (nanobody, Fab, IgG,
affibody, monobody, DARPin, miniprotein, peptide) is used only to decide which binding
region to compare.

## Why the binding region is compared separately

Antibodies are most of the known anti-TNF corpus and nearly all share the same germline
framework. Measured on the reference set, a single probe VHH matches **1,843 of ~5,000
entries** at ≥0.80 coverage and 647 of them at ≥0.70 identity — all framework background,
with only 15 genuine near-identical hits. So "is this similar to a known binder over its
whole length?" answers *yes* for almost any antibody-format design, whether or not it
shares the actual binding site.

Comparing the paratope independently fixes both error directions:

| submission | whole sequence | paratope | verdict |
|---|---|---|---|
| known framework, **new CDR3** | ~90 % identical to a known binder | no match | **passes** — a new binder, allowed |
| **known CDR3**, fresh framework | looks unrelated | matches | **flagged** — a known binder in disguise |

That discrimination is the whole point of the tool. Cut points, the calibration behind
them (1,036 labelled sequences), and the deliberate blind spots are in
[`PRESCREEN_SPEC.md`](PRESCREEN_SPEC.md); [`HOW_IT_WORKS.md`](HOW_IT_WORKS.md) is the
plain-language walkthrough.

## Install

```bash
git clone https://github.com/adaptyvbio/binder-prescreen.git && cd binder-prescreen
uv venv && uv pip install -e .     # antpack (pinned), biopython, click, loguru, pandas

# MMseqs2 — conda, or the static binary
micromamba install -c conda-forge -c bioconda mmseqs2
# or: curl -sSL https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz | tar xz
#     export PATH="$PWD/mmseqs/bin:$PATH"
```

`antpack` is pinned to `0.3.8.6.2`: 0.4 moved the numbering tools behind a license key and
imports a Qt GUI on a headless box.

## Screen something in 30 seconds, with no databases

`--skip-prior-art` turns off the public-corpus search and runs the TNF-specific arm only:
whole-sequence against the curated binder set, the paratope comparison, and the verbatim
fragment check. The reference set is inside the package, so this needs **no downloads**.

```bash
uv run prescreen examples/submissions.fasta --skip-prior-art -o report/
```

```
screened 4 sequences -> report/report.csv
  known_target_binder          3
  pass                         1
```

~25 s for 4 sequences on 8 cores. That catches a known binder, a known paratope grafted
onto a new framework, and a verbatim fragment. What it cannot tell you is whether the
sequence is *public prior art* — for that you need the corpora.

## Build the reference databases

```bash
bash scripts/build_reference_dbs.sh ~/prescreen-dbs                  # all eight arms
SKIP_USPTO=1 bash scripts/build_reference_dbs.sh ~/prescreen-dbs     # skip the patent arm
bash scripts/build_reference_dbs.sh ~/prescreen-dbs thpdb pdb        # or name the arms
export PRESCREEN_DB_ROOT=~/prescreen-dbs
```

Needs only `mmseqs`, `curl`, `python3` and network access — **no credentials for any
arm**. Each arm lands at `<root>/<arm>/<arm>` as an MMseqs2 database, which is the layout
`prescreen.refdb.resolve()` expects. Arms are built cheapest first, an arm that is already
built is skipped (so an interrupted run resumes, and `FORCE=1` rebuilds), and the 3.4 GB
patent download resumes rather than restarting.

Every arm is also **indexed** (`mmseqs createindex`). This is what makes the patent arm
usable: the same 4-sequence `--patent` batch takes **2m17s indexed against 13m17s without**,
for identical verdicts. On the `uspto` arm alone a single query goes from 4m19s to 43s; the
smaller arms roughly halve.

**It is expensive in disk.** The index carries ~860 MB of fixed overhead *per arm*, so a
112 KB arm becomes 859 MB, and it grows faster than the database: the whole set is **55 GB
indexed against 4.6 GB without**. Budget accordingly, or pass `SKIP_INDEX=1` and accept
slower searches. `INDEX_MEMORY=4G` raises the memory each index build may use (default 2G);
an arm built by an earlier run without an index is indexed in place on the next run rather
than rebuilt.

Each arm comes from a different third-party host, and any of them can be down on the day.
A failed arm is logged and the run carries on with the rest; the summary at the end names
what is missing, prints the command to retry just those arms, and exits non-zero. Nothing
is left half-built: an arm that failed does not look finished to the next run.

Counts and sizes measured on a build from scratch on 2026-10-07; the upstream sources move,
so treat them as the order of magnitude rather than a contract.

| arm | entries | database | + index | source |
|---|---|---|---|---|
| `thpdb` | 162 | 112 KB | 859 MB | THPdb FDA-approved therapeutics (Figshare) |
| `therasabdab` | 2,533 | 884 KB | 862 MB | Thera-SAbDab therapeutic antibody chains |
| `plabdab_nano` | 2,306 | 724 KB | 862 MB | PLAbDab-nano VHH / VNAR / sdAb |
| `proteinbase_public` | 3,253 | 1.2 MB | 863 MB | published Proteinbase designs (public API) |
| `plabdab` | 350,350 | 168 MB | 1.4 GB | PLAbDab paired antibodies |
| `pdb` | 1,172,061 | 457 MB | 5.7 GB | PDB seqres chains |
| `swissprot` | 575,748 | 1.1 GB | 4.5 GB | UniProtKB/Swiss-Prot |
| `uspto` | 11,173,221 | 2.8 GB | 33 GB | USPTO patent sequences |
| **total** | | **4.6 GB** | **55 GB** | |

Everything except `uspto` builds in a couple of minutes on a fast connection — PDB and
SwissProt come down as prebuilt MMseqs2 databases rather than being built locally. `uspto`
is the one that costs: a ~3.4 GB download that unpacks to ~9.3 GB, then the conversion and
`createdb`. Measured end to end at **42 minutes**, plus another 10 for its index.
It is in the default arm list, so a plain run builds it; `SKIP_USPTO=1` leaves a set that
covers everything but patents (1.7 GB of database, 14 GB indexed).

`thpdb` comes from the [Figshare deposit](https://figshare.com/articles/dataset/5198005)
rather than the THPdb web host, which is frequently unreachable. That deposit is a TSV, and
multi-chain therapeutics pack every chain into one cell, so `scripts/build/thpdb_to_fasta.py`
splits them apart — 162 distinct chains across 163 therapeutics. The other 76 THPdb entries
carry no one-letter sequence upstream (`N.A.`, or three-letter notation).

Having the patent database on disk is not the same as searching it. The patent arm is off
the *query* path by default. On an indexed root a 4-sequence batch takes 2m17s with
`--patent` against ~70 s without — the patent arm is still 130 s of the 137. Unindexed the
same batch takes 13m17s, which is what the index is for.

Every arm is verified **by a real search, not by file presence**, before the script calls
it done. That check exists because `mmseqs createdb` reports success on a file it could not
parse: the EBI patent dump is an EMBL-style flat file with no FASTA parser, and ingesting
it directly yields a database holding exactly one record that silently matches nothing.

## The full screen

```bash
export PRESCREEN_DB_ROOT=~/prescreen-dbs
uv run prescreen examples/submissions.fasta -o report/
uv run prescreen submissions.csv --id-column id --sequence-column sequence -o report/
uv run prescreen submissions.fasta --patent -o report/        # + the USPTO arm
```

~72 s for 4 sequences across the seven public non-patent arms; ~17 s if only the four
small arms are built. The search is **batched** — one MMseqs2 call per database for the
whole input file — so database load dominates and N sequences cost about what one does.
Submit your whole batch at once rather than looping.

`PRESCREEN_DB_ROOT` takes a `:`-separated list, first match wins, so a big shared mount
and a few locally built arms can be mixed:

```bash
export PRESCREEN_DB_ROOT="$HOME/prescreen-dbs:/mnt/shared/sequence-dbs"
```

### Check what was actually searched

**This is the one that bites.** An arm whose database is not on disk is *skipped
silently* — a missing arm reads exactly like a sequence with no prior art, which is the
answer a prescreen must never give by accident. Nothing warns you on stdout, but every run
records the truth in `report.json`:

```python
import json
pa = json.load(open("report/report.json"))["results"]["my_seq"]["prior_art"]
pa["arms_searched"]   # ['plabdab_nano', 'therasabdab', 'thpdb', 'proteinbase_public']
pa["arms_missing"]    # ['pdb', 'swissprot', 'plabdab'] — the verdict is weaker than it looks
```

or before you run:

```bash
uv run python -c "from prescreen import refdb; print(refdb.available())"
```

## Output

`report.csv` is one flat row per submission; `report.json` has every hit, every threshold
and every intermediate. The **verdict** is the highest-ranked flag that fired; every flag
that fired is listed, and framework-reuse observations are reported separately as
`context` and never disqualify.

| verdict | meaning |
|---|---|
| `known_target_binder` | matches a known anti-TNF binder (whole sequence **and** paratope for antibody formats) |
| `target_region_match` | carries a known anti-TNF paratope on any framework |
| `verbatim_target_fragment` | is an exact fragment of a known binder |
| `existing_design` | matches an existing published design |
| `known_sequence` | matches a public sequence (PDB / SwissProt / PLAbDab / patents) |
| `target_binder_homolog` | a close relative of a known binder (non-antibody) |
| `near_known_sequence` | a near-match to a public sequence (non-antibody) |
| `batch_duplicate` | duplicates another submission in the same batch |
| `no_reference_signal` | too short to search reliably, and found nothing |
| `pass` | none of the above |

`target_binder_status` decides the message you send back: a hit to a `tested_non_binder`
means "this exact design was already tried and did not bind", a very different thing to
tell someone than a hit to a `confirmed_binder`.

`examples/prescreen_demo.ipynb` is a runnable walkthrough of all of this.

## Use it from Claude Code

The repository ships a skill at `.claude/skills/binder-prescreen/`, so cloning is the whole
install — Claude Code picks it up from the repo root. Ask for a screen in your own words, or
invoke it by name:

```
/binder-prescreen screen submissions.fasta
```

It covers the mechanics: installing, building the reference arms, running a batch in one
call rather than a loop, and checking which arms were actually searched before trusting a
clean result. What the verdicts mean stays here and in `PRESCREEN_SPEC.md`.

## Reference data

`prescreen/data/` ships the curated TNF-α set: **4,964** unique known or claimed binders
with per-entry provenance — 28 confirmed binders, 4,901 patent/document-level claims, 35
Proteinbase designs tested and found *not* to bind — plus **2,506** deduplicated anti-TNF
IMGT CDRs as the paratope reference. By primary source: USPTO patents (4,423), PLAbDab
incl. Thera-SAbDab (479), Proteinbase (35), PDB (25, all contact-verified against P01375),
UniProt (1), THPdb (1). Named agents are resolved for adalimumab, infliximab,
certolizumab, golimumab, etanercept, ozoralizumab and others.

Composition and the 15 numbered coverage gaps:
[`prescreen/data/tnfa_coverage_gaps.md`](prescreen/data/tnfa_coverage_gaps.md). Per-arm
counts, provenance and snapshot dates: `prescreen/data/reference_manifest.json`.

`--target` points the known-binder arm at a different curated set, so the filter works for
any target, not just TNF-α.

## License

[Elastic License 2.0](LICENSE). Source-available, not OSI open-source: you may use, copy,
modify and redistribute it, but not provide it to others as a hosted or managed service,
and not circumvent its license key functionality or remove licensing notices.
