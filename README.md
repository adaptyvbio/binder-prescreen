# TNF-α binder prescreen

A sequence-level prior-art and known-binder filter for protein-design binder
competitions. 

Given a submitted sequence it answers two questions, from sequence alone:

1. **Is this already public?** MMseqs2 against PDB, SwissProt, PLAbDab (+ PLAbDab-nano
   for single-domain formats), Thera-SAbDab, THPdb, the published Proteinbase designs,
   and 10.2 M USPTO patent sequences.
2. **Is it a known TNF-α binder?** Whole-sequence *and* **binding-region** comparison
   against a curated set of 4,964 known or claimed anti-TNF binders that ships with the
   package.

It deliberately does **not** judge whether the scaffold is de novo. Reusing a published
framework is allowed in these competitions, so the molecule category (nanobody, Fab, IgG,
affibody, monobody, DARPin, miniprotein, peptide) is used only to decide which binding
region to compare. Novelty will be computed on ProteinBase using the criteria outlined [here](https://adaptyvbio.com/blog/novelty).

## Why the binding region is compared separately

Antibodies are most of the known anti-TNF corpus and nearly all share the same germline
framework. Measured on the reference set, a randomly sampled VHH matches **1,843 of ~5,000
entries** at ≥0.80 coverage and 647 of them at ≥0.70 identity. For antibodies and other molecules with designable paratope and rigid framework we need a different search. 

Comparing the paratope independently fixes both error directions:

| submission | whole sequence | paratope | verdict |
|---|---|---|---|
| known framework, **new CDR3** | ~90 % identical to a known binder | no match | **passes** — a new binder, allowed |
| **known CDR3**, fresh framework | looks unrelated | matches | **flagged** — a known binder in disguise |

That discrimination is the whole point of the tool. Cut points, the calibration behind
them (1,026 labelled sequences), and the deliberate blind spots are in
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
arm**. `uspto` is by far the largest arm (~42 min to build, 1.7 GB of the 4.6 GB total);
`SKIP_USPTO=1` leaves it out, at the cost of a screen that cannot see patent claims. Each arm lands at `<root>/<arm>/<arm>` as an MMseqs2 database, which is the layout
`prescreen.refdb.resolve()` expects.

Every arm is also **indexed** (`mmseqs createindex`). This is much faster but requires a significant amount of disk space.

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


## The full screen

```bash
export PRESCREEN_DB_ROOT=~/prescreen-dbs
uv run prescreen examples/submissions.fasta -o report/
uv run prescreen examples/submissions.csv -o report/
uv run prescreen examples/submissions.csv --no-patent -o report/   # without the USPTO arm
```

Every arm, the USPTO patent arm included, is searched on every run: a binder claimed in a
granted patent and deposited nowhere else is invisible without it. It is the slow arm
(10.2 M sequences, ~6–7 min per batch), and the cost is amortised over the whole batch —
which is the reason to submit a batch rather than one sequence at a time. `--no-patent`
drops it.

### Input format

CSV follows the [Proteinbase competition submission
template](https://proteinbase.com/templates/competition-submission-template.csv):

```csv
name,sequence,molecule_class
my-nanobody,QVQLVESGGGLVQAGGSLRLSCAAS...,nanobody
my-fab,EVQLVESGGG...VTVSS:DIQMTQSP...KVEIK,fab_kappa
```

- `name` — the submission id. Duplicates are rejected rather than silently overwritten.
- `sequence` — chains of a multi-chain entry are joined by `:`, and the two halves of the
  pipeline treat them differently on purpose. **Classification** reads the chains
  concatenated, with no linker inserted, so a Fab still types as a Fab rather than as two
  loose chains. **The searches run per chain**, because every reference is a single chain
  (median 108 aa) and a concatenated query dilutes coverage by the fraction of the molecule
  that is not the matching chain — a verbatim adalimumab Fab queried as one 432-aa string
  scores `qcov = 0.51` against its own heavy chain, under every threshold that matters.
  A chain that carries more than one variable domain — a bivalent or trivalent VHH, a
  dual-variable heavy chain, an scFv or a Fab written as one string — is split again into
  those domains, so a known binder sitting next to another domain is still found. Results
  are aggregated back to one row per submission; `prior_art_chain` / `target_chain` and
  `prior_art_domain` / `target_domain` say which input chain and which domain produced
  each hit.
- `molecule_class` — optional, and **cross-checked, never trusted**: the classifier still
  decides which region is compared, and a disagreement is reported as
  `declared_class_match=mismatch` instead of changing the comparison. A declaration the
  classifier does not know is reported as `unknown`.

Plain FASTA still works, and `--id-column` / `--sequence-column` / `--class-column`
override the column names for a CSV written to a different convention.


`PRESCREEN_DB_ROOT` takes a `:`-separated list, first match wins, so a big shared mount
and a few locally built arms can be mixed:

```bash
export PRESCREEN_DB_ROOT="$HOME/prescreen-dbs:/mnt/shared/sequence-dbs"
```

### Check what was actually searched

A missing arm reads exactly like a sequence with no prior art, which is the answer a
prescreen must never give by accident — so **a run whose reference databases are
incomplete aborts** rather than screening against whatever is mounted. `--allow-missing-arms`
overrides that when an incomplete screen is genuinely what you want; the run then warns on
stderr, and `prior_art_arms_missing` in `report.csv` names the arms left out on every row.
Either way, `report.json` records exactly what was searched:

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
Proteinbase designs tested and found *not* to bind — plus **4,075** deduplicated anti-TNF
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
