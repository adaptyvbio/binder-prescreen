# Binder prescreen

A sequence-level prior-art and known-binder filter for protein-design binder
competitions at Adaptyv. 

Given a submitted sequence it answers two questions, from sequence alone:

1. **Is this already public?** MMseqs2 against PDB, SwissProt, PLAbDab (+ PLAbDab-nano
   for single-domain formats), Thera-SAbDab, THPdb, the published Proteinbase designs,
   and 10.2 M USPTO patent sequences.
2. **Is it a known binder of the target?** Whole-sequence *and* **binding-region**
   comparison against a curated set of known or claimed binders for that target. 

This pipeline does **not** judge whether the scaffold is de novo. Reusing a published
framework is allowed in these competitions, so the molecule category (nanobody, Fab, IgG,
affibody, monobody, DARPin, miniprotein, peptide) is used only to decide which binding
region to compare. Novelty will be computed on ProteinBase using the criteria outlined [here](https://adaptyvbio.com/blog/novelty).

## Why the binding region is compared separately

Antibodies are well known-scaffolds and in contrast to de novo miniproteins only a small fraction of their sequence must be designed to alter the binding target. These are the complimentary determining regions. This region is also often called the paratope. If we would use whole-sequence similarity it would screen out many submitted designs because only the paratope was changed. 

Comparing the paratope independently fixes both error directions:

| submission | whole sequence | paratope | verdict |
|---|---|---|---|
| known framework, **new CDR3** | ~90 % identical to a known binder | no match | **passes** — a new binder, allowed |
| **known CDR3**, fresh framework | looks unrelated | matches | **flagged** — a known binder in disguise |

The calibration behind the identification of the different molecule classes and detection of the paratopes is in [`PRESCREEN_SPEC.md`](PRESCREEN_SPEC.md); [`HOW_IT_WORKS.md`](HOW_IT_WORKS.md) is a plain-language walkthrough.

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

## Screen a design in 30 seconds, with no databases

`--skip-prior-art` turns off the public-corpus search and runs the target arm only:
whole-sequence against the curated binder set, the paratope comparison, and the verbatim
fragment check. This needs **no downloads**.

```bash
uv run prescreen --list-targets                                     # registered sets
uv run prescreen examples/submissions.fasta --skip-prior-art  --target tnfa -o report/
```

Replace `--target <slug>` to pick a target other than the default. 

`examples/` carries a worked set per target — `submissions.fasta` (TNF-α) and
`submissions_glp1r.fasta` (GLP-1R); `tfr1` has none yet. Each record's description line
states the construct and the verdict it should produce, so they double as a smoke test.

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
uv run prescreen examples/submissions.fasta --target tnfa -o report/
uv run prescreen examples/submissions.csv  --target glp1r -o report/
uv run prescreen examples/submissions.csv  --target tnfa --no-patent -o report/
```

### Choosing the target


| slug | target | UniProt | references | composition |
|---|---|---|---|---|
| `tnfa` | TNF-α | [P01375](https://www.uniprot.org/uniprotkb/P01375) | 4,964 | antibody-dominated — ~3,100 Ig chains and Fvs, plus affibodies, monobodies, miniproteins, receptor-derived decoys |
| `glp1r` | GLP-1 receptor | [P43220](https://www.uniprot.org/uniprotkb/P43220) | 225 | peptide-dominated — 215 GLP-1 analogue peptides, 7 antibody chains, 3 peptide-Fc fusions |
| `tfr1` | Transferrin receptor 1 (TfR1/CD71) | [P02786](https://www.uniprot.org/uniprotkb/P02786) | 6,695 | mixed — ~3,000 claimed peptides/loops, ~1,600 Ig chains and Fvs, 150 VNARs, 39 VHHs, FN3 domains, engineered Fc/CH3 patches, miniproteins, natural ligands |

`tnfa` is the default, so omitting `--target` screens against TNF-α. 


**The thresholds were calibrated on `tnfa`, and the whole-sequence ones re-checked on
`glp1r`** (442 fixtures: 50/50 knowns caught, no false positives on 80 unrelated designs
or 6 framework-reuse chimeras), so they hold for both. The **binding-region** cut point is
still TNF-only on `glp1r`: it has 7 antibody chains, so barely any fixture exercises it.

**`tfr1` re-checked both** (927 fixtures, all eight arms): 103/105 verbatim knowns caught,
no target-arm flag on 80 unrelated designs or 358 germline light chains, and both the
whole-sequence and the binding-region cut points sit mid-plateau, so the inherited values
hold there too.

Two things to expect on `glp1r` specifically. Some of its peptides carry `X` from ChEMBL's
non-standard residues and match only where the wildcard lands. And short submissions come
back `no_reference_signal` rather than `pass` — below 25 aa the filter refuses to call
silence a pass, because a short sequence finding nothing may be a true negative or a blind
search.

Three things to expect on `tfr1`. Its 203 references under 6 aa cannot be matched by any
arm. A binder fused to a wild-type Fc is judged on the binder: the Fc's own hits are
reported as `known_carrier_prior_art` / `known_carrier_new_paratope` context, while an
engineered TfR-binding Fc is compared on its CH3 patch. And a lightly mutated copy of one
of its ~500 light-chain-only patent claims is caught less often than on `tnfa` (77-78 % at
2-5 % mutation, against 95-96 %), because a germline-common CDR-L3 cannot decide a call;
a verbatim copy still can (`verbatim_copy`).

From Python, the slug goes on the config:

```python
from prescreen import Config, screen
cfg = Config.from_env()
cfg.target = "glp1r"          # or leave it None for the default
out = screen({"sub_001": "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS"}, cfg=cfg)
print(out["target"])          # what it actually compared against
```

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
  that is not the matching chain. 
  A chain that carries more than one variable domain — a bivalent or trivalent VHH, a
  dual-variable heavy chain, an scFv or a Fab written as one string — is split again into
  those domains, so a known binder sitting next to another domain is still found. Results
  are aggregated back to one row per submission; `prior_art_chain` / `target_chain` and
  `prior_art_domain` / `target_domain` say which input chain and which domain produced
  each hit.
- `molecule_class` — optional


`PRESCREEN_DB_ROOT` takes a `:`-separated list, first match wins, so a big shared mount
and a few locally built arms can be mixed:

```bash
export PRESCREEN_DB_ROOT="$HOME/prescreen-dbs:/mnt/shared/sequence-dbs"
```

### Check what was actually searched


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

The repository ships a skill at `.claude/skills/binder-prescreen/`

Ask for a screen in your own words, or invoke it by name:

```
/binder-prescreen screen submissions.fasta
```

## Reference data

Curated known-binder sets live one directory per target under `prescreen/data/targets/`

**`tnfa`** — **4,964** unique known or claimed binders with per-entry provenance: 28
confirmed binders, 4,901 patent/document-level claims, 35 Proteinbase designs tested and
found *not* to bind, plus **4,075** deduplicated anti-TNF IMGT CDRs. By primary source:
USPTO patents (4,423), PLAbDab incl. Thera-SAbDab (479), Proteinbase (35), PDB (25, all
contact-verified against P01375), UniProt (1), THPdb (1). Named agents are resolved for
adalimumab, infliximab, certolizumab, golimumab, etanercept, ozoralizumab and others.
Composition and the 15 numbered coverage gaps:
[`tnfa_coverage_gaps.md`](prescreen/data/targets/tnfa/tnfa_coverage_gaps.md).

**`glp1r`** — **225** references: 215 GLP-1 analogue peptides (ChEMBL), 7 antibody chains
and 3 peptide-Fc fusions, with **218** CDR-table rows of which 202 are the bare peptides
carried as chain-agnostic wildcard entries. By primary source: ChEMBL (192), PDB (33).
Named agents include dulaglutide, survodutide, efinopegdutide, efsubaglutide alfa and
maridebart cafraglutide. Smaller and peptide-dominated where `tnfa` is
antibody-dominated, so it inherits the TNF-calibrated thresholds untested and carries the
two coverage limits noted above.

**`tfr1`** — **6,695** references covering the formats TfR1 is attacked with: Ig V
domains, Fvs and Fabs, VNARs, VHHs, engineered Fc/CH3 patches, Centyrin-type FN3 domains,
cystine-dense and bicyclic peptides, designed miniproteins and natural ligands, plus
**3,326** CDR-table rows annotated against PLAbDab (CDR-L3 with `--near-identity 0.80`).
53 are confirmed binders (structure, reported assay or an INN), the rest disclosed claims;
there are no tested non-binders. By primary source: JPO (3,705), USPTO (2,483), EPO (178),
KIPO (157), PLAbDab (146), PDB (17), literature (8), NCBI patents (1). Named agents include
trontinemab, pabinafusp, lepunafusp, posnafusp, delpacibart and zeleciment. One KEGG
DRUG-only reference is omitted from this copy, because KEGG may not be redistributed.


### Adding a target

A target is a directory under `prescreen/data/targets/` plus one entry in `targets.json`.

**1. The curated FASTA**, one record per known or claimed binder, provenance in the
header, pipe-separated and positional:

```
>MYT00001|class|binder_status|source|accession|named_agent
```

`class` must come from:
`fv_domain_vh`, `fv_domain_vl`, `igg_chain_heavy`, `igg_chain_light`, `igg_fv`, `fab`,
`vhh`, `vhh_putative`, `vnar`, `peptide_or_cdr` (the CDR-indexed ones), and `affibody`,
`monobody`, `miniprotein`, `receptor_derived`, `viral_decoy`, `designed_other`,
`scaffold_other`, `other`.

**2. Normalise it** — a curation track's own vocabulary rarely matches:

```bash
uv run python scripts/build/target_set_from_fasta.py \
    --fasta curated.fasta --slug mytarget \
    --out-dir prescreen/data/targets/mytarget
```

**3. Build the CDR index**, the reference set for the binding-region arm:

```bash
uv run python scripts/build/augment_binder_cdrs.py \
    --binders prescreen/data/targets/mytarget/mytarget_binders.fasta \
    --metadata prescreen/data/targets/mytarget/mytarget_binders.csv \
    --cdrs    prescreen/data/targets/mytarget/mytarget_binder_cdrs.csv \
    --prefix MYTCDRX --numbered-classes igg_chain_heavy,igg_chain_light,other
```


**3b. Annotate it against a general antibody corpus.** Without this the region arm treats
a germline CDR as evidence — and germline CDRs are the ones that recur across unrelated
antibodies:

```bash
uv run python scripts/build/annotate_cdr_corpus.py \
    --cdrs prescreen/data/targets/mytarget/mytarget_binder_cdrs.csv \
    --corpus /path/to/plabdab.fasta --jobs 24
```

This adds a `corpus_occurrences` column. A region above `region_corpus_max` (default 100)
is still reported but cannot decide a verdict on its own; a region belonging to a *named
or confirmed* binder is exempt, since it is only common because the corpus contains that
binder. An index without the column behaves as before, so this is additive.

Add `--near-identity 0.80` when the set indexes light chains. The region arm matches at
0.80 identity, so a CDR-L3 one or two edits from a germline ending is as germline as the
ending itself but rarely occurs *exactly*; this counts corpus CDR-L3s within that identity
too. Without it the TfR1 set called 345 of 358 PLAbDab light chains ending `QQRSNWPLT`
known binders, through a patent light chain ending `QQRSNWPPIT` (30 exact occurrences).

**4. Register it** in `prescreen/data/targets/targets.json`:

```json
"mytarget": {
  "name": "My target", "uniprot": "P00000", "dir": "mytarget",
  "fasta": "mytarget_binders.fasta", "metadata": "mytarget_binders.csv",
  "cdr_csv": "mytarget_binder_cdrs.csv",
  "notes": "What is in the set, and what is not."
}
```

## License

[Elastic License 2.0](LICENSE). Source-available, not OSI open-source: you may use, copy,
modify and redistribute it, but not provide it to others as a hosted or managed service,
and not circumvent its license key functionality or remove licensing notices.
