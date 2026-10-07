# TNF-α binder prescreen — specification

A sequence-level prescreen for submissions to the Anthropic × Adaptyv TNF-α binder
competition. It runs at submission time, from sequence alone, and answers two questions:

1. **Is this sequence already public prior art?** (known natural protein, patented
   sequence, or an existing deposited/published design)
2. **Is it a known binder of TNF-α, or does it carry a known anti-TNF binding region?**

It is **not** a scaffold-novelty filter. The companion novelty scales judge whether a
*scaffold* is de novo; this filter deliberately does the opposite, because
reusing a published framework is permitted in the competition. See §5 for the blind spots
that follow from that choice — above all, a novel scaffold carrying a genuinely new
paratope against the same epitope passes by design.

---

## 1. Pipeline

For each submission, in order (`prescreen.screen`):

1. **Classify** (`classify.py`) — molecule category from sequence: antibody formats
   (Fv / scFv / Fab / IgG / nanobody / VH-domain) via antpack IMGT numbering; alternative
   scaffolds (affibody / monobody / DARPin) via framework identity to the canonical
   reference, gated at 0.85; otherwise a length regime (peptide / miniprotein /
   small_protein). The category sets the `region_kind` that chooses the binding region.
2. **Extract the binding region** (`regions.py`, `cdrref.py`) — antibody CDRs (IMGT),
   projected scaffold paratope, or the whole sequence for a plain peptide/miniprotein.
3. **Prior-art arm** (`priorart.py`, `refdb.py`) — batched MMseqs2 against the public
   corpora, with a separate sensitive pass for queries < 50 aa, plus a paratope-aware
   re-check of the top hits (§3).
4. **Target-binder arm** (`target.py`) — whole-sequence search against the curated
   known-binder set, exact-substring containment, and binding-region comparison against
   the known paratopes (§3).
5. **Within-batch clustering** (`cluster.py`) — flags duplicate submissions in the batch.
6. **Verdict** (`flags.py`) — the ordered flag logic in §2.

Everything is MMseqs2 for search and Levenshtein identity for short-region comparison; no
HMMER, no profile-HMMs, no structure.

---

## 2. Flags and verdicts

Every flag that fires is reported; the **verdict** is the highest-ranked one
(`FLAG_ORDER`). Framework-reuse observations are reported separately as `context` and
never set the verdict — reusing a published scaffold is allowed.

Rank order (most to least disqualifying):

| # | flag | fires when | cut point |
|---|------|-----------|-----------|
| 1 | `known_target_binder` | **non-antibody:** whole-seq composite ≥ `tnf_known`. **antibody/scaffold:** whole-seq fident ≥ `tnf_known_identity` AND qcov ≥ `tnf_known_coverage` AND the binding region matches | 0.90 / (0.90 & 0.80) |
| 2 | `target_region_match` | binding region (CDR3 / paratope) matches a known binder: exact, or Levenshtein ≥ `region_match` on a matching chain | exact, or 0.80 |
| 3 | `verbatim_target_fragment` | sequence is an exact substring of a known binder (or contains one). Whole-sequence categories only | ≥ 8 aa |
| 4 | `existing_design` | whole-seq ≥ `known_sequence` to the Proteinbase design corpus, with paratope also matching | 0.95 |
| 5 | `known_sequence` | whole-seq ≥ `known_sequence` to any public reference (PDB, SwissProt, PLAbDab, THPdb, patents), paratope also matching | 0.95 |
| 6 | `target_binder_homolog` | **non-antibody only:** whole-seq to a known binder in [`tnf_homolog`, `tnf_known`) | 0.70–0.90 |
| 7 | `near_known_sequence` | **non-antibody only:** prior-art similarity in [`near_known_sequence`, `known_sequence`) | 0.80–0.95 |
| 8 | `batch_duplicate` | clusters with another submission in the same batch | 0.95 id / 0.80 cov |
| 9 | `no_reference_signal` | no hit on either arm, at a length where k-mer search is unreliable (< 25 aa) | — |
| 10 | `pass` | none of the above | — |

### The governing rule: antibody formats are judged on the paratope

For antibody and scaffold formats the framework is most of the sequence, and the
known-binder corpus cross-hits itself on framework alone (§3). So **whole-sequence
similarity never raises a flag for these categories by itself** — a known-binder or
known-sequence call additionally requires the binding region to match. The three
framework-only observations are recorded as `context`, not flags:

- `known_framework_new_paratope` — public framework, new binding region.
- `design_framework_reused` — same, matching the design corpus.
- `shared_framework_with_target_binder` — homolog-band similarity to a known binder with
  no region match.
- `framework_similar_to_public_sequence` / `framework_contains_reference_fragment` — a
  framework-level prior-art or substring hit on an antibody-format submission.

This is the case a whole-sequence filter gets wrong in both directions, and the reason the
paratope arm exists.

---

## 3. Why two arms, and why a paratope arm at all

The curated known-binder corpus (§4) is dominated by antibody-format sequences that share
germline frameworks. Measured on that corpus: a single probe VHH returns **1,843 hits at
≥ 0.80 coverage, 647 of them at ≥ 0.70 identity, all framework background — only 15 reach
≥ 0.90 identity.** A whole-sequence MMseqs2 flag at any identity below ~0.90 would
therefore mark almost every antibody-format submission as a known TNF binder.

The fix is to compare the **binding region** independently:

- **Antibody formats** → the submission's CDRs vs. 2,506 deduplicated known anti-TNF IMGT
  CDRs (`tnfa_binder_cdrs.csv`), matched by chain type, exact-first then Levenshtein.
  CDR1/CDR2 are germline-shared and only corroborate; the CDR3 is the determinant.
- **Alternative scaffolds** → projected paratope vs. same-class references.
- **Prior-art arm** → each top prior-art hit's own sequence (`tseq`) is re-numbered and
  its binding region compared with the query's, so a shared PLAbDab framework is reported
  as context, not as a known molecule.

Patent sequences use Markush-style claims, so some reference CDRs contain `X`; those are
excluded from identity arithmetic and matched only by a wildcard-aware exact test.

---

## 4. Reference sets

### Prior-art corpora (`reference_manifest.json`)

| arm | content | n | visibility |
|-----|---------|---|------------|
| pdb | PDB seqres protein chains | 987,345 | public |
| swissprot | SwissProt | 573,661 | public |
| plabdab | PLAbDab paired antibodies | 353,788 | public |
| plabdab_nano | VHH / VNAR / sdAb (single-domain — PLAbDab has none) | 2,228 | public |
| therasabdab | therapeutic antibody chains (incl. approved anti-TNFs) | 2,533 | public |
| thpdb | THPdb therapeutic proteins | 239 | public |
| proteinbase_public | published Proteinbase designs | 1,494 | public |
| uspto | USPTO patent protein sequences (off the submission path, §5) | 10,206,785 | public |
| proteinbase_internal | other entrants' unpublished submissions | 1,051 | **internal only** |

**Visibility is a correctness constraint.** `proteinbase_internal` and the past-submission
set may raise an internal duplicate flag for the organisers but must **never** appear as
evidence shown to a competitor. `best_prior_art` filters to public arms; organiser-only
searches require `Config(organiser_mode=True)`.

### Known-binder set (`tnfa_binders.fasta`, `.csv`)

4,964 unique known/claimed TNF-α binders with per-entry provenance: 28 confirmed binders,
4,901 patent/document-level claims, 35 tested non-binders (Proteinbase designs that did not
bind). By primary source (sums to the total; a few entries cite more than one):
USPTO patents (4,423), PLAbDab incl. Thera-SAbDab (479), Proteinbase (35), PDB (25, all
contact-verified against P01375), UniProt (1), THPdb (1). Named agents
resolved for adalimumab, infliximab, certolizumab, golimumab, etanercept, ozoralizumab and
others. 168 entries carry ambiguous (`X`) residues from Markush claims. Full composition
and the 15 numbered coverage gaps are in `tnfa_coverage_gaps.md`.

---

## 5. Cut points and calibration

Thresholds live in `config.py` and are overridable by `PRESCREEN_<FIELD>` environment
variables. Defaults:

| parameter | value | meaning |
|-----------|-------|---------|
| `tnf_known_identity` / `tnf_known_coverage` | 0.90 / 0.80 | antibody known-binder (identity AND coverage) |
| `tnf_known` | 0.90 | non-antibody known-binder (composite) |
| `region_match` | 0.80 | binding-region (CDR3 / paratope) match |
| `tnf_homolog` | 0.70 | non-antibody homolog band lower edge |
| `known_sequence` | 0.95 | public prior-art / existing-design |
| `near_known_sequence` | 0.80 | prior-art near-miss band lower edge |
| `batch_cluster_identity` / `_coverage` | 0.95 / 0.80 | within-batch duplicate |

### Calibration

A labelled set of 1,036 sequences was built (`calibrate.py`) from a stratified sample of
the curated binders: verbatim knowns, a 2–35 % random-substitution mutant series, grafts
(a known CDR3 spliced onto a different framework of the same class), framework-reuse
chimeras (a known framework carrying a scrambled region), published designs against other
targets, and shuffled controls. Pass rate by class (`fig_calibration.png`,
`calibration_table.csv`):

| class | n | passed | reading |
|-------|---|--------|---------|
| known binder | 118 | **0 %** | every verbatim known binder flagged |
| mutant 2 % / 5 % | 236 | 0 % / 1.7 % | lightly mutated knowns still caught |
| mutant 10 % / 20 % / 35 % | 354 | 12 % / 33 % / 91 % | graded — a heavily rewritten binder is a new molecule |
| **grafted known CDR3** | 65 | **6 %** | 94 % of grafts onto a new framework are caught — the case whole-sequence misses |
| **reused framework** | 65 | **55 %** | correctly passed; of the rest, 38 % flag only as in-batch duplicates and 6 % as known-binder false positives |
| unrelated design | 80 | 91 % | low false-positive rate on designs for other targets |
| shuffled control | 118 | 98 % | composition-matched negative |

The binding-region arm is the separator: region identity to the nearest known binder has
median **1.00 for grafts** (known paratope present) versus **0.40 for reused frameworks**
and **0.29 for unrelated designs**. The 0.80 region cut sits cleanly in that gap — 86 % of
grafts fall above it, ~5 % of framework-reuse and unrelated cases do.

### Evidence reported per flag

`arm`, target id, `fident`, `qcov`, `alnlen`, `qlen`, the `fident*qcov` composite, and
which length profile ran; for a region match, the reference CDR, chain, exact/identity, the
matched binder's `binder_status` and `named_agent`, and how many CDRs corroborated. A
`tested_non_binder` hit means "this exact design was already submitted and did not bind" —
a different message to the participant than a hit to a `confirmed_binder`. Snapshot date /
age per arm (from the manifest) belongs next to any "no prior art found" verdict.

---

## 6. Deliberate blind spots

1. **Scaffold de-novo-ness is not assessed.** A novel scaffold carrying a genuinely new
   paratope against the TNF-α epitope passes — by design. This filter answers "is it
   known", not "is it clever".
2. **No structure arm.** A submission that copies a known binder's fold with a rewritten
   sequence (low sequence identity, same 3-D paratope) is not caught. That is the
   scaffold novelty scale's job, not this one.
3. **Unpublished competition submissions.** The current round's submissions are not yet
   public, so a participant resubmitting their own earlier unpublished design cannot be
   caught except in organiser mode against the internal corpus.
4. **Patent arm is off the submission path.** The 10.2 M-sequence USPTO arm costs ~4 s per
   query and must run as an asynchronous batch (`Config(include_patent_arm=True)`); it also
   covers US filings only, last updated 2025-09.
5. **Database staleness.** The general and antibody arms on the mount are ~344 days old
   (built 2025-10-27). A binder published in the last year is prior art this set will not
   find — rerun `build_reference_dbs.sh` before the round opens.
6. **Short-sequence composite.** For a 10–15-mer the aligner trims mutated termini and can
   report `fident = 1.0` with `qcov < 1`; always threshold on the composite or on fident
   and qcov jointly, never on fident alone.
7. **Humanised VHHs** that read a human VH hallmark tetrad may classify as `vh_domain`
   rather than `nanobody`; both are searched against the same references, so the verdict is
   unaffected, but the reported category can differ.

---

## 7. Throughput

On the 1,036-sequence calibration batch (local, 8 threads, no patent arm): classification
13 s, target whole-sequence 5 s, target region 99 s, prior-art 142 s, clustering 3 s. The
per-arm MMseqs2 cost is dominated by PDB and SwissProt (~11–18 s per batch each) and is
near-constant in batch size, so throughput improves per-sequence as the batch grows. The
earlier multi-hour runs were resource contention with concurrent jobs, not the filter.
