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

0. **Split into search units** (`chains.py`) — a multi-chain submission joins its chains
   with `:`, and a chain carrying two or more variable domains is split again into those
   domains. Steps 1-2 and 5 read the whole submission **concatenated**, so a paired H+L
   types as a Fab; steps 3-4 work **per unit**, because every reference is one domain-sized
   chain and a concatenated query dilutes `qcov` by the fraction of the molecule that is not
   the matching part (a verbatim Fab scores 0.51 against its own heavy chain, under every
   cut in §5). Domains are found by recursing `variable_domain_spans` over the segments it
   has not covered — the library's own `split_variable_domains` is keyed by chain type and
   so caps at two, which loses the third domain of a trivalent VHH — and are kept only if
   they clear `ANTIBODY_GERMLINE_FLOOR`, which drops constant domains numbered as variable
   ones. A chain with fewer than two surviving domains is left whole. Hits are aggregated
   back to one record per submission, each carrying the chain and domain that found it.

   In the curated set, 151 of 4,964 references are genuinely multi-domain: 108 dual-variable
   heavy chains, 19 dual-variable light chains and 24 multivalent VHHs, ozoralizumab among
   them.

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

Everything is MMseqs2 for search and Levenshtein identity for short-region comparison.

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

**One exception, and why.** The CDR reference index is built only from the classes that
contributed CDRs — `fv_domain_vh/vl`, `vhh`, `vhh_putative`, `igg_chain_heavy/light`, `fab`,
`igg_fv`. A curated binder of any other class (`vnar`, `other`, `designed_other`,
`receptor_derived`) has **no CDRs in it at all**, so the region arm cannot corroborate a
match against such a reference however good the whole-sequence match is. Requiring
corroboration there meant a byte-identical resubmission of, say, the curated VNAR scored
0.40 on the region and was never called a known binder. When the matched reference's class
is absent from the index, the whole-sequence rule therefore stands on its own, and the
evidence records `paratope_corroborated: false` so the finding is not mistaken for a
confirmed one.

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
| uspto | USPTO patent protein sequences (searched on every run) | 10,206,785 | public |


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

A labelled set of 1,026 sequences was built (`calibrate.py`) from a stratified sample of
the curated binders: verbatim knowns, a 2–35 % random-substitution mutant series, grafts
(a known CDR3 spliced onto a different framework of the same class), framework-reuse
chimeras (a known framework carrying a scrambled region), published designs against other
targets, and shuffled controls. All eight arms were searched, the patent arm included.
Pass rate by class (`calibration_table.csv`):

| class | n | passed | reading |
|-------|---|--------|---------|
| known binder | 118 | **0 %** | every verbatim known binder flagged |
| mutant 2 % / 5 % | 236 | 0 % / 1 % | lightly mutated knowns still caught |
| mutant 10 % / 20 % / 35 % | 354 | 9 % / 19 % / 92 % | graded — a heavily rewritten binder is a new molecule |
| **grafted known CDR3** | 60 | **3 %** | 97 % of grafts onto a new framework are caught — the case whole-sequence misses |
| **reused framework** | 60 | **58 %** | correctly passed; of the rest, 38 % flag only as in-batch duplicates and 3 % as known-binder false positives |
| unrelated design | 80 | 0 % | **not a usable negative** — see below |
| shuffled control | 118 | 98 % | composition-matched negative |

The binding-region arm is the separator: region identity to the nearest known binder has
median **1.00 for grafts** (known paratope present) versus **0.39 for reused frameworks**
and **0.27 for unrelated designs**. The 0.80 region cut sits cleanly in that gap — 93 % of
grafts fall above it, 3 % of framework-reuse and 2 % of unrelated cases do.

**The unrelated-design row measures nothing.** `calibrate.py` draws it from the published
Proteinbase designs, which is also the `proteinbase_public` arm, so all 80 match themselves
at identity 1.00 and flag as `existing_design`. That is the correct verdict — they *are*
published designs — but it makes the row a self-identity check rather than a false-positive
rate. A real false-positive estimate needs the sampled designs held out of the design arm,
or a negative set drawn from a corpus that is not searched. Until then, the shuffled control
is the only honest negative in the set.

The patent corpus does most of the work: 645 of the 1,026 sequences take their best public
hit from `uspto`, more than every other arm combined.

**Do not read small movements in the mutant rows across runs.** Until the seed fix in
`calibrate.py` (`_stable_seed`), the per-reference mutation seed came from `hash(ref_id)`,
which Python salts per process, so each run mutated different residues: two consecutive
calibrations of the same reference set differed in 522 of 1,026 rows and the mutant pass
rates moved by up to 5 points for that reason alone. The deterministic classes — known,
graft, framework reuse, shuffled — were reproducible throughout and are the ones to compare
across runs. The mutant series is now reproducible too, from the first run generated after
that fix.

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
   known".
2. **No structure arm.** A submission that copies a known binder's fold with a rewritten
   sequence (low sequence identity, same 3-D paratope) is not caught. 
3. **Unpublished competition submissions.** The current round's submissions are not yet
   public, so a participant resubmitting their own earlier unpublished design cannot be
   caught.
4. **Database staleness.** The databases might not contain existing binders.
5. **Short-sequence composite.** For a 10–15-mer the aligner trims mutated termini and can
   report `fident = 1.0` with `qcov < 1`; always threshold on the composite or on fident
   and qcov jointly, never on fident alone.
6. **Humanised VHHs** that read a human VH hallmark tetrad may classify as `vh_domain`
   rather than `nanobody`; both are searched against the same references, so the verdict is
   unaffected, but the reported category can differ.


