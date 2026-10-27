# Binder prescreen — specification

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
   A TCR chain is split on TCR numbering instead (`tcr.tcr_domains`), never on antibody
   numbering, which would cut an alpha chain into "two domains" at its C-alpha.

   A chain with no variable domains but an IgG hinge-CH2-CH3 inside it, and at least
   `FC_CARRIER_MIN_REMAINDER` (20) residues outside it, is split into the fused binder(s)
   and the Fc (`fc_carrier_units`). Searched whole, a peptide-Fc or miniprotein-Fc fusion is
   mostly Fc and matches any Fc-bearing reference on carrier alone; split, the fused binder
   is judged on itself and the Fc on its engineered CH3 patch (step 4).

   In the curated set, 151 of 4,964 references are genuinely multi-domain: 108 dual-variable
   heavy chains, 19 dual-variable light chains and 24 multivalent VHHs, ozoralizumab among
   them.

1. **Classify** (`classify.py`) — molecule category from sequence: T-cell receptors
   first (`tcr` for a paired alpha/beta, `tcr_domain` for one V domain; `tcr.py`), when a
   V domain numbers against antpack's TCR germlines at ≥ 0.60 *and* better than against
   antibody ones — checked first because a TCR chain otherwise force-fits an antibody
   germline at 0.55-0.66 and types as an Fv; antibody formats
   (Fv / scFv / Fab / IgG / nanobody / VH-domain) via antpack IMGT numbering; alternative
   scaffolds (affibody / monobody / DARPin / Centyrin-type FN3 / engineered Fc,
   `fc_tfr`) via framework identity to the canonical reference, gated at 0.85 — `fc_tfr`
   only for a query no longer than 1.15× an Fc, so a whole IgG heavy chain stays an
   antibody; otherwise a length regime (peptide / miniprotein / small_protein). The category sets the `region_kind` that chooses the binding region.
2. **Extract the binding region** (`regions.py`, `cdrref.py`) — antibody CDRs (IMGT),
   TCR CDRs (IMGT, keyed `A`/`B` rather than `H`/`L`; CDR1/CDR2 are compared against the
   same TCR chain only), projected scaffold paratope, or the whole sequence for a plain
   peptide/miniprotein. Every CDR3 - antibody or TCR, either chain - is compared against
   every TCR CDR3: a TCR CDR3 lifted into an antibody, or onto the other TCR chain, is
   still the known paratope. On VDJdb 8 of 179k CDR3s of other specificities come within
   0.80 of a `kras-pmhc` reference on the opposite chain, against 4,348 on their own.
   An engineered Fc's paratope is the residues that *differ* from wild-type IgG1 at the
   TfR-binding CH3 patch positions, so a wild-type Fc projects to the empty paratope `-`.
3. **Prior-art arm** (`priorart.py`, `refdb.py`) — batched MMseqs2 against the public
   corpora, with a separate sensitive pass for queries < 50 aa, plus a paratope-aware
   re-check of the top hits (§3).
4. **Target-binder arm** (`target.py`) — whole-sequence search against the curated
   known-binder set, exact-substring containment, and binding-region comparison against
   the known paratopes (§3). Containment needs 8 literal residues, or 6 for a reference
   with no `X` that the antibody corpus has never seen (`corpus_occurrences` = 0), which
   is what lets a 7-aa TfR peptide such as T7 be found inside a larger design. A query
   unit typed `fc_tfr` has its whole-sequence hit re-checked on the CH3 patch: unless
   query and reference carry the same engineered patch (identity ≥ 0.8), the hit is
   explained by the carrier and is reported as `known_carrier_new_paratope` context
   instead of a known-binder call.
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
| 1 | `known_target_binder` | **non-antibody:** whole-seq composite ≥ `target_known`. **antibody/scaffold:** whole-seq fident ≥ `target_known_identity` AND qcov ≥ `target_known_coverage` AND the binding region matches | 0.90 / (0.90 & 0.80) |
| 2 | `target_region_match` | binding region (CDR3 / paratope) matches a known binder: exact, or Levenshtein ≥ `region_match` on a matching chain | exact, or 0.80 |
| 3 | `verbatim_target_fragment` | sequence is an exact substring of a known binder (or contains one). Whole-sequence categories only | ≥ 8 aa |
| 4 | `existing_design` | whole-seq ≥ `known_sequence` to the Proteinbase design corpus, with paratope also matching | 0.95 |
| 5 | `known_sequence` | whole-seq ≥ `known_sequence` to any public reference (PDB, SwissProt, PLAbDab, THPdb, patents), paratope also matching | 0.95 |
| 6 | `target_binder_homolog` | **non-antibody only:** whole-seq to a known binder in [`target_homolog`, `target_known`) | 0.70–0.90 |
| 7 | `near_known_sequence` | **non-antibody only:** prior-art similarity in [`near_known_sequence`, `known_sequence`) | 0.80–0.95 |
| 8 | `batch_duplicate` | clusters with another submission in the same batch | 0.95 id / 0.80 cov |
| 9 | `no_reference_signal` | no hit on either arm, at a length where k-mer search is unreliable (< 25 aa) | — |
| 10 | `pass` | none of the above | — |

**Index coverage is a correctness property, not a detail.** The region arm compares a
submission's CDR3 against this index and nothing else, so a binding region missing from it
cannot be recognised however exactly it is copied — the arm reports a low identity because
it has nothing to match, and that silence reads as "the paratope is new". The index
originally held only the conventional antibody classes, leaving 1,609 curated binding
regions unrepresented, 1,554 of them `peptide_or_cdr` entries that *are* claimed loops.
`scripts/build/augment_binder_cdrs.py` adds them: bare loops with no parent chain are
indexed as wildcard entries matching any region on any chain, which is the right semantics
for a claimed loop (it is prior art wherever it reappears), and length-normalised identity
keeps that honest since a 40-aa entry cannot score highly against a 14-aa CDR3.

Where a class still contributes no CDRs, a whole-sequence match to it cannot be
corroborated, so the known-binder rule stands on the whole-sequence evidence alone and the
evidence records `paratope_corroborated: false`.

A region carried by more than `region_corpus_max` (100) antibodies in a general corpus is
germline-common and cannot decide a call, however well it matches; the match is reported
as `target_region_match_not_decisive` context. For CDR-L3 the count can include corpus
endings within the region-match identity (`annotate_cdr_corpus.py --near-identity`), since
that is how the arm matches; the TfR1 index is annotated this way. One exception keeps a
resubmitted reference flagged: a whole-sequence hit at ≥ 0.99 identity and coverage whose
region matches *exactly* is `known_target_binder` with `verbatim_copy: true`. A
framework-reuse chimera cannot reach it, because a swapped CDR3 is never an exact match.

A TCR is judged on its *pair*. One TCR CDR3 carries a lot of germline - JDI's
`ALSGPSGAGSYQLT` is TRAV19 joined to the TRAJ28 segment `GAGSYQLT` - and patent listings
hold public CDR3s: US11458191 pair 16's beta `ASSLALSYEQY` is also an EBV EBNA4 TCR's.
Measured on VDJdb (2026-06-03) against the `kras-pmhc` set's 83 CDR3s: 2.4 % of 179k
CDR3s of other specificities are within 0.80 of one (16 exactly), and of 31,708 paired
TCRs 41 have both CDR3s within 0.80 of *some* reference but only 3 of the *same*
receptor - one of them TCR1 itself, listed again under an off-target peptide. So when the
best region match is a TCR CDR3 (`prescreen._tcr_pairing`):

- a paired query decides only if its two CDR3s are both within `region_match` of one
  reference receptor (`CdrReference.match_pair`, chains swapped allowed; the index's
  `receptor_ids` say which CDR3s belong together), or its exact match is a confirmed or
  named binder's;
- a lone CDR3 - a single TCR domain, or a TCR CDR3 carried by an antibody - decides only
  against a confirmed or named binder, or when it is the CDR3 of the very receptor its
  whole-sequence hit is (a lightly mutated patent V domain is a copy, not a lone CDR3).

Anything else is `target_region_match_not_decisive` context, with `pair_identity` and
`pair_receptor` where a pair was compared.

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

- **Antibody formats** → the submission's CDRs vs. 4,075 deduplicated known anti-TNF IMGT
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


### Known-binder set (`data/targets/<slug>/<slug>_binders.fasta`, `.csv`)

One per target; the figures below are the TNF-α set, which is also the set every cut
point in this document was measured on. A target registered later inherits those cut
points untested until `calibrate.py` is run over its own references.

`glp1r` has been through that run (2026-10-12, 442 fixtures): the **whole-sequence** cut
points reproduce their TNF behaviour on it — 50/50 knowns caught, no false positives on
80 unrelated designs or 6 framework-reuse chimeras, and `target_known` = 0.90 sitting
mid-plateau — so they are now measured on both sets. `region_match` is not: only ~10
fixtures and 6 grafts exercise the region arm there, because 215 of its 225 references
are peptides with no CDR to compare. That value remains TNF-only evidence.

`tfr1` has been through it too (2026-10-19, 927 fixtures, on the full 6,696-reference set
this copy is one KEGG reference short of): 103/105 verbatim knowns caught (the two misses
are 5-aa references), no target-arm flag on 80 unrelated designs or on 358 germline light
chains, `target_known` mid-plateau from 0.86 to 0.94 and `region_match` flat from 0.70 to
0.84 - so both are measured on this set as well. `known_target_binder` at 2/5/10/20/35 %
mutation is 78/77/51/0/0 % against TNF's 96/95/54/0/0: the set carries ~500
light-chain-only patent claims, and a lightly mutated one whose CDR-L3 is germline-common
cannot be told from a germline light chain.

`kras-pmhc` (2026-10-27, 144 references, all eight arms) is TCR-only. The standard run (909
fixtures) catches 101/101 knowns and 100/99/84/51/1 % at 2/5/10/20/35 % mutation, with no
target-arm flag on 80 unrelated designs, 101 shuffled controls or 61 framework-reuse
chimeras; 11 of its 61 grafts are flagged, the other 50 moving a patent claim's CDR3 into
a different single V domain, which the lone-CDR3 rule (§2) reports as context. 532 TCR
fixtures all behave as specified: CDR3 pairs grafted onto unrelated TCR frameworks, with
1-2 edits, swapped between chains, a confirmed binder's CDR3 as an antibody CDR-H3 and bare
CDR3 peptides are flagged; unrelated PDB TCRs, JDI's framework with foreign CDR3s, grafts
near on one chain, a patent CDR3 alone, and VDJdb TCRs of other specificities are not. The
inherited cut points are kept rather than re-measured; the TCR-specific evidence is the
VDJdb measurement behind the pair rule.


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
| `region_corpus_max` | 100 | a region carried by this many sequences in a general antibody corpus is germline, and cannot decide a match on its own (a *named or confirmed* binder's region is exempt) |
| `verbatim_fragment_min_length` | 20 | an exact copy of a *named or confirmed* binder at least this long inside a larger construct is flagged even for antibody formats, where a shorter containment is only shared framework |
| `target_known_identity` / `target_known_coverage` | 0.90 / 0.80 | antibody known-binder (identity AND coverage) |
| `target_known` | 0.90 | non-antibody known-binder (composite) |
| `region_match` | 0.80 | binding-region (CDR3 / paratope) match |
| `target_homolog` | 0.70 | non-antibody homolog band lower edge |
| `known_sequence` | 0.95 | public prior-art / existing-design |
| `near_known_sequence` | 0.80 | prior-art near-miss band lower edge |
| `batch_cluster_identity` / `_coverage` | 0.95 / 0.80 | within-batch duplicate |

### Calibration

A labelled set of 1,028 sequences was built (`calibrate.py`) from a stratified sample of
the curated binders: verbatim knowns, a 2–35 % random-substitution mutant series, grafts
(a known CDR3 spliced onto a different framework of the same class), framework-reuse
chimeras (a known framework carrying a scrambled region), published designs against other
targets, and shuffled controls. All eight arms were searched, the patent arm included.
Pass rate by class (`calibration_table.csv`):

| class | n | passed | reading |
|-------|---|--------|---------|
| known binder | 118 | **0 %** | every verbatim known binder flagged |
| mutant 2 % / 5 % | 236 | 0 % / 1 % | lightly mutated knowns still caught |
| mutant 10 % / 20 % / 35 % | 354 | 12 % / 30 % / 96 % | graded — a heavily rewritten binder is a new molecule |
| **grafted known CDR3** | 61 | **3 %** | 59 of 61 grafts onto a new framework are caught — the case whole-sequence misses |
| **reused framework** | 61 | **57 %** | correctly passed; a further 33 % flag only as in-batch duplicates. The remaining 6 are a fixture limitation, not a filter error — see below |
| unrelated design | 80 | 0 % | **not a usable negative** — see below |
| shuffled control | 118 | 98 % | composition-matched negative |

The binding-region arm is the separator: region identity to the nearest known binder has
median **1.00 for grafts** (known paratope present) versus **0.44 for reused frameworks**
and **0.30 for unrelated designs**. The 0.80 region cut sits cleanly in that gap — 97 % of
grafts fall above it, 8 % of framework-reuse and 1 % of unrelated cases do.

**The five flagged framework-reuse rows are a fixture limitation.** Every one of them has
`fident = 1.00`, `qcov = 1.00` and region identity 1.00 against a curated binder, and every
one is built on a multi-domain host — three bivalent, one trivalent, one carrying four
variable domains. `fixtures.focus_span` returns only the H-chain CDR3, so the chimera scrambles
one paratope and leaves the others verbatim; the construct still *contains* an intact known
binder and the filter is right to say so. They are the same 8 % that sit above the 0.80
region cut. Until the fixture scrambles every domain, this class carries no false-positive
rate — the number it would yield is a property of the fixture, not of the filter.

**The unrelated-design row measures nothing.** `calibrate.py` draws it from the published
Proteinbase designs, which is also the `proteinbase_public` arm, so all 80 match themselves
at identity 1.00 and flag as `existing_design`. That is the correct verdict — they *are*
published designs — but it makes the row a self-identity check rather than a false-positive
rate. A real false-positive estimate needs the sampled designs held out of the design arm,
or a negative set drawn from a corpus that is not searched. Until then, the shuffled control
is the only honest negative in the set.

The patent corpus does most of the work: 630 of the 1,026 sequences take their best public
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
7. **TCR references are mostly CDR3s, and a lone CDR3 is weak evidence.** Receptors known
   only by their CDR3s (IEDB, VDJdb, a few patents) are found through them, not on whole V
   domain. A single TCR domain, or an antibody carrying a TCR CDR3, is flagged only
   against a confirmed or named binder's CDR3: a patent-claimed CDR3 alone is context,
   because patents hold public CDR3s shared with TCRs of other specificities.


