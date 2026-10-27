# Changelog

## 0.5.0 - 2026-10-27

### Added
- **A KRAS-G12D / HLA-A*11:01 pMHC known-binder set** (`--target kras-pmhc`): 51 TCRs
  from IEDB, VDJdb, the PDB and patent sequence listings (US11458191/US11382954,
  WO2016085904, WO2019112941, US12624084), as 61 V domains and 83 CDR3s. Binders shown
  only for G12V, and patent TCRs claimed for G12V/G12D, are kept as claims. Calibrated
  with all eight arms over 909 fixtures and 532 TCR-specific ones: every known, every
  CDR3-pair graft and 1-2-edit variant flagged, and no target-arm flag on unrelated
  designs, unrelated PDB TCRs, framework-reuse chimeras or VDJdb TCRs of other
  specificities.
- **T-cell receptors are a category** (`tcr`, `tcr_domain`). A V domain that numbers
  against antpack's TCR germlines at >= 0.60 and better than against antibody ones is a
  TCR; its CDRs are keyed `A`/`B` and compared against TCR CDRs. A single-chain
  Valpha-linker-Vbeta is split into its two domains for the searches, as an scFv is.
- **A TCR is judged on its pair.** A TCR CDR3 match decides only if both of the query's
  CDR3s match one reference receptor (`receptor_ids` in the CDR index), or the exact
  match is a confirmed or named binder's; a lone CDR3 decides only against a confirmed or
  named binder, or as the CDR3 of the V domain the query copies. Otherwise it is
  `target_region_match_not_decisive` context. On VDJdb, 2.4 % of unrelated CDR3s come
  within 0.80 of a reference CDR3 - patents hold public ones - but only 3 of 31,708
  paired TCRs match one receptor on both chains.
- **TCR CDR3s are found wherever they are put.** Every CDR3 - an antibody's, or either TCR
  chain's - is compared against the CDR3s of both TCR chains, so a TCR CDR3 used as a
  CDR-H3 (a TCR-mimic) or moved to the other chain is still the known paratope. No other
  target's CDR index holds TCR rows, so nothing changes for them.
- `scripts/build/tcr_prior_art_to_target.py` builds a TCR target set from a merged
  prior-art table; `examples/submissions_kras-pmhc.fasta` is a worked set.

### Fixed
- **A TCR was typed as an antibody.** antpack fits a TCR chain to antibody germlines at
  0.55-0.66, above the antibody floor, so a beta chain with its constant domain was an
  `fv` and an alpha chain a `single_domain_antibody`, with CDRs compared as H and L. No
  reference in the TNF-alpha, GLP-1R or TfR1 sets changes category (11,885 checked), and
  their region caches are regenerated with identical content.

## 0.4.1 - 2026-10-19

### Fixed
- **Every wild-type-Fc fusion was called a known sequence.** 0.4.0 splits an Fc off its
  fused binder and searches the parts separately, and a wild-type IgG1 Fc is in PDB
  verbatim, so the prior-art arms called any new peptide-, miniprotein- or VHH-Fc fusion
  `known_sequence` on the carrier alone. A wild-type Fc carrier's prior-art hits are now
  left out of the call and reported as `known_carrier_prior_art` context; an engineered
  Fc keeps them, since its CH3 patch is the binder, and a bare Fc submitted on its own
  is still a known sequence.

## 0.4.0 - 2026-10-19

### Added
- **A TfR1/CD71 known-binder set** (`--target tfr1`): 6,695 references across the formats
  the target is attacked with - Ig V domains, VNAR, VHH, engineered Fc/CH3 patches,
  Centyrin-type FN3, cystine-dense and bicyclic peptides, miniproteins and natural
  ligands. One KEGG DRUG-only reference is omitted from this copy and five others are
  attributed to their patent source instead, because KEGG may not be redistributed.
  Calibrated over 927 fixtures with all eight arms: 103/105 verbatim knowns caught, no
  target-arm flag on 80 unrelated designs or 358 germline light chains, and the inherited
  cut points sit mid-plateau, so they hold.
- **Centyrin-type FN3 domains and engineered Fc are scaffold categories** (`centyrin`,
  `fc_tfr`), each with a projected paratope. An engineered Fc's paratope is the residues
  that differ from wild-type IgG1 at the TfR-binding CH3 patch, so a wild-type Fc has none.
- **Very short references can be found inside a design.** Containment needed 8 residues,
  which put 7-aa TfR peptides such as T7 out of reach. A reference with no `X` that the
  antibody corpus has never seen now needs 6.

### Fixed
- **An Fc fusion could be called a known binder on its Fc alone.** A peptide-Fc or
  miniprotein-Fc construct is mostly Fc, and Fc is ~95% conserved, so whole-sequence it
  matched any Fc-bearing reference - an anti-TNF IgG heavy chain included - as a known
  binder. The Fc is now split off and searched on its own, and a bare-Fc hit only counts
  when both sides carry the same engineered CH3 patch; otherwise it is reported as
  `known_carrier_new_paratope` context. On the TNF-alpha set this clears a false
  `known_target_binder` on efinopegdutide, and removes a carrier-only
  `target_binder_homolog` flag from dulaglutide and efsubaglutide alfa.
- **A verbatim copy of a reference was missed when its CDR3 is germline-common.** Such a
  region cannot decide a known-binder call on its own, which is right for a framework that
  shares it and wrong for the reference itself. A whole-sequence hit at >= 0.99 identity
  and coverage whose region also matches exactly is now `known_target_binder`, reported
  with `verbatim_copy: true`. On the TNF-alpha set this newly flags one PLAbDab light
  chain that is in the set verbatim.

### Changed
- `annotate_cdr_corpus.py --near-identity 0.80` also counts corpus CDR-L3s within the
  region-match identity, not only exact hits. Run on the TfR1 set: without it, a patent
  light chain ending `QQRSNWPPIT` (30 exact occurrences, 0.80 from the germline
  `QQRSNWPLT`) called 345 of 358 PLAbDab light chains known TfR1 binders; with it, none.
  The TNF-alpha index is unchanged.

## 0.3.1 - 2026-10-15

### Fixed
- **Patent references written as Markush claims could never match.** A reference like
  `HXEGTXTSDXX`, where `X` means any residue, was compared literally, so the containment
  check could only find it in a sequence that contains the letter `X` - which no
  submission does. `X` in a reference is now a wildcard, as it already was for CDRs, and
  the report shows the residues that actually matched. A wildcard match counts only by
  the reference residues it agreed with literally, which must reach the same 8-residue
  floor as any reference - otherwise a Markush framework written as runs of `X` would
  match any short sequence whole.

## 0.3.0 - 2026-10-14

### Fixed
- **Designed antibodies with a shortened CDR2 were not recognised as antibodies.** One
  length band required every CDR to be at least 3 residues, so an Fv whose CDR2 had been
  truncated to 1-2 residues was typed as a small protein: the paratope arm never ran on
  it, and it was flagged instead on universal framework it happened to contain. The
  numbering itself was never the problem - these score 0.945-1.000 identity with all four
  conserved IMGT anchors. CDR1 and CDR2 length is now advisory, as it always should have
  been; the anchors and the identity floors are what exclude a non-antibody. This also
  corrects a shark VNAR in the curated set, a class that has no CDR2 by construction.
- **Framework fragments in the curated set could flag any antibody.** USPTO sequence
  listings enumerate framework segments as separate entries, so the TNF-alpha set held
  `GTLVTVSS` (heavy-chain FR4, present in 25% of all antibodies), VH FR1 and VL FR3 as if
  they were binders. References that common no longer support a containment match, unless
  they belong to a named or confirmed binder.

### Changed
- A sequence that could not be numbered now says why, in `report.csv` and in the record's
  context, instead of silently becoming a small protein. A category decided on length
  alone no longer reports confidence 1.0.

## 0.2.1 - 2026-10-14

### Fixed
- **A known binder fused to a carrier domain escaped the screen.** A verbatim GLP-1
  agonist joined to a nanobody by a linker passed both arms: the construct types as an
  antibody, so the paratope check compared the nanobody's CDRs and vetoed the
  known-binder call, while the agonist sat inside it already found and recorded as
  context. An exact copy of a named or confirmed binder of at least 20 residues is now
  flagged wherever it is carried.

## 0.2.0 - 2026-10-12

### Fixed
- **A germline light chain could fail a design on its own.** A Fab whose every
  heavy-chain CDR had been replaced by poly-glycine was still reported as a known
  TNF-α binder, because its ordinary kappa CDR-L3 exactly matched a patent-claimed
  region and outscored the chain that actually carried the design. Region matches are
  now weighed by how common the matched region is across a general antibody corpus, so
  a germline-encoded CDR is reported but cannot decide a verdict on its own.

### Added
- **GLP-1 receptor target set** (`--target glp1r`): 225 curated references, selectable
  alongside TNF-α from a registry of known-binder sets. Its whole-sequence cut points
  were calibrated on it; the binding-region cut point is still inherited from TNF-α.
