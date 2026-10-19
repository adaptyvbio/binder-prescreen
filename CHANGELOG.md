# Changelog

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
