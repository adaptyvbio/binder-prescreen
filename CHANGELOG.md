# Changelog

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
