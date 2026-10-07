# Known TNF-alpha binder reference set — coverage and gaps

Target: human TNF-alpha, UniProt **P01375** (soluble form = precursor residues 77–233).
Receptors used as receptor-derived references: TNFR1 **P19438**, TNFR2 **P20333**.

Purpose: ground truth for the "is this submission a known TNF-alpha binder, or does it carry a known
anti-TNF binding region on some other framework?" arm of the competition prescreen. This set deliberately
says nothing about whether a scaffold is de novo.

## What is in the set

`tnfa_binders.csv` / `tnfa_binders.fasta` — 4964 unique sequences (exact-duplicate collapsed, provenance merged).

| class | n |
|---|---|
| peptide_or_cdr | 1626 |
| fv_domain_vh | 1335 |
| fv_domain_vl | 1040 |
| vhh | 278 |
| igg_chain_heavy | 232 |
| igg_chain_light | 228 |
| affibody | 84 |
| other | 44 |
| receptor_derived | 33 |
| vhh_putative | 20 |
| miniprotein | 11 |
| monobody | 10 |
| designed_other | 9 |
| fab | 8 |
| igg_fv | 2 |
| vnar | 2 |
| scaffold_other | 1 |
| viral_decoy | 1 |

Binder status (how strong the TNF-binding evidence is):

| status | n | meaning |
|---|---|---|
| confirmed_binder | 28 | in a deposited complex with TNF-alpha (contact-verified), or a natural TNF receptor ectodomain / approved Fc-fusion |
| disclosed_binder_claim | 4901 | disclosed in a patent or publication whose subject is an anti-TNF-alpha binder; binding asserted at document level, not verified per sequence |
| tested_non_binder | 35 | previously submitted TNF-alpha-targeting designs on Proteinbase that were experimentally tested and did **not** bind — kept because they are the sequences most likely to be resubmitted |

Source databases (a sequence can have several):

| source | n sequences |
|---|---|
| USPTO patent sequences | 4620 |
| PLAbDab | 479 |
| Proteinbase | 35 |
| RCSB PDB | 25 |
| UniProtKB | 2 |
| THPdb | 1 |

`tnfa_binder_cdrs.csv` / `tnfa_binder_cdrs.fasta` — 4075 unique known anti-TNF binding regions:
2521 IMGT CDRs (AntPack numbering) extracted from the antibody-class sequences (1161 CDR3, 810 CDR1, 550 CDR2),
plus 1554 claimed loops from the `peptide_or_cdr`, `vnar`, `other` and `designed_other` references, which the
original index left out. The latter carry region `*` and chain `*` where no parent chain could be assigned, and
are matched against a query's region whatever chain it was numbered at — see
`scripts/build/augment_binder_cdrs.py`.
This is the file that answers "known anti-TNF binding region carried on any framework" — a submission can have a
completely foreign framework and still be caught by a CDR3 match.

## Verification actually performed

* Every PDB-derived chain was checked for real contact with a TNF chain: heavy-atom pairs within 5 A, computed from
  the deposited coordinates. Two candidates were **excluded** on this basis or on biology:
  21TV serum albumin (0 contacts — it is bound by the albumin-binding VHH arm of ozoralizumab, not by TNF) and
  4Y6O SIRT2 (a deacetylase bound to a myristoylated TNF-derived peptide substrate, not a TNF binder).
* Etanercept was identified rather than assumed: THPdb Th1005 (467 aa) is 100% identical to P20333 residues 23–257
  over full query coverage, and its C-terminal 232 aa are 99.1% identical to IgG1 heavy constant P01857.
* Cross-corroboration: the adalimumab heavy chain from 3WD5 is 99.0% identical to a USPTO patent heavy chain, and the
  TNF30 VHH from 8Z8M is 100% identical to a USPTO patent VHH — two independent sources agree.
* Negative control: streptococcal protein G B1 domain (1pga) returns zero hits against the packed DB.


## Class assignment method (so the next phase can re-derive it)

* PDB-, Proteinbase-, UniProt-, THPdb- and marker-asserted classes win: a chain curated as VHH/VNAR/Fab/affibody/
  monobody/receptor-derived keeps that label.
* Antibody classes otherwise come from AntPack IMGT numbering (V-domain identity >= 0.7 required).
* VHH vs VH is decided on the IMGT FR2 hallmark tetrad (positions 42, 49, 50, 52; human VH consensus V-G-L-W).
  >= 3 deviations -> `vhh`; exactly 2 -> `vhh_putative`; <= 1 -> `fv_domain_vh` / `igg_chain_heavy`. The tetrad is
  reported per entry in `vhh_hallmark_tetrad_imgt_42_49_50_52`. The split is sharply bimodal on this corpus
  (1498 sequences at VGLW, 275 at FERF/FGRF-type), which is why the rule is usable. Humanised nanobodies are the
  known failure mode: TNF30/ozoralizumab (8Z8M) and llama VHH2 (5M2J) both read VGLW and are only labelled `vhh`
  because the PDB assertion overrides the tetrad. Patent-derived humanised VHHs may therefore sit in
  `fv_domain_vh`.
* `affibody` = alignment to Staphylococcus protein A (P38507) over >= 40 aa; `monobody` = alignment to human
  fibronectin (P02751) over >= 70 aa; `receptor_derived` = >= 90% identity to a TNFR1/TNFR2 ectodomain over >= 60 aa.
* Sequences >= 90% identical to TNF-alpha itself over >= 80% query coverage (94 entries) and Ig-constant-only
  sequences with no V-domain (33 entries) were **removed** from the set: the first are the antigen, the second carry
  no paratope and would produce false flags.

## Gaps — binders known to exist whose sequences are not in this set

1. **Anthropic x Adaptyv 2026 competition submissions are not published yet.** The competition page lists no designs,
   and only 35 designs sit under Proteinbase target `human-tnfa` (all from one author, all BoltzGen, all
   experimentally non-binding). A participant resubmitting their own unpublished earlier submission cannot be caught
   by this file. This is the single largest coverage gap. The local Nov-2025 Proteinbase dumps contain all 35 but
   carry no target annotation, so the live API is the only route to the labels.
2. **Affilin (ubiquitin-based artificial binding protein) anti-TNF variant "10F"** — PMID 22363609,
   doi:10.1371/journal.pone.0031298. The open-access PDF describes the variant only as six randomised positions on
   human ubiquitin F45W plus a D58/Y59 insertion; no machine-readable sequence is printed and the complex is not
   deposited under P01375. Not reconstructed, by policy.
3. **Computationally designed anti-TNF peptides** — helical peptides targeting TNF-alpha (PMID 24038781,
   doi:10.1002/anie.201305963) and d-peptide TNF-alpha inhibitors (PMID 31102258, doi:10.1002/1873-3468.13444).
   Both closed access; sequences not retrievable from an allowlisted source.
4. **Published anti-TNF-alpha affibody clones** — PMID 19545237 (doi:10.1042/BA20090085), PMID 19576305, PMID
   33069194. Closed access. The affibody class is represented by 84 sequences from the USPTO patent
   "Staphylococcus protein A domain mutants that bind to TNF-alpha", but the specific published clones may differ
   from the patented ones.
5. **WP9QY-class TNFR1-mimetic peptides** were not individually identified. They may be among the 1626 `peptide_or_cdr` entries drawn from TNF-titled patents, but no entry is name-matched to them.
6. **Named clinical-stage agents without an individually name-matched sequence:** afelimomab, CDP571,
   placulumab (AT001/V565), DLX105/ESBA105, AZD9773, onercept, lenercept, pegsunercept. Several are very likely
   inside the ESBATech "Stable and soluble antibodies inhibiting TNF alpha" and TNFR patent families already in the
   set, and the TNFR-derived agents are covered in substance by the P19438/P20333 ectodomain entries — but the
   mapping from INN to sequence was not made.
7. **Thera-SAbDab was not accessed directly** (opig.stats.ox.ac.uk). The 7 Thera-SAbDab-derived rows embedded in the
   local PLAbDab snapshot (file dated 2024-08-29) served as the proxy: adalimumab, certolizumab, golimumab,
   infliximab, licaminlimab, remtolumab. Anti-TNF therapeutics that received an INN after Aug 2024 are missing.
   Ozoralizumab is present from the PDB rather than from Thera-SAbDab.
8. **PLAbDab assignment is document-level.** A patent titled "Anti-TNF-alpha antibodies" can also list control or
   unrelated antibodies, so a minority of the 479 PLAbDab-derived sequences may not be anti-TNF. Of 125 candidate
   documents, 86 were kept and 39 dropped: 35 as targeting something else (TL1A, OX40, GITR, CD40, 4-1BB/TNFRSF9,
   LT-beta-R, IL-12/IL-23, IL-17, GM-CSF, anti-idiotypic, anti-TNFR2) and 4 as too ambiguous to call
   ("Antibody formulations and methods of making same", "Molecules with extended half-lives",
   "Multivalent and monovalent multispecific complexes", "Human anti-self antibodies with high specificity from
   phage display libraries").
9. **Patent coverage is US-only** (EBI patentdata `uspto_prt`, LAST_UPDATE 2025-09). EP, WO, JP and CN sequence
   listings are not covered.
10. **168 entries contain ambiguous residues** (X/B/Z/J/U/O) from Markush-style patent claims; they are flagged in
    `has_ambiguous_residues` and should be excluded from identity arithmetic or treated as wildcards.
11. **44 entries are class `other`** — sequences from TNF-titled patents that are neither an Ig V-domain nor a
    recognised scaffold nor a short peptide. Some may be non-binder constructs.
12. **No de novo designed TNF-alpha binder exists in the Protein Design Archive.** The 1963-entry PDA snapshot
    contains exactly two TNF-related designs and both target receptors (9CU8 TNFR2_mb1 monobody, 9F7Z anti-4-1BB
    DARPin), not TNF-alpha. Likewise PubMed returns no anti-TNF-alpha DARPin. So the de novo arm of the prior art is
    close to empty — which is itself the finding, not a hole in the retrieval.
13. **TNFR-targeting agents are deliberately out of scope.** Anti-TNFR1/anti-TNFR2 antibodies, the TNFR2_mb1
    monobody, and agonistic TNFRSF binders were excluded. A submission that binds the receptor rather than the
    cytokine will not be flagged by this set.
14. **Non-protein binders out of scope**: small molecules (e.g. the 2AZ5/6X8x/9OJx fragment series co-crystallised
    with TNF) and aptamers.
15. One upstream annotation error is carried with a note: PDB 5M2J entity 2 is described as "Anti-(ED-B) scFV" but
    the entry title and citation identify it as llama VHH2, and it contacts TNF. It is included as a VHH.

## Honest assessment of completeness

* **Approved and clinical antibody biologics: essentially complete.** All five approved anti-TNF agents are present
  (adalimumab, infliximab, certolizumab pegol, golimumab, etanercept), plus ozoralizumab, licaminlimab and
  remtolumab, each from at least one and usually two independent sources.
* **Patent-disclosed binders: good breadth, document-level precision.** 221 TNF-titled US patent families contribute
  4620 sequences spanning VHH, VH/VL, full IgG chains, affibodies, FN3 monobodies, receptor fusions and ~1600 short
  peptide/CDR entries. Per-sequence target attribution inside a patent was not attempted.
* **Structurally confirmed binders: complete for P01375 as of this PDB snapshot.** All 52 P01375 entries were
  enumerated and every non-TNF protein chain in them was contact-checked; the remaining 31 entries are TNF alone,
  TNF with a small molecule, or TNF-derived peptides.
* **Non-antibody scaffold binders: partial.** Affibody and FN3-monobody classes are represented from patents; the
  tetranectin-CTLD antagonist (3L9J) and the poxvirus 2L decoy (3IT8) are present; the Affilin and the designed
  peptide literature are not.
* **Prior competition designs: weakest arm,** limited to the 35 public `human-tnfa` Proteinbase designs.
