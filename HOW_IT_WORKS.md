# How the TNF-α prescreen works

## The question it answers

When a submission arrives, the filter decides — from sequence alone — whether
it is **already known**. "Known" means one of two things:

1. **Already public** — the sequence (or a near-identical one) is in a protein database, a
   patent, or an existing published design.
2. **Already a TNF-α binder** — it matches a known anti-TNF binder, either over the whole
   sequence or just in its **binding region** (the paratope).

It does **not** ask whether the scaffold is new. Reusing a published framework is allowed
in the competition, so a brand-new paratope on an old scaffold is a legitimate entry and
passes. This is the one thing a naive "BLAST it and see if it's similar to a known binder"
check gets wrong, and the reason the filter is built the way it is.

## Why whole-sequence similarity is not enough

Antibodies — the bulk of known anti-TNF binders — nearly all share the same germline
framework. Most of an antibody's sequence is that shared scaffolding; only the CDR loops
(and above all CDR3) actually touch the antigen. Measured on the reference set, a single
probe nanobody matches **1,843 of ~5,000 entries** at high coverage on framework alone —
only 15 of those are genuine near-identical hits.

So we need to compare the **binding region separately**:

- If a submission reuses a **known framework but carries a new CDR3**, its whole sequence
  looks ~90% identical to a known binder, but its paratope does not match → **passes**
  (new binder, allowed).
- If a submission puts a **known anti-TNF CDR3 onto a fresh framework**, its whole sequence
  looks unrelated, but its paratope matches → **flagged** (known binder in disguise).

That discrimination is the core of the filter.

## The pipeline

1. **Classify** the submission by format (nanobody, Fab, IgG, affibody, monobody, DARPin,
   miniprotein, peptide …) from sequence. The format decides which binding region to read.
2. **Extract the binding region** — antibody CDRs by IMGT numbering, a projected paratope
   for alternative scaffolds, or the whole sequence for a short peptide.
3. **Prior-art search** — one batched MMseqs2 search per public database (PDB, SwissProt,
   PLAbDab + PLAbDab-nano for single-domain formats, Thera-SAbDab, THPdb, published
   Proteinbase designs; USPTO patents). Short sequences get a second, sensitive
   pass so a 12-mer peptide copied from a paper is not silently missed. Each top hit's own
   sequence is re-numbered and its paratope compared with the submission's — so a shared
   framework is reported as context, not as a known molecule.
4. **Known-binder search** — the submission against the curated TNF-α binder set, both
   whole-sequence and region-level, plus an exact-substring check for verbatim copies.
5. **Within-batch duplicates** — clusters repeated submissions in the same batch.
6. **Verdict** — the single most disqualifying flag that fired. Everything that fired is
   listed; framework-reuse observations are reported as *context* and never disqualify.

## What you get back

One row per submission with the verdict, every flag, the matched reference and its
provenance (which database, which named agent, and whether that reference is a *confirmed
binder*, a *patent claim*, or a design that was *tested and did not bind*), plus the
identity/coverage numbers behind each call. A hit to a tested-non-binder design tells a
participant "this exact design was already tried and failed" — a different message from a
hit to a confirmed binder.

## The verdicts

Ordered from most to least disqualifying (the verdict is the first that fires):

| verdict | meaning |
|---------|---------|
| `known_target_binder` | matches a known anti-TNF binder (whole-sequence **and** paratope for antibody formats; whole-sequence for peptides/miniproteins) |
| `target_region_match` | carries a known anti-TNF paratope (CDR3 / projected paratope) on any framework |
| `verbatim_target_fragment` | is an exact fragment of a known binder (short peptides) |
| `existing_design` | matches an existing published design |
| `known_sequence` | matches a public sequence (PDB / SwissProt / PLAbDab / patents) |
| `target_binder_homolog` | a close relative of a known binder (non-antibody) |
| `near_known_sequence` | a near-match to a public sequence (non-antibody) |
| `batch_duplicate` | duplicates another submission in the same batch |
| `no_reference_signal` | too short to search reliably and found nothing |
| `pass` | none of the above — not known, submit for assay |

## What is in the reference database

4,964 unique known or claimed TNF-α binders, each with traceable provenance: 28 confirmed
binders, 35 designs tested and found not to bind, and 4,901 patent/document-level claims.
Plus 4,075 deduplicated anti-TNF binding regions as the paratope reference. By format:

| category | n | what it is |
|----------|---|-----------|
| antibody Fv / CDR peptides | 1,626 | CDR and peptide fragments |
| antibody VH domains | 1,335 | heavy variable domains |
| antibody VL domains | 1,040 | light variable domains |
| nanobody / VHH | 298 | camelid single-domain |
| IgG heavy / light chains | 460 | full antibody chains |
| affibody | 84 | protein-A Z-domain scaffolds |
| receptor-derived | 33 | TNFR ectodomains (etanercept-like) |
| monobody (FN3) | 10 | fibronectin domain scaffolds |
| miniprotein / designed | 20 | de novo designs |
| VNAR | 2 | shark single-domain |

By primary source (a few entries cite more than one): USPTO patents (4,423), PLAbDab
antibodies (479), Proteinbase designs (35), PDB structures (25, all contact-verified
against TNF-α), UniProt (1), THPdb (1).

## Example known binders, by category

Verbatim submission of any of these — or a close variant — is what the filter is built to catch. Identifiers are the stable `id` column of `tnfa_binders.csv`.

**Approved / clinical antibodies** (searched as heavy + light chains, and as structural Fv)

| agent | id | format | evidence |
|-------|-----|--------|----------|
| Adalimumab (Humira) | TNFB01956 | Fv (VL) | PLAbDab + patent |
| Infliximab (Remicade) | TNFB01787 | Fv (VL) | PLAbDab + patent |
| Certolizumab (Cimzia) | TNFB01823 | Fv (VL) | PLAbDab + patent |
| Golimumab (Simponi) | TNFB00009/00010 | Fv, confirmed | PDB 5YOY |

**Nanobodies / VHH** (single-domain; also searched against the dedicated PLAbDab-nano arm)

| id | note | sequence (start) |
|----|------|------------------|
| TNFB00125 | Ozoralizumab, confirmed binder (PDB) | `QVQLVESGGGLVQ…` |
| TNFB00013 | anti-TNF VHH, confirmed binder (PDB) | `QVQLVESGGGLVQPGGSLRLSCAASGFTFSNYWMY…` |
| TNFB00309 | VNAR, confirmed binder (PDB) | `ARVDQTPQTITKETGESLTINCVLRDSNC…` |

**Receptor-derived** (TNF binds its receptor; etanercept is a TNFR2–Fc fusion)

| id | note |
|----|------|
| TNFB03288 | Etanercept (Enbrel), confirmed — TNFR2 ectodomain + IgG1 Fc, 467 aa |
| TNFB03263 | TNFR ectodomain, confirmed binder (PDB + patent) |

**Alternative scaffolds** (patent-disclosed anti-TNF binders)

| id | scaffold | sequence (start) |
|----|----------|------------------|
| TNFB03146 | affibody (protein-A Z-domain) | `VDNKFNKEAAWAPFEIQHLPNLNHPQNDAFIDSLTDDPSQSANLLAEAKK` |
| TNFB03230 | monobody (FN3 / fibronectin) | `SPPKDLVVTEVTEETVNLAWDNEMRVTEYLVVYTPT…` |

**Peptides**

| id | note | sequence |
|----|------|----------|
| TNFB03773 | confirmed binder, 12-mer (PDB) | `ACPPCLWQVLCG` |

**Designs tested and found not to bind** 
A hit here means "already tried, did not work"

| id | note |
|----|------|
| TNFB03240 | Proteinbase miniprotein design, tested non-binder, 81 aa |
| TNFB03241 | Proteinbase miniprotein design, tested non-binder, 81 aa |

## Examples (from calibration)

| submission type | verdict | why |
|-----------------|---------|-----|
| exact adalimumab / known VHH / affibody | `known_target_binder` | whole-sequence + paratope both match |
| known anti-TNF CDR3 on a foreign framework | `target_region_match` | paratope matches, whole-sequence does not |
| known framework carrying a new CDR3 | `pass` | framework reuse is allowed; paratope is new |
| unrelated design for another target | `existing_design` | already published on Proteinbase |

On the 1,026-sequence calibration set: known binders flagged 100 %, grafted known
paratopes 97 %, reused frameworks correctly pass, shuffled controls pass 98 %. The
unrelated designs are drawn from the published-design corpus the screen searches, so they
all match themselves — see the note in `PRESCREEN_SPEC.md` §5.
