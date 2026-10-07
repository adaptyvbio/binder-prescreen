"""The target-binder arm: how close is a submission to a *known binder of this target*.

This is the part the existing novelty service has no equivalent for. Its question is not
"is this protein known" but "has someone already made this binder", and it is answered on
two levels:

**Whole sequence** — batched MMseqs2 search against the curated binder FASTA, plus an
exact-substring check. The substring check exists because the cheapest copy is a verbatim
one: a 25-aa peptide lifted from a paper, or a CDR loop pasted into a new framework, is a
substring of a reference even when no alignment is reported.

**Binding region** — the submission's focus region (CDR3 for antibody formats, projected
paratope for alternative scaffolds, whole sequence otherwise) compared against the
reference regions of *compatible* references, using the library's own region-identity
functions. This is what catches a known anti-TNF paratope grafted onto a fresh framework:
whole-sequence identity to the reference is low, region identity is near 1.

Region extraction over the reference set is cached to JSON beside the FASTA, so a batch
run pays for it once. The key is the FASTA's checksum *and* a fingerprint of the code that
produced the annotations: the classifier is an input to the cached result, so a cache keyed
on the FASTA alone keeps serving the old categories after the classifier changes, and the
reference side of the comparison silently disagrees with the query side.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from . import mmseqs
from .cdrref import CdrReference
from .classify import ANTIBODY_CATEGORIES, classify
from .compat import proteintyper
from .regions import FOCUS, extract, focus_region, identity_fn_for

MIN_SUBSTRING = 8


def known_binder_rank(row: dict, identity_cut: float, coverage_cut: float) -> tuple:
    """Sort key that picks the hit the known-binder rule is actually tested on.

    ``flags.evaluate`` calls a submission a known binder on a conjunction —
    ``fident >= tnf_known_identity`` AND ``qcov >= tnf_known_coverage`` — but the hit it
    reads was, until this key existed, chosen by ``argmax(fident * qcov)``. Those are
    different orderings: a hit at 0.92/0.85 (composite 0.78) loses to one at 0.85/0.95
    (composite 0.81), and the survivor fails the identity cut, so a real known binder goes
    unflagged on evidence that was never the strongest under the rule.

    Preferring a hit that clears both cuts, and breaking ties on the composite, makes the
    selection agree with the rule. With both cuts at 0.0 every row clears and the key
    degenerates to the plain composite ordering.

    The cost of this choice: when a rule-clearing hit has a lower composite than some
    other alignment, the reported ``target_similarity`` / ``target_hit`` is that
    rule-clearing hit rather than the single strongest alignment.
    """
    clears = (
        float(row.get("fident") or 0.0) >= identity_cut
        and float(row.get("qcov") or 0.0) >= coverage_cut
    )
    return (clears, float(row.get("similarity_check") or 0.0))


def _checksum(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()[:16]


def _annotator_fingerprint() -> str:
    """Fingerprint of the modules whose output ``ensure_annotations`` caches.

    Hashing the source rather than a hand-maintained version number means a change to the
    classifier cannot be forgotten: there is no bump to omit. A comment-only edit
    invalidates the cache too, which costs a few seconds of re-annotation — far cheaper
    than comparing a query annotated by the new classifier against references annotated
    by the old one.
    """
    h = hashlib.sha256()
    for mod in ("classify.py", "regions.py"):
        h.update((Path(__file__).parent / mod).read_bytes())
    return h.hexdigest()[:8]


@dataclass
class RegionHit:
    """Best region-level match of a query against the reference set."""

    # ``None`` means the region arm could not run — the query's focus region would not
    # extract, or no compatible reference existed to compare it with. That is not the same
    # as 0.0 ("compared against real references, and different"), and a caller that reads
    # the first as the second reports a check that never ran as a clean result.
    identity: float | None
    reference_id: str | None
    reference_region: str | None
    query_region: str
    chain: str | None
    label: str
    exact: bool = False
    binder_status: str | None = None
    named_agent: str | None = None
    citation: str | None = None
    supporting_exact: int = 0
    per_region: dict | None = None

    def as_dict(self) -> dict:
        return {
            "region_identity": None if self.identity is None else round(self.identity, 4),
            "region_compared": self.identity is not None,
            "region_exact": self.exact,
            "region_reference": self.reference_id,
            "region_reference_seq": self.reference_region,
            "region_query_seq": self.query_region,
            "region_chain": self.chain,
            "region_label": self.label,
            # Provenance of the matched reference, so a reviewer sees what the match
            # means: a confirmed binder, a patent claim, or a design already tested and
            # found not to bind.
            "binder_status": self.binder_status,
            "named_agent": self.named_agent,
            "citation": self.citation,
            "supporting_exact_regions": self.supporting_exact,
            "per_region": self.per_region or {},
        }


class TargetReference:
    """A curated set of known binders for one target, with cached region annotations."""

    #: Reference classes whose binding region is an IMGT CDR, and which are therefore
    #: represented by the CDR reference rather than by per-reference annotation.
    ANTIBODY_REF_CLASSES = frozenset(
        {
            "fv_domain_vh",
            "fv_domain_vl",
            "igg_chain_heavy",
            "igg_chain_light",
            "igg_fv",
            "fab",
            "vhh",
            "vhh_putative",
            "vnar",
            "peptide_or_cdr",
        }
    )

    def __init__(
        self,
        fasta: str | Path,
        metadata: str | Path | None = None,
        cache_path: str | Path | None = None,
        cdr_csv: str | Path | None = None,
    ):
        self.fasta = Path(fasta)
        if not self.fasta.exists():
            raise FileNotFoundError(f"target reference FASTA not found: {self.fasta}")
        self.sequences = mmseqs.read_fasta(self.fasta)
        self.checksum = _checksum(self.fasta)
        self.metadata = self._load_metadata(metadata)
        self.cdrs = CdrReference(cdr_csv) if cdr_csv and Path(cdr_csv).exists() else None
        self.cache_path = (
            Path(cache_path)
            if cache_path
            else self.fasta.with_suffix(
                f".regions.{self.checksum}.{_annotator_fingerprint()}.json"
            )
        )
        self.annotations: dict = {}

    @staticmethod
    def _load_metadata(metadata: str | Path | None) -> dict:
        if not metadata:
            return {}
        import csv

        out: dict = {}
        with Path(metadata).open(newline="") as fh:
            for row in csv.DictReader(fh):
                key = row.get("id") or row.get("name") or row.get("sequence_id")
                if key:
                    out[key] = row
        return out

    # -- reference annotation -------------------------------------------------------
    @staticmethod
    def parse_header(ref_id: str) -> dict:
        """Fields encoded in the FASTA header.

        ``TNFB00001|class|binder_status|source|accession[|named_agent]`` — the curation
        track's stable header format, so provenance travels with the sequence and does
        not depend on a separate metadata join.
        """
        parts = ref_id.split("|")
        keys = ("id", "class", "binder_status", "source", "accession", "named_agent")
        return {k: (parts[i] if i < len(parts) else "") for i, k in enumerate(keys)}

    @classmethod
    def reference_class(cls, ref_id: str) -> str | None:
        return cls.parse_header(ref_id).get("class") or None

    def _needs_annotation(self, ref_id: str) -> bool:
        """Whether this reference has to be numbered to be comparable.

        With a CDR reference loaded, the antibody-class references are already
        represented there in deduplicated, chain-typed form, so numbering all of them
        again would cost minutes and add nothing. Only the alternative-scaffold and
        miniprotein references need their paratopes projected.
        """
        if not self.cdrs:
            return True
        cls = self.reference_class(ref_id)
        if cls is None:
            return True
        return cls not in self.ANTIBODY_REF_CLASSES

    def ensure_annotations(self, rebuild: bool = False) -> dict:
        """Classify and region-annotate the reference sequences (cached on disk)."""
        if self.annotations and not rebuild:
            return self.annotations
        wanted = {k: v for k, v in self.sequences.items() if self._needs_annotation(k)}
        if self.cache_path.exists() and not rebuild:
            self.annotations = json.loads(self.cache_path.read_text())
            if set(self.annotations) == set(wanted):
                return self.annotations
        out: dict = {}
        for ref_id, seq in wanted.items():
            cl = classify(seq)
            regions = extract(seq, cl.category, cl.region_kind)
            out[ref_id] = {
                "category": cl.category,
                "region_kind": cl.region_kind,
                "length": cl.length,
                "regions": regions,
            }
        self.annotations = out
        try:
            self.cache_path.write_text(json.dumps(out))
        except OSError:
            pass
        return out

    def compatible(self, category: str, region_kind: str) -> list:
        """Reference ids whose regions are comparable with a query of this category."""
        self.ensure_annotations()
        if region_kind == "antibody_cdrs":
            keep = set(ANTIBODY_CATEGORIES)
        elif region_kind == "scaffold_paratope":
            keep = {category}
        else:
            keep = None
        return [
            ref_id
            for ref_id, ann in self.annotations.items()
            if keep is None or ann["category"] in keep
        ]

    # -- whole-sequence arm ---------------------------------------------------------
    def search_whole(
        self,
        query_fasta: str | Path,
        workdir: str | Path,
        settings: dict | None = None,
        mmseqs_bin: str = "mmseqs",
        threads: int = 0,
        identity_cut: float = 0.0,
        coverage_cut: float = 0.0,
    ) -> dict:
        """Best whole-sequence hit per query against the reference FASTA.

        Queries are routed by length into the same two parameter sets the public
        databases use (:func:`prescreen.mmseqs.settings_for_length`), so a 12-aa peptide
        entry is not silently missed here either. Pass ``settings`` to force one set.

        Returns ``{query_id: row}`` where ``row`` carries MMseqs2 fields plus
        ``similarity_check = fident * qcov`` (the composite the deployed service uses).
        """
        from .refdb import LENGTH_SPLIT

        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        queries = mmseqs.read_fasta(query_fasta)
        groups = {"short": {}, "standard": {}}
        for name, seq in queries.items():
            groups["short" if len(seq) < LENGTH_SPLIT else "standard"][name] = seq

        best: dict = {}
        for label, group in groups.items():
            if not group:
                continue
            gfasta = mmseqs.write_fasta(group, workdir / f"target_q_{label}.fasta")
            out = workdir / f"target_hits.{label}.m8"
            mmseqs.easy_search(
                gfasta,
                self.fasta,
                out,
                settings=settings or mmseqs.settings_for_length(
                    1 if label == "short" else 10**6
                ),
                mmseqs_bin=mmseqs_bin,
                threads=threads,
                tmp_dir=workdir / "tmp",
            )
            for row in mmseqs.parse_m8(out):
                row["similarity_check"] = row["fident"] * row["qcov"]
                row["profile"] = label
                hdr = self.parse_header(row["target"])
                row["binder_status"] = hdr.get("binder_status") or ""
                row["named_agent"] = hdr.get("named_agent") or ""
                row["reference_class"] = hdr.get("class") or ""
                # Whether the CDR index holds anything for this reference's class. When it
                # does not, the region arm cannot confirm a match against this reference,
                # and the flag logic must not treat its silence as evidence of novelty.
                row["region_reference_available"] = bool(
                    self.cdrs is None
                    or row["reference_class"] in self.cdrs.covered_classes
                )
                prev = best.get(row["query"])
                if prev is None or known_binder_rank(
                    row, identity_cut, coverage_cut
                ) > known_binder_rank(prev, identity_cut, coverage_cut):
                    best[row["query"]] = row
        return best

    def substring_match(self, seq: str) -> dict | None:
        """Exact containment of the query (or a reference) in the other, either direction."""
        if len(seq) < MIN_SUBSTRING:
            return None
        for ref_id, ref_seq in self.sequences.items():
            if seq in ref_seq:
                return {
                    "reference": ref_id,
                    "direction": "query_in_reference",
                    "matched_length": len(seq),
                }
            if len(ref_seq) >= MIN_SUBSTRING and ref_seq in seq:
                return {
                    "reference": ref_id,
                    "direction": "reference_in_query",
                    "matched_length": len(ref_seq),
                }
        return None

    # -- region arm -----------------------------------------------------------------
    def compare_cdrs(self, query_regions: dict) -> dict:
        """Match an antibody-format query's CDRs against the known-paratope reference.

        Returns ``{"focus": CdrMatch|None, "focus_chain": str|None,
        "supporting_exact": int, "per_region": {...}}``, where ``supporting_exact``
        counts how many of CDR1/2/3 matched exactly. CDR1 and CDR2 are germline-shared
        across unrelated antibodies, so they corroborate a CDR3 match but must not carry a
        flag on their own.
        """
        matches = self.cdrs.match_chains(query_regions)
        focus, best_chain = None, None
        supporting, per_region = 0, {}
        for chain, labels in matches.items():
            for label, match in labels.items():
                per_region[f"{chain}:{label}"] = match.as_dict()
                if match.exact:
                    supporting += 1
                if label == FOCUS and (focus is None or match.identity > focus.identity):
                    focus, best_chain = match, chain
        return {
            "focus": focus,
            "focus_chain": best_chain,
            "supporting_exact": supporting,
            "per_region": per_region,
        }

    def compare_region(
        self, query_regions: dict, category: str, region_kind: str
    ) -> RegionHit:
        """Best focus-region identity of a query against compatible references.

        Antibody formats are compared against the deduplicated CDR reference when one is
        loaded; alternative scaffolds against the projected paratopes of same-class
        references; plain sequences against reference focus regions and whole sequences.
        """
        if region_kind == "antibody_cdrs" and self.cdrs:
            result = self.compare_cdrs(query_regions)
            focus = result["focus"]
            if focus is None:
                chain, q_focus = focus_region(query_regions)
                # No CDR of this query matched anything in the reference index. If the
                # query had no extractable focus region at all, the arm did not run.
                return RegionHit(
                    None if not q_focus else 0.0, None, None, q_focus, chain, FOCUS
                )
            return RegionHit(
                identity=focus.identity,
                reference_id=focus.reference_id,
                reference_region=focus.reference_seq,
                query_region=focus.query,
                chain=result["focus_chain"],
                label=FOCUS,
                exact=focus.exact,
                binder_status=focus.status,
                named_agent=focus.agents,
                citation=focus.citation,
                supporting_exact=result["supporting_exact"],
                per_region=result["per_region"],
            )

        cdr_novelty = proteintyper("cdr_novelty")
        identity = identity_fn_for(category)
        chain, q_focus = focus_region(query_regions)
        if not q_focus:
            return RegionHit(None, None, None, "", None, FOCUS)

        candidates = self.compatible(category, region_kind)
        if not candidates:
            # No reference of a comparable category exists — e.g. a monobody query when
            # every curated monobody reference failed its own framework gate. Comparing
            # against nothing scores 0.0, which reads as "checked and clean".
            return RegionHit(None, None, None, q_focus, chain, FOCUS)
        hit_regions, hit_ids = [], []
        for ref_id in candidates:
            ann = self.annotations[ref_id]
            if region_kind == "whole":
                # A plain sequence is compared against both the reference's own focus
                # region (catching a loop lifted out of a known binder) and the whole
                # reference sequence.
                ref_chain, ref_focus = focus_region(ann["regions"])
                refs = {"W": {FOCUS: ref_focus or self.sequences[ref_id]}}
            else:
                refs = ann["regions"]
            hit_regions.append(refs)
            hit_ids.append(ref_id)

        if region_kind == "whole":
            query_for_compare = {"W": {FOCUS: q_focus}}
        else:
            query_for_compare = query_regions

        chains, _ = cdr_novelty.compare_regions(
            query_for_compare,
            hit_regions,
            hit_meta=hit_ids,
            focus=FOCUS,
            identity_fn=identity,
        )
        best = RegionHit(0.0, None, None, q_focus, chain, FOCUS)
        for chain_key, labels in chains.items():
            entry = labels.get(FOCUS)
            if not entry:
                continue
            if entry["identity"] > best.identity:
                best = RegionHit(
                    identity=entry["identity"],
                    reference_id=entry.get("closest_hit"),
                    reference_region=entry.get("closest_hit_cdr"),
                    query_region=entry["query"],
                    chain=chain_key,
                    label=FOCUS,
                )
        return best
