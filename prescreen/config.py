"""Cut points, database locations and tunables, all overridable without a code edit.

Every threshold here is a *prescreen* threshold: its job is to triage submissions into
"needs a human look" versus "carry on", not to assign a novelty level. They are
deliberately separate from the novelty scales, whose cut points answer a different
question (is the scaffold de novo) and are inherited from the published four-level
scale.

Defaults are calibrated values; see ``PRESCREEN_SPEC.md`` for the data behind each one.
Override at runtime with ``Config.from_env()`` (``PRESCREEN_*`` variables) or by passing a
``Config`` into :func:`prescreen.screen`.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field, asdict

DB_ROOT_ENV = "PRESCREEN_DB_ROOT"


def resolve_db_root(explicit: str | None = None) -> str | None:
    """A ``:``-separated list of existing database roots, or ``None``.

    Arm lookup itself lives in :mod:`prescreen.refdb` (``roots`` argument →
    ``$PRESCREEN_DB_ROOT`` → its own fallbacks), so this only has to produce the root
    list to hand over. Several roots are normal: the verified databases sit on one mount
    and the arms built for this filter in another directory.
    """
    from . import refdb

    candidates = []
    for value in (explicit, os.environ.get(DB_ROOT_ENV)):
        if value:
            candidates.extend(value.split(":"))
    candidates.extend(refdb.FALLBACK_ROOTS)
    roots = [c for c in candidates if c and os.path.isdir(c)]
    seen, ordered = set(), []
    for root in roots:
        if root not in seen:
            seen.add(root)
            ordered.append(root)
    return ":".join(ordered) or None


@dataclass
class Config:
    """Prescreen configuration."""

    # --- prior-art arm (whole-sequence identity x query coverage) ---
    known_sequence: float = 0.95
    near_known_sequence: float = 0.80

    # --- target-binder arm (whole sequence) ---
    # Cut points measured by the curation track on the known-binder corpus: the antibody
    # set cross-hits itself at the framework level (a probe VHH returns 647 hits >= 0.70
    # fident, pure Ig background), so a known-binder call needs identity AND coverage
    # both high. tnf_known_identity requires fident; tnf_known_coverage requires qcov.
    tnf_known_identity: float = 0.90
    tnf_known_coverage: float = 0.80
    # Above this whole-sequence composite a non-antibody binder (miniprotein, peptide,
    # affibody) is treated as the known molecule; the antibody formats are judged on the
    # region instead.
    tnf_known: float = 0.90
    tnf_homolog: float = 0.70

    # --- target-binder arm (binding region: CDR3 / scaffold paratope) ---
    # A known anti-target paratope on any framework: exact CDR3 match, or CDR3 Levenshtein
    # identity >= region_match with matching chain type. CDR1/CDR2 are germline-shared and
    # only corroborate.
    region_match: float = 0.80
    region_min_length: int = 5

    # --- paratope-aware prior art ---
    # A whole-sequence match to a public reference only means the *molecule* is known if
    # the binding region matches too; below this the match is a shared framework, which
    # is explicitly allowed in these competitions.
    prior_art_region_known: float = 0.90
    prior_art_top_hits: int = 3

    # --- within-batch duplicate clustering ---
    batch_cluster_identity: float = 0.95
    batch_cluster_coverage: float = 0.80

    # --- search behaviour ---
    db_root: str | None = None  # ':'-separated list of database roots
    short_sequence_length: int = 50  # mirrors refdb.LENGTH_SPLIT
    mmseqs_bin: str = "mmseqs"
    threads: int = 0  # 0 = let mmseqs decide
    # The patent arm (10.2M sequences) costs seconds per query and belongs in an
    # asynchronous batch, not in the submission path.
    include_patent_arm: bool = False
    # Organiser mode adds the internal Proteinbase corpus — other entrants' unpublished
    # submissions. Never enable it for output a competitor will see.
    organiser_mode: bool = False
    target_name: str = "TNF-alpha"
    target_fasta: str | None = None
    target_metadata: str | None = None
    target_cdr_csv: str | None = None
    extras: dict = field(default_factory=dict)

    @classmethod
    def from_env(cls) -> "Config":
        """Build a config, letting ``PRESCREEN_<FIELD>`` environment variables win."""
        kwargs: dict = {}
        for name, f in cls.__dataclass_fields__.items():
            raw = os.environ.get(f"PRESCREEN_{name.upper()}")
            if raw is None:
                continue
            if f.type in ("float", float):
                kwargs[name] = float(raw)
            elif f.type in ("int", int):
                kwargs[name] = int(raw)
            else:
                kwargs[name] = raw
        cfg = cls(**kwargs)
        cfg.db_root = resolve_db_root(cfg.db_root)
        return cfg

    def as_dict(self) -> dict:
        return asdict(self)
