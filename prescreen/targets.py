"""The roster of curated known-binder sets, one per target.

The prior-art arm is target-agnostic - it searches the same public corpora whatever is
being screened - so the only thing that changes between targets is the curated set the
*known-binder* arm is built from: a FASTA, an optional provenance CSV and an optional
CDR index. This module is the one place those three paths are chosen.

A set is picked by slug (``--target tnfa``) and read out of ``data/targets/targets.json``.
There is deliberately no way to pass a loose FASTA instead: a set is three files that
have to agree with each other, and the earlier flags for overriding them individually
made it possible to screen one target's sequences against another target's paratopes.
Adding a target means registering it, which is a four-step recipe in the README.

Registering a target does not calibrate it. The cut points in :mod:`prescreen.config`
were measured on the TNF-alpha set, and a new target inherits them untested until
``calibrate.py`` is run over its own references; ``notes`` is the place to say so.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

#: Where the registry and the per-target directories live.
TARGETS_DIR = Path(__file__).parent / "data" / "targets"
REGISTRY_FILE = TARGETS_DIR / "targets.json"

#: Below this a reference is unreachable by both search paths: MMseqs2's short profile
#: indexes 6-mers (``refdb.PROFILE_SHORT``) and the exact-substring fallback has its own
#: floor (``target.MIN_SUBSTRING``). References this short are carried in the set but
#: cannot be matched, which is worth a warning rather than a silent gap.
MIN_SEARCHABLE_LENGTH = 6


class UnknownTarget(ValueError):
    """Raised for a slug that is not in the registry, naming the ones that are."""


@dataclass(frozen=True)
class TargetSet:
    """One target's curated known-binder set."""

    slug: str
    name: str
    fasta: Path
    uniprot: str | None = None
    metadata: Path | None = None
    cdr_csv: Path | None = None
    notes: str = ""

    def describe(self) -> str:
        label = f"{self.slug} ({self.name})"
        return f"{label} - {self.uniprot}" if self.uniprot else label


def _load() -> dict:
    if not REGISTRY_FILE.exists():
        return {"default": "tnfa", "targets": {}}
    return json.loads(REGISTRY_FILE.read_text())


def _build(slug: str, entry: dict) -> TargetSet:
    """Turn one registry entry into a :class:`TargetSet`, resolving its paths.

    Paths are relative to the entry's own directory under ``data/targets/``, so the
    registry stays readable and a target directory can be copied between checkouts
    without rewriting it.
    """
    base = TARGETS_DIR / (entry.get("dir") or slug)

    def _path(key: str) -> Path | None:
        value = entry.get(key)
        return (base / value) if value else None

    fasta = _path("fasta")
    if fasta is None:
        raise ValueError(f"target {slug!r} in {REGISTRY_FILE} has no 'fasta' entry")
    return TargetSet(
        slug=slug,
        name=entry.get("name") or slug,
        fasta=fasta,
        uniprot=entry.get("uniprot"),
        metadata=_path("metadata"),
        cdr_csv=_path("cdr_csv"),
        notes=entry.get("notes") or "",
    )


def registry() -> dict[str, TargetSet]:
    """Every registered target, keyed by slug."""
    data = _load()
    return {slug: _build(slug, entry) for slug, entry in (data.get("targets") or {}).items()}


def available() -> list[str]:
    """The registered slugs, for CLI choices and request validation."""
    return sorted(registry())


def default_slug() -> str:
    data = _load()
    return data.get("default") or "tnfa"


def resolve(slug: str | None = None) -> TargetSet:
    """The registered set for ``slug``, or the registry's default."""
    sets = registry()
    slug = slug or default_slug()
    if slug not in sets:
        raise UnknownTarget(
            f"unknown target {slug!r}; registered: {', '.join(available()) or '(none)'}. "
            "Add one with the recipe in the README, under 'Adding a target'."
        )
    return sets[slug]


def check(target: TargetSet) -> list[str]:
    """Warnings worth printing before a run: missing files, unsearchable references."""
    problems: list[str] = []
    if not target.fasta.exists():
        problems.append(f"known-binder FASTA not found: {target.fasta}")
        return problems
    if target.metadata is not None and not target.metadata.exists():
        problems.append(f"metadata CSV not found: {target.metadata}")
    if target.cdr_csv is not None and not target.cdr_csv.exists():
        problems.append(
            f"CDR index not found: {target.cdr_csv} - the binding-region arm will fall "
            "back to whole-sequence annotation for antibody references"
        )
    short = _count_short(target.fasta)
    if short:
        problems.append(
            f"{short} reference(s) shorter than {MIN_SEARCHABLE_LENGTH} aa cannot be "
            "matched by either the MMseqs2 short profile or the substring fallback"
        )
    return problems


def _count_short(fasta: Path) -> int:
    n, length = 0, 0
    with fasta.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if 0 < length < MIN_SEARCHABLE_LENGTH:
                    n += 1
                length = 0
            else:
                length += len(line.strip())
    if 0 < length < MIN_SEARCHABLE_LENGTH:
        n += 1
    return n
