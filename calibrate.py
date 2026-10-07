"""Build a labelled test set, run the prescreen over it, and locate the cut points.

Labels (constructed, not assumed):

``known``              a curated known binder, verbatim
``mutant_<pct>``       a known binder with pct% random substitutions
``graft``              a known binder's focus region on a different framework of the
                       same category  -> must be flagged
``framework_reuse``    a known binder's framework with a foreign focus region
                       -> must NOT be flagged
``unrelated_design``   a published design that is not in the target binder set
``shuffled``           a known binder with its residues shuffled (composition control)

Run with the curated reference set:

    python calibrate.py --target data/tnfa_binders.fasta --metadata data/tnfa_binders.csv \
        --unrelated path/to/unrelated_designs.fasta \
        --out calibration
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
from pathlib import Path

os.environ.setdefault("KMP_AFFINITY", "disabled")
sys.path.insert(0, str(Path(__file__).parent))

from prescreen import Config, mmseqs, screen, to_rows  # noqa: E402
from prescreen import fixtures  # noqa: E402
from prescreen.classify import classify  # noqa: E402
from prescreen.regions import extract, focus_region  # noqa: E402

MUTANT_FRACTIONS = (0.02, 0.05, 0.10, 0.20, 0.35)


def sample_references(
    references: dict, per_class: int = 8, seed: int = 0
) -> dict:
    """Stratified sample of references by header class, so the calibration set is a few
    hundred sequences rather than tens of thousands and every class is represented."""
    rng = random.Random(seed)
    by_class: dict = {}
    for ref_id, seq in references.items():
        parts = ref_id.split("|")
        cls = parts[1] if len(parts) > 1 else "unknown"
        by_class.setdefault(cls, []).append((ref_id, seq))
    out: dict = {}
    for cls, items in by_class.items():
        rng.shuffle(items)
        for ref_id, seq in items[:per_class]:
            out[ref_id] = seq
    return out


def build_labelled(
    references: dict,
    unrelated: dict,
    n_unrelated: int = 60,
    per_class: int = 8,
    seed: int = 0,
) -> tuple[dict, dict]:
    """Return ``(sequences, labels)`` for the calibration set.

    References are sampled stratified by class (``per_class`` each) to keep the labelled
    set to a few hundred sequences. Grafts and framework-reuse chimeras are built within
    the sample; whole-sequence identity for each arm is still measured against the FULL
    reference set at screen time.
    """
    rng = random.Random(seed)
    references = sample_references(references, per_class=per_class, seed=seed)
    seqs: dict = {}
    labels: dict = {}

    ann = {}
    for ref_id, seq in references.items():
        cl = classify(seq)
        regions = extract(seq, cl.category, cl.region_kind)
        chain, focus = focus_region(regions)
        ann[ref_id] = {
            "category": cl.category,
            "region_kind": cl.region_kind,
            "focus": focus,
        }

    def add(name, seq, label, note=""):
        if not seq or len(seq) < 5:
            return
        seqs[name] = seq
        labels[name] = {"label": label, "note": note}

    # 1. verbatim knowns and their mutant series
    for ref_id, seq in references.items():
        add(f"known__{ref_id}", seq, "known", ref_id)
        span = None
        for frac in MUTANT_FRACTIONS:
            n = max(1, int(round(frac * len(seq))))
            add(
                f"mut{int(frac * 100):02d}__{ref_id}",
                fixtures.mutant(seq, n, seed=hash(ref_id) % 10000, protect=span),
                f"mutant_{int(frac * 100):02d}",
                ref_id,
            )
        add(
            f"shuf__{ref_id}",
            "".join(rng.sample(seq, len(seq))),
            "shuffled",
            ref_id,
        )

    # 2. chimeras: same category, different framework
    by_cat: dict = {}
    for ref_id, a in ann.items():
        by_cat.setdefault((a["category"], a["region_kind"]), []).append(ref_id)
    for (category, region_kind), ids in by_cat.items():
        if region_kind == "whole" or len(ids) < 2:
            continue
        for donor, host in zip(ids, ids[1:] + ids[:1]):
            if donor == host:
                continue
            grafted = fixtures.graft(
                references[host], references[donor], category, region_kind
            )
            add(
                f"graft__{donor}_onto_{host}",
                grafted,
                "graft",
                f"donor={donor};host={host}",
            )
            foreign = ann[donor]["focus"]
            # Framework reuse: keep the host's framework, install a region from a
            # *different* category member scrambled so it is not any known paratope.
            scrambled = "".join(rng.sample(foreign, len(foreign))) if foreign else ""
            reuse = fixtures.framework_reuse(
                references[host], scrambled, category, region_kind
            )
            add(
                f"fwreuse__{host}",
                reuse,
                "framework_reuse",
                f"host={host}",
            )

    # 3. unrelated published designs
    pool = [k for k in unrelated if k not in references]
    rng.shuffle(pool)
    for name in pool[:n_unrelated]:
        add(f"unrel__{name}", unrelated[name], "unrelated_design", name)

    return seqs, labels


def sweep(rows: list, labels: dict) -> list:
    """Sweep single-metric thresholds and report counts per label class."""
    grid = [round(x / 100, 2) for x in range(50, 100, 2)] + [0.995]
    metrics = ("target_region_identity", "target_similarity", "prior_art_similarity")
    out = []
    for metric in metrics:
        for cut in grid:
            counts: dict = {}
            for row in rows:
                label = labels[row["id"]]["label"]
                value = float(row.get(metric) or 0.0)
                key = (label, "above" if value >= cut else "below")
                counts[key] = counts.get(key, 0) + 1
            entry = {"metric": metric, "cutoff": cut}
            for (label, side), n in counts.items():
                entry[f"{label}__{side}"] = n
            out.append(entry)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--unrelated", required=True)
    ap.add_argument("--out", default="calibration")
    ap.add_argument("--n-unrelated", type=int, default=60)
    ap.add_argument("--per-class", type=int, default=8)
    ap.add_argument("--skip-prior-art", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    references = mmseqs.read_fasta(args.target)
    unrelated = mmseqs.read_fasta(args.unrelated)
    seqs, labels = build_labelled(
        references, unrelated, n_unrelated=args.n_unrelated, per_class=args.per_class
    )
    print(f"labelled set: {len(seqs)} sequences from {len(references)} references")

    screened = screen(
        seqs,
        cfg=Config.from_env(),
        target_fasta=args.target,
        target_metadata=args.metadata,
        workdir=outdir / "work",
        skip_prior_art=args.skip_prior_art,
    )
    rows = to_rows(screened)
    for row in rows:
        row["label"] = labels[row["id"]]["label"]
        row["label_note"] = labels[row["id"]]["note"]

    table = outdir / "calibration_table.csv"
    with table.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    sweep_rows = sweep(rows, labels)
    sweep_path = outdir / "threshold_sweep.csv"
    keys = sorted({k for r in sweep_rows for k in r})
    with sweep_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(sweep_rows)

    (outdir / "timing.json").write_text(json.dumps(screened["timing"], indent=1))

    by_label: dict = {}
    for row in rows:
        by_label.setdefault(row["label"], {}).setdefault(row["verdict"], 0)
        by_label[row["label"]][row["verdict"]] += 1
    print(json.dumps(by_label, indent=1, sort_keys=True))
    print("timing:", screened["timing"])
    print("wrote", table, "and", sweep_path)


if __name__ == "__main__":
    main()
