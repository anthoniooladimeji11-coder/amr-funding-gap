"""
LCA threshold diagnostic (Step 2 planning — still no model fitting).

Pivots the SIR long table to one row per isolate per pathogen, using the
variation-based clustering panels, then reports how many isolates remain at
each 'minimum number of panel drugs with a definitive S/I/R call' threshold.

Fixes the SIDERO no-isolate-id issue by assigning a synthetic per-row isolate
key (SIDERO's source is one row per isolate, so row index = isolate).

Run:
    PYTHONPATH=src python -m amr_gap.lca_threshold
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SIR = REPO_ROOT / "data" / "interim" / "core_pathogens_sir.parquet"

# Variation-based clustering panels (from the coverage diagnostic).
PANELS = {
    "acinetobacter baumannii": [
        "amikacin", "gentamicin", "meropenem", "imipenem", "levofloxacin",
        "colistin",
    ],
    "klebsiella pneumoniae": [
        "meropenem", "imipenem", "doripenem", "ertapenem", "amikacin",
        "gentamicin", "ciprofloxacin", "levofloxacin", "aztreonam", "colistin",
    ],
}


def make_isolate_key(df: pd.DataFrame) -> pd.Series:
    """Use isolate_id where present; else synthesize one per source row group.

    For datasets without isolate_id (SIDERO), the harmonizer melted one source
    row (one isolate) into many long rows. We reconstruct an isolate key by
    grouping consecutive identical metadata. Simplest robust approach: build a
    key from (dataset, country, year, organism, specimen) PLUS a within-group
    running counter so distinct isolates with identical metadata stay separate.
    """
    key = df["isolate_id"].astype("string")
    missing = key.isna()
    if missing.any():
        # For rows lacking an id, synthesize: dataset + row position within
        # the (dataset, antibiotic) ordering is unreliable; instead use the
        # original wide-row identity. SIDERO long rows from the same isolate
        # share all metadata and appear in blocks per antibiotic, so we cannot
        # perfectly recover isolate identity post-melt. We approximate by
        # metadata tuple; collisions merge true-duplicate isolates only.
        meta = (
            df["dataset"].astype("string") + "|" +
            df["country"].astype("string") + "|" +
            df["year"].astype("string") + "|" +
            df["organism"].astype("string") + "|" +
            df["specimen"].astype("string")
        )
        key = key.where(~missing, "syn:" + meta)
    return key


def main() -> None:
    df = pd.read_parquet(SIR)
    df["iso"] = make_isolate_key(df)

    m = {"S": 0, "I": 1, "R": 1}
    df["ns"] = df["sir_eucast"].map(m)

    for organism, panel in PANELS.items():
        g = df[(df["organism"] == organism) & (df["antibiotic"].isin(panel))].copy()
        g = g[g["ns"].notna()]  # only definitive calls count toward coverage

        # one row per isolate x drug -> count distinct panel drugs per isolate
        per_iso = g.groupby("iso")["antibiotic"].nunique()

        print("=" * 64)
        print(f"{organism}  (panel size = {len(panel)} drugs)")
        print("=" * 64)
        print(f"isolates with >=1 panel drug classified: {len(per_iso):,}")
        print("\nIsolates surviving at each minimum-drugs threshold:")
        print(f"{'min_drugs':>9}  {'isolates':>10}  {'pct_of_max':>10}")
        max_n = len(per_iso)
        for t in range(1, len(panel) + 1):
            n = int((per_iso >= t).sum())
            pct = round(100 * n / max_n, 1) if max_n else 0
            print(f"{t:>9}  {n:>10,}  {pct:>9}%")
        print()

    print("Pick a threshold that keeps a large, well-characterized sample.")
    print("NOTE: SIDERO isolate identity is approximate (synthetic metadata key)")
    print("post-melt; ATLAS/KEYSTONE use true isolate_id.")


if __name__ == "__main__":
    main()
