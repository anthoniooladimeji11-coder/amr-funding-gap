"""
LCA coverage diagnostic (Step 2 planning — does NOT fit any model).

Purpose: choose the clustering drug panel and the per-isolate minimum-drugs
threshold from real numbers, before building the LCA. Reports, per pathogen:
  1. how isolate identity is available per dataset (can we pivot to one row
     per isolate?),
  2. per-drug: coverage (how many isolates tested) and non-susceptibility rate
     among classified calls (drugs with ~0% or ~100% NS, or very low coverage,
     are poor for clustering),
  3. the distribution of "number of panel drugs with a definitive S/I/R call"
     per isolate (drives the minimum-drugs threshold).

Run:
    PYTHONPATH=src python -m amr_gap.lca_diagnostic
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SIR = REPO_ROOT / "data" / "interim" / "core_pathogens_sir.parquet"

PANEL = [
    "cefiderocol", "meropenem", "imipenem", "doripenem", "ertapenem",
    "colistin", "amikacin", "gentamicin", "tobramycin", "ciprofloxacin",
    "levofloxacin", "ceftazidime-avibactam", "ceftolozane-tazobactam",
    "aztreonam", "trimethoprim-sulfamethoxazole",
]


def main() -> None:
    df = pd.read_parquet(SIR)
    df = df[df["antibiotic"].isin(PANEL)].copy()

    # Binary non-susceptible flag from EUCAST calls; U -> missing (NA).
    # S -> 0 (susceptible), I/R -> 1 (non-susceptible), U -> NA.
    m = {"S": 0, "I": 1, "R": 1}
    df["ns"] = df["sir_eucast"].map(m)  # U and anything else -> NaN

    print("=" * 70)
    print("ISOLATE IDENTITY CHECK (can we pivot to one row per isolate?)")
    print("=" * 70)
    for ds, g in df.groupby("dataset"):
        n_rows = len(g)
        n_id = g["isolate_id"].notna().sum()
        n_unique = g["isolate_id"].nunique(dropna=True)
        print(f"{ds:10s}  rows={n_rows:>9,}  with isolate_id={n_id:>9,}  "
              f"unique_ids={n_unique:>9,}")
    print("\nNote: datasets WITHOUT isolate_id need a synthetic key "
          "(row-group) before pivoting. We'll address in the LCA build.\n")

    for organism, gorg in df.groupby("organism"):
        print("=" * 70)
        print(f"PATHOGEN: {organism}   (panel rows: {len(gorg):,})")
        print("=" * 70)

        # Per-drug coverage + NS rate among classified (non-missing) calls.
        rows = []
        for drug, gd in gorg.groupby("antibiotic"):
            classified = gd["ns"].notna().sum()
            total = len(gd)
            ns_rate = gd["ns"].mean()  # ignores NaN
            rows.append((drug, total, classified,
                         round(100 * classified / total, 1) if total else 0,
                         round(100 * ns_rate, 1) if classified else float("nan")))
        cov = pd.DataFrame(
            rows,
            columns=["drug", "rows", "classified", "classified_pct", "ns_rate_pct"],
        ).sort_values("classified", ascending=False)
        print("\nPer-drug coverage and non-susceptibility (among classified):")
        print(cov.to_string(index=False))
        print("\n  -> Good clustering drugs: high 'classified' AND ns_rate "
              "not near 0 or 100 (i.e. real variation).")

    print("\nDONE. Use this to pick (a) the clustering panel and (b) the "
          "min-drugs-per-isolate threshold in the LCA step.")


if __name__ == "__main__":
    main()
