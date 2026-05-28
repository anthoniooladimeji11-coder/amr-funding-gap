"""
Co-resistance diagnostic (Step 2 planning — previews cluster structure).

For each pathogen, at the chosen min-drugs threshold, this:
  1. pivots to one row per isolate (binary non-susceptible matrix),
  2. shows per-drug non-susceptibility rate (on the thresholded sample),
  3. shows the pairwise co-resistance correlation matrix (do drugs move in
     blocks? -> LCA will find clean classes),
  4. lists the most common observed resistance PATTERNS (a preview of the
     latent classes the LCA should recover).

Run:
    PYTHONPATH=src python -m amr_gap.coresistance
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SIR = REPO_ROOT / "data" / "interim" / "core_pathogens_sir.parquet"

PANELS = {
    "acinetobacter baumannii": (
        ["amikacin", "gentamicin", "meropenem", "imipenem", "levofloxacin",
         "colistin"],
        3,  # min drugs
    ),
    "klebsiella pneumoniae": (
        ["meropenem", "imipenem", "doripenem", "ertapenem", "amikacin",
         "gentamicin", "ciprofloxacin", "levofloxacin", "aztreonam", "colistin"],
        4,  # min drugs
    ),
}


def isolate_matrix(df: pd.DataFrame, organism: str, panel: list[str],
                   min_drugs: int) -> pd.DataFrame:
    g = df[(df["organism"] == organism) & (df["antibiotic"].isin(panel))].copy()
    m = {"S": 0, "I": 1, "R": 1}
    g["ns"] = g["sir_eucast"].map(m)
    g = g[g["ns"].notna()]
    # one row per isolate, columns = drugs, values = 0/1 (NaN if untested)
    wide = g.pivot_table(index="iso", columns="antibiotic", values="ns",
                         aggfunc="max")
    # keep isolates with at least min_drugs non-missing
    wide = wide[wide.notna().sum(axis=1) >= min_drugs]
    return wide


def main() -> None:
    df = pd.read_parquet(SIR)
    df["iso"] = df["isolate_id"]  # now guaranteed unique from harmonizer

    for organism, (panel, min_drugs) in PANELS.items():
        wide = isolate_matrix(df, organism, panel, min_drugs)
        print("=" * 68)
        print(f"{organism}  | isolates={len(wide):,} | min_drugs={min_drugs}")
        print("=" * 68)

        print("\nPer-drug non-susceptibility rate (thresholded sample):")
        rate = (wide.mean(numeric_only=True) * 100).round(1).sort_values(
            ascending=False)
        print(rate.to_string())

        print("\nPairwise co-resistance (Pearson corr of NS, blank=too sparse):")
        corr = wide.corr(min_periods=200).round(2)
        print(corr.to_string())

        print("\nTop 10 observed resistance patterns "
              "(1=non-susceptible, 0=susceptible, '-'=not tested):")
        pat = wide.reindex(columns=panel)
        as_str = pat.apply(
            lambda r: "".join("-" if pd.isna(v) else str(int(v)) for v in r),
            axis=1,
        )
        top = as_str.value_counts().head(10)
        # header showing drug order
        print("  pattern order:", " ".join(d[:4] for d in panel))
        for patt, n in top.items():
            print(f"  {patt}   {n:>7,}  ({round(100*n/len(wide),1)}%)")
        print()

    print("Block structure in the correlation matrix => LCA will find clean "
          "phenotypes. A few dominant patterns => low class count likely.")


if __name__ == "__main__":
    main()
