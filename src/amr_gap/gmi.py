"""
Step 4b - Global Misalignment Index (GMI).

Compares each pathogen's BURDEN share (from GRAM 2019 attributable deaths/DALYs)
against its FUNDING share (from the FAP). The misalignment per pathogen is
funding_share - burden_share:
  > 0  -> RELATIVELY OVERFUNDED for its burden
  < 0  -> RELATIVELY UNDERFUNDED
  ~ 0  -> aligned

The GMI itself is a single summary: sum of absolute misalignments across
pathogens, or half of it (so it lies in [0, 1] as a misallocation index).

Variants (cartesian product, all reported):
  - burden metric: deaths (HEADLINE) | dalys (SENSITIVITY)
  - funding case:  base (HEADLINE) | sens_proportional | sens_equal

Both inputs (FAP + burden weights) are at the same resolution: GLOBAL per
pathogen (D-005, D-016). Funding is aggregated across all FAP years; burden is
2019. We document this temporal mismatch as a stated assumption — the misalign-
ment compares cumulative R&D effort against the burden it is meant to address.

Run:
    PYTHONPATH=src python -m amr_gap.gmi
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
PROC = REPO_ROOT / "data" / "processed"
EXT = REPO_ROOT / "data" / "external"
FAP_PATH = PROC / "fap_pathogen_year.parquet"
BURDEN_PATH = EXT / "burden_weights.csv"

# Map FAP genus -> the species name the burden file uses.
GENUS_TO_SPECIES = {
    "escherichia":    "escherichia coli",
    "klebsiella":     "klebsiella pneumoniae",
    "acinetobacter":  "acinetobacter baumannii",
    "pseudomonas":    "pseudomonas aeruginosa",
    "staphylococcus": "staphylococcus aureus",
    "streptococcus":  "streptococcus pneumoniae",
}


def shares(values: pd.Series) -> pd.Series:
    """Normalize to shares summing to 1 (or 100%)."""
    t = values.sum()
    return values / t if t > 0 else values * 0.0


def main() -> None:
    fap = pd.read_parquet(FAP_PATH)
    burden = pd.read_csv(BURDEN_PATH, comment="#")
    assert (burden["status"] == "OK").all(), "burden weights have un-verified rows"

    # ---- Funding side: aggregate FAP across all years per pathogen ----
    fund = fap.groupby("pathogen")[
        ["base_usd", "sens_proportional_usd", "sens_equal_usd"]].sum()
    fund["species"] = fund.index.map(GENUS_TO_SPECIES)
    fund = fund.reset_index().rename(columns={"pathogen": "genus"})

    # ---- Burden side: pivot to one row per species, columns = metric ----
    bw = burden[burden["gbd_region"] == "Global"].copy()
    burden_wide = bw.pivot_table(index="pathogen", columns="metric",
                                 values="value", aggfunc="first")
    burden_wide.index.name = "species"
    burden_wide = burden_wide.reset_index()

    # ---- Join on species ----
    g = fund.merge(burden_wide, on="species", how="inner").set_index("species")
    if len(g) != len(fund):
        missing = set(fund["species"]) - set(g.index)
        raise SystemExit(f"missing burden rows for: {missing}")

    print("=== Pathogen inputs (USD totals + burden) ===")
    show = g[["genus", "base_usd", "deaths", "dalys"]].copy()
    show["base_usd_M"] = (show["base_usd"] / 1e6).round(0)
    show["deaths_k"] = (show["deaths"] / 1e3).round(0)
    show["dalys_k"] = (show["dalys"] / 1e3).round(0)
    print(show[["genus", "base_usd_M", "deaths_k", "dalys_k"]].to_string())

    # ---- Compute shares and misalignment for every (funding, burden) combo ----
    out_rows = []
    fund_cases = {
        "base":         "base_usd",
        "sens_prop":    "sens_proportional_usd",
        "sens_equal":   "sens_equal_usd",
    }
    burden_metrics = ["deaths", "dalys"]
    summary = {}
    for fcase, fcol in fund_cases.items():
        for bm in burden_metrics:
            f_share = shares(g[fcol])
            b_share = shares(g[bm])
            mis = f_share - b_share  # positive = relatively overfunded
            gmi = float(mis.abs().sum() / 2.0)  # half-sum-abs in [0, 1]
            summary[(fcase, bm)] = gmi
            for sp in g.index:
                out_rows.append({
                    "pathogen": sp,
                    "genus": g.loc[sp, "genus"],
                    "funding_case": fcase,
                    "burden_metric": bm,
                    "funding_usd": float(g.loc[sp, fcol]),
                    "burden_value": float(g.loc[sp, bm]),
                    "funding_share": float(f_share[sp]),
                    "burden_share": float(b_share[sp]),
                    "misalignment": float(mis[sp]),  # f_share - b_share
                })

    long = pd.DataFrame(out_rows)
    PROC.mkdir(parents=True, exist_ok=True)
    out_path = PROC / "gmi.parquet"
    long.to_parquet(out_path, index=False)

    # ---- Headline: base funding vs deaths burden ----
    head = long[(long["funding_case"] == "base") &
                (long["burden_metric"] == "deaths")].set_index("pathogen")
    head = head.sort_values("misalignment")
    print("\n=== HEADLINE: base funding vs attributable deaths (2019) ===")
    pct = (head[["funding_share", "burden_share", "misalignment"]] * 100).round(1)
    pct.columns = ["funding %", "burden %", "f-b pp"]
    print(pct.to_string())
    print(f"\nGMI (half-sum-abs misalignment) = {summary[('base','deaths')]*100:.1f}%  "
          "(0% = perfect alignment; 100% = fully misaligned)")

    print("\n=== GMI under all variants ===")
    grid = pd.DataFrame(
        {bm: {fc: summary[(fc, bm)] * 100 for fc in fund_cases}
         for bm in burden_metrics}).round(1)
    grid.index.name = "funding_case"
    print(grid.to_string(), "  (units: percentage points / 2 of total share)")

    print(f"\nwrote {out_path.relative_to(REPO_ROOT)}  ({len(long)} rows: "
          "6 pathogens x 3 funding cases x 2 burden metrics)")
    print("\nUNDER = misalignment < 0 (burden share exceeds funding share);")
    print("OVER  = misalignment > 0. Magnitudes in percentage points.")


if __name__ == "__main__":
    main()
