"""
Step 4a - Funding Allocation Profile (FAP), the funding half of the GMI.

Builds a pathogen x year table of R&D funding (USD) attributed via the Hub
'Categories' coded field, in three columns:
  - base_usd               : funding from projects naming the genus specifically
  - sens_proportional_usd  : base + this genus's PROPORTIONAL share of the broad
                             'Gram negative' pool (proportional to base shares)
  - sens_equal_usd         : base + an EQUAL share of the broad pool (split evenly
                             across all Gram-negative genera present)

Decisions encoded (see DECISIONS.md D-005, D-015):
  - Attribution via Categories lines "... / Gram negative / <Genus> spp."
  - Genus level: 'Klebsiella spp.' ~ K. pneumoniae; 'Acinetobacter spp.' ~ A. baumannii
    (genus->species proxy, declared as an assumption).
  - Year = Start Year. Multi-pathogen projects split USD equally among named genera.
  - BASE CASE is the headline; the two sensitivity variants are a robustness band
    around the arbitrary redistribution of the broad pool.

Run:
    PYTHONPATH=src python -m amr_gap.fap
"""

from __future__ import annotations

from pathlib import Path
import re
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
HUB = REPO_ROOT / "data" / "raw" / "Projects.xlsx"
PROC = REPO_ROOT / "data" / "processed"

# Genera we model (key) and the label they carry in the Categories field (value).
TARGET_GENERA = {
    "klebsiella": "Klebsiella spp.",
    "acinetobacter": "Acinetobacter spp.",
}
# Broad pool marker: a Categories agent line that names the Gram-negative group
# but NOT a specific genus (i.e. ends at 'Gram negative' or 'Other Gram negative').
BROAD_MARKERS = ["gram negative / other gram negative",
                 "/ bacteria / gram negative"]  # line ending at the group level


def agent_lines(cell: str) -> list[str]:
    return [ln.strip() for ln in str(cell).split("\n")
            if "infectious agent" in ln.lower()]


def named_genera(cell: str) -> list[str]:
    """Return which TARGET genera are named in this project's Categories."""
    low = str(cell).lower()
    out = []
    for key in TARGET_GENERA:
        if key in low and "infectious agent" in low:
            # ensure the genus appears on an infectious-agent line
            if any(key in ln.lower() for ln in agent_lines(cell)):
                out.append(key)
    return out


def is_broad_gramneg(cell: str) -> bool:
    """True if project targets Gram-negatives broadly (group-level, no specific
    target genus among ours). Used for the sensitivity pool."""
    lines = [ln.lower() for ln in agent_lines(cell)]
    has_gn = any("gram negative" in ln for ln in lines)
    return has_gn


def main() -> None:
    df = pd.read_excel(HUB, sheet_name="data")
    df["amt"] = pd.to_numeric(df["Amount USD"], errors="coerce").fillna(0.0)
    df["yr"] = pd.to_numeric(df["Start Year"], errors="coerce").astype("Int64")
    df["cats"] = df["Categories"].astype("string").fillna("")

    df["genera"] = df["cats"].map(named_genera)
    df["n_named"] = df["genera"].map(len)
    df["broad"] = df["cats"].map(is_broad_gramneg)

    # ----- BASE CASE: projects naming a target genus; split USD equally if it
    # names more than one of our genera. -----
    base_rows = []
    for _, r in df[df["n_named"] > 0].iterrows():
        share = r["amt"] / r["n_named"]
        for g in r["genera"]:
            base_rows.append((g, r["yr"], share, 1))
    base = pd.DataFrame(base_rows, columns=["pathogen", "year", "usd", "proj"])
    base = base.groupby(["pathogen", "year"], dropna=True).agg(
        base_usd=("usd", "sum"), n_projects=("proj", "sum")).reset_index()

    # ----- BROAD POOL by year: projects that target Gram-negatives broadly but
    # do NOT name one of our genera (avoid double counting). -----
    broad_mask = df["broad"] & (df["n_named"] == 0)
    broad_by_year = (df[broad_mask].groupby("yr")["amt"].sum()
                     .rename("broad_usd").reset_index()
                     .rename(columns={"yr": "year"}))

    # ----- Sensitivity allocations -----
    # base shares per pathogen-year (for proportional split)
    fap = base.copy()
    fap = fap.merge(broad_by_year, on="year", how="left")
    fap["broad_usd"] = fap["broad_usd"].fillna(0.0)

    # proportional: each genus gets broad_usd * (its base / total base that year)
    yr_base_total = fap.groupby("year")["base_usd"].transform("sum")
    prop_share = np.where(yr_base_total > 0, fap["base_usd"] / yr_base_total, 0.0)
    fap["sens_proportional_usd"] = fap["base_usd"] + fap["broad_usd"] * prop_share

    # equal: broad_usd split evenly across the number of target genera present
    n_genera_year = fap.groupby("year")["pathogen"].transform("nunique")
    fap["sens_equal_usd"] = fap["base_usd"] + fap["broad_usd"] / n_genera_year

    fap = fap.drop(columns=["broad_usd"]).sort_values(["pathogen", "year"])

    PROC.mkdir(parents=True, exist_ok=True)
    out = PROC / "fap_pathogen_year.parquet"
    fap.to_parquet(out, index=False)

    # ----- Reporting -----
    print("=== FAP summary (totals across all years) ===")
    tot = fap.groupby("pathogen")[
        ["base_usd", "sens_proportional_usd", "sens_equal_usd"]].sum()
    print((tot / 1e6).round(1).to_string(), "  (USD millions)")
    print(f"\nbroad Gram-negative pool (no target genus named), total: "
          f"${df[broad_mask]['amt'].sum()/1e6:,.1f}M across "
          f"{int(broad_mask.sum()):,} projects")
    print(f"\nyears covered: {int(fap['year'].min())}–{int(fap['year'].max())}")
    print("\nrows per pathogen:")
    print(fap.groupby("pathogen").size().to_string())
    print(f"\nwrote {out.relative_to(REPO_ROOT)}  ({len(fap)} pathogen-year rows)")
    print("\nNOTE: base_usd is the headline; the two sens_* columns are a "
          "robustness band around the arbitrary broad-pool redistribution.")


if __name__ == "__main__":
    main()
