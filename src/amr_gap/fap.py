"""
Step 4a - Funding Allocation Profile (FAP), the funding half of the GMI.

Extended to the SIX leading GRAM pathogens. Builds a pathogen x year table of
R&D funding (USD) attributed via the Hub 'Categories' coded field, in three
columns:
  - base_usd               : funding from projects naming the genus specifically
  - sens_proportional_usd  : base + this genus's PROPORTIONAL share of its
                             Gram-class broad pool (proportional to base shares)
  - sens_equal_usd         : base + an EQUAL share of its Gram-class broad pool
                             (split evenly across target genera in that class)

KEY DESIGN: Gram-negative and Gram-positive broad pools are kept SEPARATE.
Gram-negative-targeted funding redistributes only among Gram-negative genera;
Gram-positive only among Gram-positive. Pooling them would let
Staphylococcus-targeted money redistribute onto Acinetobacter, which is wrong.

Decisions encoded (DECISIONS.md D-005, D-015):
  - Attribution via Categories lines "Infectious Agent / Bacteria / Gram <X> / <Genus> spp."
  - Genus level: e.g. Klebsiella spp. ~ K. pneumoniae (genus->species proxy, stated).
    Streptococcus spp. is the loosest fit (genus includes non-pneumoniae spp.).
  - Year = Start Year. Multi-genus projects split USD equally among named genera.
  - BASE CASE is headline; sens_* columns are a robustness band.

Run:
    PYTHONPATH=src python -m amr_gap.fap
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
HUB = REPO_ROOT / "data" / "raw" / "Projects.xlsx"
PROC = REPO_ROOT / "data" / "processed"

# Genera we model and the Gram class each belongs to. The key is the lowercase
# search token used inside Categories agent lines; the value is the Gram class
# that drives which broad pool the genus inherits from.
TARGETS = {
    "escherichia":     "negative",
    "klebsiella":      "negative",
    "acinetobacter":   "negative",
    "pseudomonas":     "negative",
    "staphylococcus":  "positive",
    "streptococcus":   "positive",
}


def agent_lines(cell: str) -> list[str]:
    return [ln.strip() for ln in str(cell).split("\n")
            if "infectious agent" in ln.lower()]


def named_genera(cell: str) -> list[str]:
    """Which TARGET genera are named on an Infectious Agent line in this project."""
    lines = agent_lines(cell)
    out = []
    for key in TARGETS:
        if any(key in ln.lower() for ln in lines):
            out.append(key)
    return out


def has_broad_gram(cell: str, gram: str) -> bool:
    """True if any Infectious Agent line contains 'gram <gram>' (case-insensitive)."""
    return any(f"gram {gram}" in ln.lower() for ln in agent_lines(cell))


def main() -> None:
    df = pd.read_excel(HUB, sheet_name="data")
    df["amt"] = pd.to_numeric(df["Amount USD"], errors="coerce").fillna(0.0)
    df["yr"] = pd.to_numeric(df["Start Year"], errors="coerce").astype("Int64")
    df["cats"] = df["Categories"].astype("string").fillna("")

    df["genera"] = df["cats"].map(named_genera)
    df["n_named"] = df["genera"].map(len)
    df["broad_neg"] = df["cats"].map(lambda c: has_broad_gram(c, "negative"))
    df["broad_pos"] = df["cats"].map(lambda c: has_broad_gram(c, "positive"))

    # ----- BASE CASE: per-genus funding; split USD equally if multiple named. ---
    base_rows = []
    for _, r in df[df["n_named"] > 0].iterrows():
        share = r["amt"] / r["n_named"]
        for g in r["genera"]:
            base_rows.append((g, r["yr"], share, 1))
    base = pd.DataFrame(base_rows, columns=["pathogen", "year", "usd", "proj"])
    base = base.groupby(["pathogen", "year"], dropna=True).agg(
        base_usd=("usd", "sum"), n_projects=("proj", "sum")).reset_index()
    # add Gram class for downstream pool routing
    base["gram"] = base["pathogen"].map(TARGETS)

    # ----- BROAD POOLS by year, per Gram class. A project enters its Gram-class
    # broad pool iff (a) it has the Gram-class line AND (b) it does NOT name any
    # target genus (avoid double counting). -----
    pool_by_year = {}
    for gram, col in (("negative", "broad_neg"), ("positive", "broad_pos")):
        mask = df[col] & (df["n_named"] == 0)
        pool_by_year[gram] = (df[mask].groupby("yr")["amt"].sum()
                              .rename(f"broad_{gram}_usd").reset_index()
                              .rename(columns={"yr": "year"}))
        print(f"  broad Gram-{gram} pool: {int(mask.sum()):,} projects, "
              f"${df[mask]['amt'].sum()/1e6:,.1f}M")

    # ----- Sensitivity allocations (per Gram class) -----
    fap = base.copy()
    # attach this row's Gram-class broad pool for its year
    for gram in ("negative", "positive"):
        fap = fap.merge(pool_by_year[gram], on="year", how="left")
    fap["broad_negative_usd"] = fap["broad_negative_usd"].fillna(0.0)
    fap["broad_positive_usd"] = fap["broad_positive_usd"].fillna(0.0)
    # each row gets only its own Gram class's pool
    fap["broad_for_row"] = np.where(fap["gram"] == "negative",
                                    fap["broad_negative_usd"],
                                    fap["broad_positive_usd"])

    # proportional split: within each (year, gram), share of base * pool
    grp = ["year", "gram"]
    yr_gram_base = fap.groupby(grp)["base_usd"].transform("sum")
    prop_share = np.where(yr_gram_base > 0, fap["base_usd"] / yr_gram_base, 0.0)
    fap["sens_proportional_usd"] = fap["base_usd"] + fap["broad_for_row"] * prop_share

    # equal split: divide the pool by the number of target genera in this Gram class
    n_target_in_class = fap.groupby(grp)["pathogen"].transform("nunique")
    fap["sens_equal_usd"] = (fap["base_usd"]
                             + fap["broad_for_row"] / n_target_in_class)

    fap = (fap.drop(columns=["broad_negative_usd", "broad_positive_usd",
                             "broad_for_row"])
              .sort_values(["pathogen", "year"]))

    PROC.mkdir(parents=True, exist_ok=True)
    out = PROC / "fap_pathogen_year.parquet"
    fap.to_parquet(out, index=False)

    # ----- Reporting -----
    print("\n=== FAP summary (totals across all years, USD millions) ===")
    tot = fap.groupby(["gram", "pathogen"])[
        ["base_usd", "sens_proportional_usd", "sens_equal_usd"]].sum()
    print((tot / 1e6).round(1).to_string())

    print(f"\nyears covered: {int(fap['year'].min())}-{int(fap['year'].max())}")
    print("\nrows per pathogen:")
    print(fap.groupby("pathogen").size().to_string())
    print(f"\nwrote {out.relative_to(REPO_ROOT)}  ({len(fap)} pathogen-year rows)")
    print("\nNOTE: base_usd is the headline; sens_* are a robustness band.")
    print("Streptococcus spp. genus->species proxy is the loosest; flag in writeup.")


if __name__ == "__main__":
    main()
