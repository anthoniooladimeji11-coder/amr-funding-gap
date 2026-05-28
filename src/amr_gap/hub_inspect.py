"""
Hub funding data inspection (Step 4 planning).

Answers the questions that determine GMI granularity:
  1. Do our two pathogens (K. pneumoniae, A. baumannii) appear in the agent
     fields, and how are they labelled?
  2. Is there any field resembling DRUG CLASS / antibiotic target? (We expect
     not -> GMI funding side likely pathogen-level only.)
  3. How much funding is attributable to a SPECIFIC pathogen vs broad / platform
     / multi-pathogen research? (drives the 'unattributable bucket' rule.)
  4. Year coverage and total amount, for the 5-year lag and FAP shares.

Run:
    PYTHONPATH=src python -m amr_gap.hub_inspect
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
HUB = REPO_ROOT / "data" / "raw" / "Projects.xlsx"


def show_values(df, col, n=20):
    if col not in df.columns:
        print(f"  [{col}] NOT PRESENT")
        return
    vals = df[col].dropna()
    uniq = vals.unique()
    print(f"  [{col}] {len(uniq)} unique, {df[col].isna().sum():,} missing")
    for v in list(uniq)[:n]:
        print("      ", repr(v)[:90])
    if len(uniq) > n:
        print(f"       ... (+{len(uniq)-n} more)")


def main() -> None:
    df = pd.read_excel(HUB, sheet_name="data")
    print(f"Hub projects: {len(df):,} rows, {len(df.columns)} cols\n")

    print("=== AMOUNTS & YEARS ===")
    for c in ["Amount USD", "Amount EUR", "Start Year", "End Year"]:
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            print(f"  {c}: non-null={s.notna().sum():,}  "
                  f"sum={s.sum():,.0f}  min={s.min()}  max={s.max()}")
    print()

    print("=== PATHOGEN ATTRIBUTION FIELDS ===")
    for c in ["Infectious Agent", "Individual Infectious Agent", "Disease"]:
        show_values(df, c, n=25)
        print()

    print("=== SEARCH FOR OUR TWO PATHOGENS (any agent field) ===")
    agent_cols = [c for c in ["Infectious Agent", "Individual Infectious Agent",
                              "Disease", "Title", "Abstract"] if c in df.columns]
    for patho in ["klebsiella", "acinetobacter"]:
        hits = pd.Series(False, index=df.index)
        for c in agent_cols:
            hits = hits | df[c].astype("string").str.lower().str.contains(
                patho, na=False)
        amt = pd.to_numeric(df.get("Amount USD"), errors="coerce")
        print(f"  '{patho}': {hits.sum():,} projects mention it; "
              f"funding(USD)={amt[hits].sum():,.0f}")
    print()

    print("=== RESEARCH AREA / SECTOR / PRODUCT (is there a DRUG-CLASS field?) ===")
    for c in ["Sector", "Sector Subcategory", "Research Area",
              "Research Subcategory", "Project Group", "Product Name",
              "Categories"]:
        show_values(df, c, n=20)
        print()

    print("DONE. Key reads: (a) which agent field names our pathogens cleanly;")
    print("(b) whether any column encodes drug/antibiotic class; (c) how much")
    print("funding is pathogen-specific vs broad -> sets GMI granularity + bucket.")


if __name__ == "__main__":
    main()
