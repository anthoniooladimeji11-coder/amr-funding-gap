"""
Hub Categories inspection - extended to the 6 GRAM pathogens (Step 4b planning).

Verifies, before extending the FAP, that each of the six leading pathogens has
clean attribution in the Hub Categories coded field, and quantifies the two
broad pools (Gram-negative and Gram-positive) separately.

Run:
    PYTHONPATH=src python -m amr_gap.hub_categories_6
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
HUB = REPO_ROOT / "data" / "raw" / "Projects.xlsx"

# Genera we model and the lowercase search keys used inside Categories lines.
TARGETS = {
    "escherichia": "Escherichia spp.",
    "staphylococcus": "Staphylococcus spp.",
    "klebsiella": "Klebsiella spp.",
    "acinetobacter": "Acinetobacter spp.",
    "streptococcus": "Streptococcus spp.",
    "pseudomonas": "Pseudomonas spp.",
}

# Gram class of each (decides which broad pool it inherits from).
GRAM_CLASS = {
    "escherichia": "negative", "klebsiella": "negative",
    "acinetobacter": "negative", "pseudomonas": "negative",
    "staphylococcus": "positive", "streptococcus": "positive",
}


def agent_lines(cell: str):
    return [ln.strip() for ln in str(cell).split("\n")
            if "infectious agent" in ln.lower()]


def names_target(cell, key):
    return any(key in ln.lower() for ln in agent_lines(cell))


def names_any_target(cell):
    return any(any(k in ln.lower() for k in TARGETS) for ln in agent_lines(cell))


def is_broad(cell, gram):
    """True if Categories has a Gram-<gram> agent line."""
    return any(f"gram {gram}" in ln.lower() for ln in agent_lines(cell))


def main():
    df = pd.read_excel(HUB, sheet_name="data")
    df["cats"] = df["Categories"].astype("string").fillna("")
    df["amt"] = pd.to_numeric(df["Amount USD"], errors="coerce").fillna(0)

    print("=== Per-pathogen attribution via Categories ===")
    rows = []
    for key, label in TARGETS.items():
        m = df["cats"].apply(names_target, key=key)
        rows.append((label, GRAM_CLASS[key], int(m.sum()), df.loc[m, "amt"].sum()))
    out = pd.DataFrame(rows, columns=["genus", "gram", "n_projects", "usd"])
    out["usd_millions"] = (out["usd"] / 1e6).round(1)
    print(out.drop(columns="usd").to_string(index=False))

    print("\n=== Broad pools (Gram-class projects NOT naming any target genus) ===")
    any_target = df["cats"].apply(names_any_target)
    for gram in ("negative", "positive"):
        broad = df["cats"].apply(is_broad, gram=gram) & (~any_target)
        print(f"  Gram-{gram} broad pool: {int(broad.sum()):,} projects, "
              f"${df.loc[broad, 'amt'].sum()/1e6:,.1f}M")

    print("\nNOTE: 'Streptococcus spp.' may be a coarse genus including non-AMR-priority")
    print("species (e.g. Strep agalactiae, Group A). The GRAM burden figure for")
    print("'Streptococcus pneumoniae' is species-specific. The genus->species proxy")
    print("is a STATED ASSUMPTION (consistent with D-005); worth flagging in writeup.")


if __name__ == "__main__":
    main()
