"""
Step 4c - Pipeline cross-flag for the GMI.

For each pathogen, counts late-stage therapeutic candidates the Hub is funding
('Therapeutics / Development' OR 'Therapeutics / Approval & Post approval' in
the Categories coded field). Merges with the GMI table to produce, per pathogen:
  - funding share (from FAP)
  - burden share (from burden weights)
  - misalignment (funding - burden, percentage points)
  - n_pipeline_projects: late-stage therapeutic projects naming this pathogen
  - pipeline_share: share of total late-stage pipeline across the 6 pathogens
  - QUADRANT label combining funding-vs-burden with pipeline-thin-vs-rich.

Definitions (D-018 - see DECISIONS.md):
  - 'Late-stage' = Development + Approval (WHO/PEW pipeline standard).
  - 'Pipeline-thin' vs 'pipeline-rich': split at the median pipeline_share.
  - 'Under/overfunded' relative to burden: split at zero on misalignment.
  - QUADRANTS:
      NEGLECTED                  : underfunded AND pipeline-thin
      SERVED DESPITE LOW FUNDING : underfunded AND pipeline-rich
      INVESTMENT-TRANSLATION GAP : overfunded AND pipeline-thin
      WELL-RESOURCED             : overfunded AND pipeline-rich

Caveats (stated in writeup):
  - Project-count proxies candidate-count (one drug can have multiple projects).
    Affects all pathogens similarly -> fine for relative comparison.
  - Multi-pathogen projects: counted once per named target genus (matches the
    FAP attribution rule).

Run:
    PYTHONPATH=src python -m amr_gap.pipeline
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
HUB = REPO_ROOT / "data" / "raw" / "Projects.xlsx"
PROC = REPO_ROOT / "data" / "processed"
GMI_PATH = PROC / "gmi.parquet"

# Same six genera as the FAP/GMI.
TARGETS = ["escherichia", "klebsiella", "acinetobacter", "pseudomonas",
           "staphylococcus", "streptococcus"]

GENUS_TO_SPECIES = {
    "escherichia":    "escherichia coli",
    "klebsiella":     "klebsiella pneumoniae",
    "acinetobacter":  "acinetobacter baumannii",
    "pseudomonas":    "pseudomonas aeruginosa",
    "staphylococcus": "staphylococcus aureus",
    "streptococcus":  "streptococcus pneumoniae",
}

# Late-stage tokens we look for in the Categories field. These match the
# Research Subcategory coded lines: "3200 Research Area / Therapeutics /
# Development" and "3300 Research Area / Therapeutics / Approval & Post approval".
LATE_STAGE_TOKENS = (
    "therapeutics / development",
    "therapeutics / approval & post approval",
    "therapeutics / approval and post approval",  # spelling variant safety
)


def agent_lines(cell: str) -> list[str]:
    return [ln.strip() for ln in str(cell).split("\n")
            if "infectious agent" in ln.lower()]


def named_genera(cell: str) -> list[str]:
    lines = agent_lines(cell)
    return [g for g in TARGETS if any(g in ln.lower() for ln in lines)]


def is_late_stage(cell: str) -> bool:
    low = str(cell).lower()
    return any(tok in low for tok in LATE_STAGE_TOKENS)


def main() -> None:
    # ---- Hub: count late-stage therapeutic projects per pathogen ----
    df = pd.read_excel(HUB, sheet_name="data")
    df["cats"] = df["Categories"].astype("string").fillna("")
    df["genera"] = df["cats"].map(named_genera)
    df["late_stage"] = df["cats"].map(is_late_stage)

    # quick visibility on the late-stage filter
    n_late = int(df["late_stage"].sum())
    n_total = len(df)
    print(f"Hub projects: {n_total:,}; late-stage therapeutics: "
          f"{n_late:,} ({100*n_late/n_total:.1f}%)")

    # explode named genera and count late-stage among them
    rows = []
    for _, r in df[df["late_stage"]].iterrows():
        for g in r["genera"]:
            rows.append(g)
    if not rows:
        raise SystemExit("No late-stage projects matched any target genus - "
                         "check token spellings against the Categories field.")
    pipe = (pd.Series(rows, name="genus").value_counts()
            .rename("n_pipeline_projects").reset_index())
    pipe["species"] = pipe["genus"].map(GENUS_TO_SPECIES)

    # Ensure every target appears (zero if absent)
    for g in TARGETS:
        if g not in pipe["genus"].values:
            pipe = pd.concat([pipe, pd.DataFrame(
                [{"genus": g, "n_pipeline_projects": 0,
                  "species": GENUS_TO_SPECIES[g]}])], ignore_index=True)
    pipe["pipeline_share"] = pipe["n_pipeline_projects"] / pipe["n_pipeline_projects"].sum()

    # ---- GMI headline row (base funding vs deaths) ----
    gmi = pd.read_parquet(GMI_PATH)
    head = gmi[(gmi["funding_case"] == "base") &
               (gmi["burden_metric"] == "deaths")].copy()

    out = head.merge(pipe[["species", "n_pipeline_projects", "pipeline_share"]],
                     left_on="pathogen", right_on="species", how="left").drop(
        columns=["species"])

    # Quadrant labels: split pipeline_share at median; misalignment at 0
    pip_med = out["pipeline_share"].median()
    def quadrant(r):
        under = r["misalignment"] < 0
        pipe_rich = r["pipeline_share"] > pip_med
        if under and not pipe_rich: return "NEGLECTED"
        if under and pipe_rich:     return "SERVED DESPITE LOW FUNDING"
        if not under and not pipe_rich: return "INVESTMENT-TRANSLATION GAP"
        return "WELL-RESOURCED"
    out["quadrant"] = out.apply(quadrant, axis=1)

    PROC.mkdir(parents=True, exist_ok=True)
    out_path = PROC / "pipeline_cross.parquet"
    out.to_parquet(out_path, index=False)

    # ---- Report ----
    print("\n=== Pipeline counts (late-stage therapeutics: Development + Approval) ===")
    show = out[["pathogen", "n_pipeline_projects", "pipeline_share"]].copy()
    show["pipeline_share_pct"] = (show["pipeline_share"] * 100).round(1)
    print(show[["pathogen", "n_pipeline_projects", "pipeline_share_pct"]]
          .sort_values("n_pipeline_projects", ascending=False).to_string(index=False))
    print(f"\n  median pipeline_share (split threshold) = {pip_med*100:.1f}%")

    print("\n=== GMI + Pipeline cross-flag (base funding vs deaths burden) ===")
    cols = ["pathogen", "funding_share", "burden_share", "misalignment",
            "pipeline_share", "quadrant"]
    rep = out[cols].copy()
    rep["funding %"] = (rep["funding_share"] * 100).round(1)
    rep["burden %"]  = (rep["burden_share"] * 100).round(1)
    rep["f-b pp"]    = (rep["misalignment"] * 100).round(1)
    rep["pipeline %"] = (rep["pipeline_share"] * 100).round(1)
    rep = rep[["pathogen", "funding %", "burden %", "f-b pp",
               "pipeline %", "quadrant"]].sort_values("f-b pp")
    print(rep.to_string(index=False))

    print(f"\nwrote {out_path.relative_to(REPO_ROOT)}  ({len(out)} pathogens)")
    print("\nReadings:")
    print("  NEGLECTED                  -> action: invest urgently (underfunded + thin pipeline)")
    print("  SERVED DESPITE LOW FUNDING -> action: maintain, monitor for pipeline attrition")
    print("  INVESTMENT-TRANSLATION GAP -> action: investigate why high spend hasn't yielded candidates")
    print("  WELL-RESOURCED             -> baseline reference category")


if __name__ == "__main__":
    main()
