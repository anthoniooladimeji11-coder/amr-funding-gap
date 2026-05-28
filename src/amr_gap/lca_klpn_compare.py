"""
Klebsiella K-selection helper (Step 2 decision support).

Fits the Klebsiella LCA at K=4, 5, 6 and prints, for each:
  - fit metrics (BIC, AIC, assignment entropy, mean max-prob),
  - class sizes,
  - class profiles: P(non-susceptible | class) per drug.

Use this to decide the number of phenotypes by PARSIMONY + INTERPRETABILITY,
not BIC alone: does each extra class reveal a genuinely new, nameable
phenotype, or just split an existing one into trivial variants?

Run:
    PYTHONPATH=src python -m amr_gap.lca_klpn_compare
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from stepmix.stepmix import StepMix

REPO_ROOT = Path(__file__).resolve().parents[2]
SIR = REPO_ROOT / "data" / "interim" / "core_pathogens_sir.parquet"

RANDOM_STATE = 42
N_INIT = 10
ORG = "klebsiella pneumoniae"
PANEL = ["meropenem", "imipenem", "doripenem", "ertapenem", "amikacin",
         "gentamicin", "ciprofloxacin", "levofloxacin", "aztreonam", "colistin"]
MIN_DRUGS = 4


def build_matrix(df):
    g = df[(df["organism"] == ORG) & (df["antibiotic"].isin(PANEL))].copy()
    g["ns"] = g["sir_eucast"].map({"S": 0, "I": 1, "R": 1})
    g = g[g["ns"].notna()]
    wide = g.pivot_table(index="isolate_id", columns="antibiotic",
                         values="ns", aggfunc="max").reindex(columns=PANEL)
    return wide[wide.notna().sum(axis=1) >= MIN_DRUGS]


def main():
    df = pd.read_parquet(SIR)
    wide = build_matrix(df)
    X = wide.to_numpy(dtype=float)
    print(f"Klebsiella isolates: {len(wide):,}  | panel: {len(PANEL)} drugs\n")

    for k in (4, 5, 6):
        m = StepMix(n_components=k, measurement="binary_nan", n_init=N_INIT,
                    random_state=RANDOM_STATE, max_iter=300, verbose=0)
        m.fit(X)
        proba = m.predict_proba(X)
        ent = float((-(np.clip(proba, 1e-12, 1) *
                       np.log(np.clip(proba, 1e-12, 1))).sum(1) /
                     np.log(k)).mean())
        labels = m.predict(X)
        lab = wide.copy()
        lab["c"] = labels
        prof = (lab.groupby("c")[PANEL].mean() * 100).round(0).astype(int)
        sizes = lab["c"].value_counts().sort_index()
        pct = (sizes / len(lab) * 100).round(1)

        print("=" * 72)
        print(f"K = {k}   BIC={m.bic(X):,.0f}  AIC={m.aic(X):,.0f}  "
              f"entropy={ent:.3f}  mean_maxprob={proba.max(1).mean():.3f}")
        print("=" * 72)
        sztbl = pd.DataFrame({"n": sizes, "pct": pct})
        print("class sizes:\n" + sztbl.to_string())
        print("\nP(non-susceptible | class), %:")
        print(prof.to_string())
        print()

    print("DECISION GUIDE: pick the smallest K where every class is a distinct,")
    print("nameable phenotype (e.g. susceptible / ESBL-FQ / carbapenem-R / XDR).")
    print("If K=6's extra classes are just minor splits, prefer K=4 or K=5.")


if __name__ == "__main__":
    main()
