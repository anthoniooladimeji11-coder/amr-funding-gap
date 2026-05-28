"""
Step 2 - Latent Class Analysis (multidrug-resistance phenotyping).

Fits, separately for each pathogen, a latent class model over the binary
non-susceptibility matrix (S=0, I/R=1, untested=missing). Selects the number of
classes K in [2..6] by BIC, reports class profiles (per-class probability of
non-susceptibility to each drug) and assignment certainty (entropy), and writes
each isolate's most-likely class back to disk.

Design decisions (see DECISIONS.md):
  - Separate model per pathogen (different panels, different biology).
  - Binary encoding (I collapsed into non-susceptible).
  - Variation-based clustering panels; min-drugs thresholds Ab=3, Kp=4.
  - Missing data handled by stepmix (full-information ML), not dropped.
  - Fixed random_state + multiple inits for reproducibility.

Run:
    PYTHONPATH=src python -m amr_gap.lca                 # both pathogens, K=2..6
    PYTHONPATH=src python -m amr_gap.lca --kmax 7
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from stepmix.stepmix import StepMix

REPO_ROOT = Path(__file__).resolve().parents[2]
SIR = REPO_ROOT / "data" / "interim" / "core_pathogens_sir.parquet"
OUT_DIR = REPO_ROOT / "data" / "interim"
PROC_DIR = REPO_ROOT / "data" / "processed"

RANDOM_STATE = 42
N_INIT = 10  # multiple initializations -> stable solution

PANELS = {
    "acinetobacter baumannii": (
        ["amikacin", "gentamicin", "meropenem", "imipenem", "levofloxacin",
         "colistin"],
        3,
    ),
    "klebsiella pneumoniae": (
        ["meropenem", "imipenem", "doripenem", "ertapenem", "amikacin",
         "gentamicin", "ciprofloxacin", "levofloxacin", "aztreonam", "colistin"],
        4,
    ),
}


def build_matrix(df: pd.DataFrame, organism: str, panel: list[str],
                 min_drugs: int) -> pd.DataFrame:
    g = df[(df["organism"] == organism) & (df["antibiotic"].isin(panel))].copy()
    m = {"S": 0, "I": 1, "R": 1}
    g["ns"] = g["sir_eucast"].map(m)
    g = g[g["ns"].notna()]
    wide = g.pivot_table(index="isolate_id", columns="antibiotic",
                         values="ns", aggfunc="max")
    wide = wide.reindex(columns=panel)
    wide = wide[wide.notna().sum(axis=1) >= min_drugs]
    return wide


def entropy_of_assignment(proba: np.ndarray) -> float:
    """Normalized entropy of posterior class probabilities (0=certain,1=uniform).

    Returns mean across isolates of 1 - (relative entropy). Higher 'certainty'
    = better-separated classes. We report mean posterior max-prob too.
    """
    eps = 1e-12
    p = np.clip(proba, eps, 1)
    k = p.shape[1]
    ent = -(p * np.log(p)).sum(axis=1) / np.log(k)  # in [0,1], 0=certain
    return float(ent.mean())


def fit_pathogen(wide: pd.DataFrame, kmax: int) -> tuple[int, pd.DataFrame, dict]:
    """Fit K=2..kmax, pick best by BIC. Returns (best_k, labeled_df, report)."""
    # stepmix expects a numpy array; missing as np.nan is supported.
    X = wide.to_numpy(dtype=float)

    results = {}
    for k in range(2, kmax + 1):
        model = StepMix(
            n_components=k,
            measurement="binary_nan",   # binary indicators with missing values
            n_init=N_INIT,
            random_state=RANDOM_STATE,
            max_iter=300,
            verbose=0,
        )
        model.fit(X)
        bic = model.bic(X)
        aic = model.aic(X)
        proba = model.predict_proba(X)
        ent = entropy_of_assignment(proba)
        mean_maxp = float(proba.max(axis=1).mean())
        results[k] = dict(model=model, bic=bic, aic=aic,
                          entropy=ent, mean_maxp=mean_maxp)
        print(f"    K={k}:  BIC={bic:,.0f}  AIC={aic:,.0f}  "
              f"assign_entropy={ent:.3f}  mean_maxprob={mean_maxp:.3f}")

    best_k = min(results, key=lambda k: results[k]["bic"])
    best = results[best_k]
    print(f"  -> best K by BIC = {best_k}")

    # Class profiles: P(non-susceptible | class) per drug.
    model = best["model"]
    labels = model.predict(model_X := wide.to_numpy(dtype=float))
    proba = model.predict_proba(model_X)

    labeled = wide.copy()
    labeled["lca_class"] = labels
    labeled["lca_maxprob"] = proba.max(axis=1)

    # Profile = mean NS per drug within each assigned class (observed),
    # which is interpretable and robust to report alongside model params.
    profile = labeled.groupby("lca_class")[wide.columns.tolist()].mean()
    sizes = labeled["lca_class"].value_counts().sort_index()

    report = dict(results=results, best_k=best_k, profile=profile, sizes=sizes)
    return best_k, labeled, report


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--kmax", type=int, default=6)
    args = p.parse_args()

    df = pd.read_parquet(SIR)
    PROC_DIR.mkdir(parents=True, exist_ok=True)

    for organism, (panel, min_drugs) in PANELS.items():
        short = organism.split()[0][:2] + organism.split()[1][:2]  # e.g. 'acba'
        print("=" * 64)
        print(f"LCA: {organism}  (panel={len(panel)}, min_drugs={min_drugs})")
        print("=" * 64)

        wide = build_matrix(df, organism, panel, min_drugs)
        print(f"  isolates: {len(wide):,}  | drugs: {len(panel)}")
        print(f"  fitting K=2..{args.kmax} ...")

        best_k, labeled, report = fit_pathogen(wide, args.kmax)

        print(f"\n  Class sizes (K={best_k}):")
        print(report["sizes"].to_string())
        print(f"\n  Class profiles  P(non-susceptible | class), by drug:")
        print((report["profile"] * 100).round(1).to_string())

        out = PROC_DIR / f"lca_{short}.parquet"
        labeled.reset_index().to_parquet(out, index=False)
        print(f"\n  wrote {out.relative_to(REPO_ROOT)}  "
              f"({len(labeled):,} isolates labeled)\n")

    print("Compare profiles to the co-resistance preview: do classes match the "
          "expected phenotypes (susceptible / FQ-ESBL / carbapenem-R / XDR)?")


if __name__ == "__main__":
    main()
