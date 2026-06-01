"""
Step 3b - Fit the random-walk Beta-binomial Stan model and project XDR
prevalence to 2030 for K. pneumoniae and A. baumannii.

Fits separately per pathogen using cmdstanpy on year-only data through 2024
(2025 excluded as incomplete reporting). Projects 2025-2030.

Writes per-pathogen posterior samples:
    data/processed/xdr_projection_{acba|klpn}.parquet
with columns:
    year, observed_pct, p_mean, p_lo80, p_hi80, p_lo95, p_hi95, kind
where kind in {historical, projected}.

Also writes a combined summary:
    data/processed/xdr_projection_summary.csv

Run:
    PYTHONPATH=src python -m amr_gap.xdr_project
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from cmdstanpy import CmdStanModel

REPO_ROOT = Path(__file__).resolve().parents[2]
PROC = REPO_ROOT / "data" / "processed"
STAN_FILE = REPO_ROOT / "src" / "amr_gap" / "stan" / "xdr_rw.stan"

PATHOGENS = {
    "acba": "acinetobacter baumannii",
    "klpn": "klebsiella pneumoniae",
}

FIT_YEAR_END = 2024   # last year used in fit; 2025 excluded as incomplete
HORIZON_END  = 2030   # project to this year inclusive


def fit_one(short: str, organism: str, model: CmdStanModel) -> pd.DataFrame:
    ts = pd.read_parquet(PROC / f"xdr_timeseries_{short}_yearonly.parquet")
    ts["year"] = ts["year"].astype(int)
    # Restrict to complete years; exclude 2025 (partial reporting).
    obs = ts[ts["year"] <= FIT_YEAR_END].sort_values("year").reset_index(drop=True)
    years_hist = obs["year"].to_list()
    n = obs["n_total"].astype(int).to_list()
    y = obs["n_xdr"].astype(int).to_list()

    last = years_hist[-1]
    years_pred = list(range(last + 1, HORIZON_END + 1))
    T = len(years_hist)
    T_pred = len(years_pred)

    print(f"\n=== {organism} ===")
    print(f"  fitting on {T} years ({years_hist[0]}-{last}); "
          f"projecting {T_pred} years ({years_pred[0]}-{years_pred[-1]})")

    fit = model.sample(
        data={"T": T, "T_pred": T_pred, "n": n, "y": y},
        chains=4, parallel_chains=4,
        iter_warmup=1000, iter_sampling=2000,
        seed=42, show_progress=False,
    )

    # Diagnostics
    diag = fit.diagnose()
    if "no problems detected" not in diag.lower():
        print("  Stan diagnostics flagged issues:\n", diag)
    else:
        print("  Stan diagnostics: no problems detected.")

    draws = fit.draws_pd()
    print(f"  posterior draws: {len(draws):,}  | sigma (RW step SD on logit): "
          f"mean={draws['sigma'].mean():.3f}  sd={draws['sigma'].std():.3f}")

    # Assemble per-year posterior summaries.
    rows = []
    for i, yr in enumerate(years_hist, start=1):
        col = f"p_hist[{i}]"
        s = draws[col]
        rows.append(_summarize(yr, "historical", s, obs.iloc[i-1]))
    for i, yr in enumerate(years_pred, start=1):
        col = f"p_pred[{i}]"
        s = draws[col]
        rows.append(_summarize(yr, "projected", s, None))

    out = pd.DataFrame(rows)
    out_path = PROC / f"xdr_projection_{short}.parquet"
    out.to_parquet(out_path, index=False)
    print(f"  wrote {out_path.relative_to(REPO_ROOT)}")

    # Also save the raw posterior draws of p_pred for downstream uncertainty
    # propagation (used by the projected GMI). One column per projected year.
    pred_cols = [f"p_pred[{i}]" for i in range(1, T_pred + 1)]
    pred_draws = draws[pred_cols].copy()
    pred_draws.columns = [str(y) for y in years_pred]
    draws_path = PROC / f"xdr_projection_{short}_draws.parquet"
    pred_draws.to_parquet(draws_path, index=False)
    print(f"  wrote {draws_path.relative_to(REPO_ROOT)} "
          f"({len(pred_draws):,} draws x {T_pred} years)")

    # Show the headline result row by row
    show = out[["year", "kind", "observed_pct", "p_mean", "p_lo80", "p_hi80"]].copy()
    for c in ("observed_pct", "p_mean", "p_lo80", "p_hi80"):
        show[c] = (show[c] * 100).round(1) if c != "observed_pct" else show[c].round(1)
    print("\n  year | kind        | obs%  | proj % (80% UI)")
    for _, r in show.iterrows():
        obs_str = f"{r['observed_pct']:>5.1f}" if pd.notna(r["observed_pct"]) else "  -  "
        print(f"  {int(r['year'])} | {r['kind']:<11s} | {obs_str} | "
              f"{r['p_mean']:>4.1f}  ({r['p_lo80']:>4.1f}, {r['p_hi80']:>4.1f})")

    out["pathogen"] = organism
    return out


def _summarize(year, kind, s, obs_row):
    p_mean = float(s.mean())
    p_lo80, p_hi80 = (float(np.quantile(s, q)) for q in (0.10, 0.90))
    p_lo95, p_hi95 = (float(np.quantile(s, q)) for q in (0.025, 0.975))
    return {
        "year": int(year),
        "kind": kind,
        "observed_pct": (float(obs_row["n_xdr"] / obs_row["n_total"] * 100)
                         if obs_row is not None else np.nan),
        "p_mean": p_mean,
        "p_lo80": p_lo80, "p_hi80": p_hi80,
        "p_lo95": p_lo95, "p_hi95": p_hi95,
    }


def main():
    model = CmdStanModel(stan_file=STAN_FILE)
    print(f"compiled: {STAN_FILE.relative_to(REPO_ROOT)}")

    all_results = []
    for short, organism in PATHOGENS.items():
        all_results.append(fit_one(short, organism, model))

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(PROC / "xdr_projection_summary.csv", index=False)
    print(f"\nwrote {PROC / 'xdr_projection_summary.csv'}")
    print("\nDone. Next: feed the 2030 projected prevalences into the projected GMI.")


if __name__ == "__main__":
    main()
