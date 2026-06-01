"""
Step 3c - Projected Global Misalignment Index for 2030.

Combines the Bayesian XDR projection draws (Kp, Ab) with the current FAP funding
shares and GRAM 2019 burden weights to produce a projected GMI with uncertainty
propagated from the projection posterior.

Modelling assumption (Option 1, stated openly in DECISIONS.md D-019):
  For each projected pathogen, attributable burden scales with projected XDR
  prevalence relative to its 2019 prevalence:
       burden(X, 2030) = burden(X, 2019) * [prev(X, 2030) / prev(X, 2019)]
  Pathogens not projected (E. coli, S. aureus, S. pneumoniae, P. aeruginosa)
  hold their 2019 burdens (and so are also stated as an assumption).
  Funding shares are held at current cumulative values (no projected FAP).

Uncertainty propagation:
  For each posterior draw of projected prevalence (Kp_2030, Ab_2030), recompute
  per-pathogen projected burdens, shares, misalignments, and the GMI. This
  produces a posterior over the 2030 GMI from which we report mean + 80% UI.

Outputs:
    data/processed/gmi_projected_2030.parquet  (per-pathogen, per-draw)
    Console: current vs projected GMI table + per-pathogen misalignment changes.

Run:
    PYTHONPATH=src python -m amr_gap.gmi_projected
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
PROC = REPO_ROOT / "data" / "processed"
EXT = REPO_ROOT / "data" / "external"

GMI_PATH = PROC / "gmi.parquet"               # current GMI (from Step 4b)
BURDEN_PATH = EXT / "burden_weights.csv"
KLPN_DRAWS = PROC / "xdr_projection_klpn_draws.parquet"
ACBA_DRAWS = PROC / "xdr_projection_acba_draws.parquet"

# 2019 XDR prevalence anchors (observed) used for the scaling ratio.
# These come from xdr_timeseries_*_yearonly.parquet 2019 rows:
#   klpn 2019: 967 / 6671 = 14.5%
#   acba 2019: 2373 / 3621 = 65.5%
KLPN_PREV_2019 = 967 / 6671
ACBA_PREV_2019 = 2373 / 3621

# Map projected genus -> species name as used in burden_weights and gmi tables
PROJECTED = {
    "klebsiella":    ("klebsiella pneumoniae",  KLPN_DRAWS, KLPN_PREV_2019),
    "acinetobacter": ("acinetobacter baumannii", ACBA_DRAWS, ACBA_PREV_2019),
}

PROJ_YEAR = "2030"


def main():
    # ----- Current GMI (headline: base funding vs attributable deaths) -----
    gmi = pd.read_parquet(GMI_PATH)
    cur = gmi[(gmi["funding_case"] == "base") &
              (gmi["burden_metric"] == "deaths")].copy()
    cur = cur.set_index("pathogen")
    cur_burden = cur["burden_value"].copy()
    cur_funding_share = cur["funding_share"].copy()  # held constant for projection

    print("=== Current (2019 burden, cumulative funding) ===")
    show = cur[["funding_share", "burden_share", "misalignment"]] * 100
    show.columns = ["funding %", "burden %", "f-b pp"]
    show = show.round(1)
    print(show.to_string())
    current_gmi = float(cur["misalignment"].abs().sum() / 2)
    print(f"  current GMI = {current_gmi*100:.1f}%")

    # ----- Load projection draws -----
    klpn_d = pd.read_parquet(KLPN_DRAWS)[PROJ_YEAR].values
    acba_d = pd.read_parquet(ACBA_DRAWS)[PROJ_YEAR].values
    n_draws = min(len(klpn_d), len(acba_d))
    klpn_d, acba_d = klpn_d[:n_draws], acba_d[:n_draws]
    print(f"\nUsing {n_draws:,} posterior draws to propagate uncertainty.")

    # ----- Build per-draw projected burdens -----
    # Scale Kp and Ab burdens by projected_prev / 2019_prev; others unchanged.
    pathogens = list(cur_burden.index)
    projected = np.tile(cur_burden.values.astype(float),
                        (n_draws, 1))  # shape (n_draws, n_pathogens)

    kp_idx = pathogens.index("klebsiella pneumoniae")
    ab_idx = pathogens.index("acinetobacter baumannii")
    projected[:, kp_idx] *= klpn_d / KLPN_PREV_2019
    projected[:, ab_idx] *= acba_d / ACBA_PREV_2019

    # Burden shares per draw
    totals = projected.sum(axis=1, keepdims=True)
    burden_share_proj = projected / totals          # (n_draws, n_pathogens)

    # Misalignment per pathogen per draw (funding share is fixed)
    fund_share = cur_funding_share.values
    mis_proj = fund_share[None, :] - burden_share_proj  # (n_draws, n_pathogens)

    # Projected GMI per draw = half-sum-abs
    gmi_proj = np.abs(mis_proj).sum(axis=1) / 2  # (n_draws,)

    # ----- Summaries -----
    print("\n=== Projected 2030 (XDR-scaled Kp & Ab burdens; others held; "
          "funding held) ===")
    mis_mean = mis_proj.mean(axis=0)
    mis_lo80 = np.quantile(mis_proj, 0.10, axis=0)
    mis_hi80 = np.quantile(mis_proj, 0.90, axis=0)
    bsh_mean = burden_share_proj.mean(axis=0)

    rep = pd.DataFrame({
        "pathogen": pathogens,
        "funding %":           (fund_share * 100).round(1),
        "burden % (2019)":     (cur["burden_share"].values * 100).round(1),
        "burden % (2030 mean)":(bsh_mean * 100).round(1),
        "f-b pp (2019)":       (cur["misalignment"].values * 100).round(1),
        "f-b pp (2030 mean)":  (mis_mean * 100).round(1),
        "f-b 2030 80% UI lo":  (mis_lo80 * 100).round(1),
        "f-b 2030 80% UI hi":  (mis_hi80 * 100).round(1),
    })
    print(rep.to_string(index=False))

    print(f"\n=== GMI: current vs projected 2030 ===")
    print(f"  GMI current (2019 burden):  {current_gmi*100:.1f}%")
    print(f"  GMI 2030 (mean):            {gmi_proj.mean()*100:.1f}%")
    print(f"  GMI 2030 (80% UI):          "
          f"{np.quantile(gmi_proj, 0.10)*100:.1f}%  to  "
          f"{np.quantile(gmi_proj, 0.90)*100:.1f}%")
    print(f"  GMI 2030 (95% UI):          "
          f"{np.quantile(gmi_proj, 0.025)*100:.1f}%  to  "
          f"{np.quantile(gmi_proj, 0.975)*100:.1f}%")

    # ----- Save per-draw projected misalignment table -----
    long = []
    for d in range(n_draws):
        for i, p in enumerate(pathogens):
            long.append({
                "draw": d, "pathogen": p,
                "funding_share": float(fund_share[i]),
                "burden_share_2030": float(burden_share_proj[d, i]),
                "misalignment_2030": float(mis_proj[d, i]),
            })
    out_df = pd.DataFrame(long)
    out_path = PROC / "gmi_projected_2030.parquet"
    out_df.to_parquet(out_path, index=False)
    print(f"\nwrote {out_path.relative_to(REPO_ROOT)} "
          f"({len(out_df):,} rows: {n_draws} draws x {len(pathogens)} pathogens)")

    # ----- Honest framing reminder -----
    print("\nAssumption (D-019, state in report): the four unprojected pathogens "
          "(E. coli, S. aureus, S. pneumoniae, P. aeruginosa) hold their 2019 "
          "burdens. The projected GMI reflects only Kp + Ab dynamics.")


if __name__ == "__main__":
    main()
