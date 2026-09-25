"""
Figure 2 for Section 3b: XDR prevalence with Bayesian projection to 2030.

Two-panel time series:
  - Panel A: Acinetobacter baumannii XDR (2004-2030)
  - Panel B: Klebsiella pneumoniae XDR (2012-2030)

Each panel shows:
  - observed XDR % as dots (size proportional to n_isolates)
  - posterior smoothed fit as a line with 80% UI band (historical)
  - projection 2025-2030 as line continuing into a widening fan
  - vertical dashed line at 2024 separating fit from projection

Run:
    PYTHONPATH=src python -m amr_gap.figure2
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[2]
PROC = REPO_ROOT / "data" / "processed"
OUT = PROC / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PATHOGENS = [
    ("acba", "Acinetobacter baumannii", "A. baumannii", 2004),
    ("klpn", "Klebsiella pneumoniae",    "K. pneumoniae", 2012),
]

CUTOFF_YEAR = 2024  # boundary between fit and projection


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=200, sharey=False)

    for ax, (short, full, label, start_year) in zip(axes, PATHOGENS):
        # Load projection summary (per-year mean and UIs, both historical and projected)
        proj = pd.read_parquet(PROC / f"xdr_projection_{short}.parquet")
        # Load the year-only series for sample sizes
        ts = pd.read_parquet(PROC / f"xdr_timeseries_{short}_yearonly.parquet")
        ts["year"] = ts["year"].astype(int)

        # Convert posterior shares to percentages
        for c in ("p_mean", "p_lo80", "p_hi80", "p_lo95", "p_hi95"):
            proj[c] = proj[c] * 100

        hist = proj[proj["kind"] == "historical"].sort_values("year")
        pred = proj[proj["kind"] == "projected"].sort_values("year")
        # Bridge: prepend the last historical row to projection so the line connects
        bridge = pd.concat([hist.tail(1), pred], ignore_index=True)

        # ----- Posterior median line + 80% band (historical) -----
        ax.fill_between(hist["year"], hist["p_lo80"], hist["p_hi80"],
                        color="#3b6db5", alpha=0.18, linewidth=0,
                        label="Posterior 80% UI (fit)")
        ax.plot(hist["year"], hist["p_mean"],
                color="#1f3c70", linewidth=1.8, label="Posterior mean")

        # ----- Projection: line + widening 80% fan -----
        ax.fill_between(bridge["year"], bridge["p_lo80"], bridge["p_hi80"],
                        color="#c0392b", alpha=0.16, linewidth=0,
                        label="Projection 80% UI")
        ax.plot(bridge["year"], bridge["p_mean"],
                color="#c0392b", linewidth=1.8, linestyle="-",
                label="Projection mean")

        # ----- Observed dots (size proportional to n_isolates) -----
        obs_in_range = ts[ts["year"] <= CUTOFF_YEAR].copy()
        obs_in_range["pct_xdr"] = obs_in_range["n_xdr"] / obs_in_range["n_total"] * 100
        # Size: square-root scaling, between 18 and 110
        n = obs_in_range["n_total"].values
        sizes = 18 + 92 * (np.sqrt(n) - np.sqrt(n.min())) / (np.sqrt(n.max()) - np.sqrt(n.min()) + 1e-9)
        ax.scatter(obs_in_range["year"], obs_in_range["pct_xdr"],
                   s=sizes, c="#222", alpha=0.78, edgecolors="white",
                   linewidths=0.8, zorder=4, label="Observed (size = n)")

        # ----- Vertical separator -----
        ax.axvline(CUTOFF_YEAR + 0.5, linestyle="--", color="#888",
                   linewidth=1.0, alpha=0.6, zorder=2)
        ax.text(CUTOFF_YEAR + 0.7, ax.get_ylim()[1] * 0.96
                if ax.get_ylim()[1] > 0 else 5,
                "projection\n→",
                fontsize=8.5, color="#666", va="top", ha="left")

        # ----- Labelling -----
        ax.set_title(f"{label}", fontsize=12, loc="left", style="italic")
        ax.set_xlabel("Year", fontsize=10)
        ax.set_ylabel("XDR prevalence (%)", fontsize=10)
        ax.set_xlim(start_year - 0.5, 2030.5)
        ax.set_ylim(0, max(90, ax.get_ylim()[1]))
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(True, alpha=0.25, linestyle=":")
        ax.set_axisbelow(True)

    # Shared legend at the bottom (one figure)
    handles, labels = axes[0].get_legend_handles_labels()
    # Deduplicate while preserving order
    seen = set(); h2, l2 = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l); h2.append(h); l2.append(l)
    fig.legend(h2, l2, loc="lower center", ncol=5, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("Figure 2. XDR prevalence with Bayesian random-walk projection to 2030",
                 fontsize=11.5, x=0.06, ha="left", y=0.99)

    plt.tight_layout(rect=(0, 0.04, 1, 0.96))
    out_png = OUT / "fig2_projection.png"
    out_pdf = OUT / "fig2_projection.pdf"
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out_png.relative_to(REPO_ROOT)}")
    print(f"wrote {out_pdf.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
