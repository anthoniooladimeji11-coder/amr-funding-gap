"""
Figure 1 for the report: GMI Quadrant Chart.

Axes:  burden share (x) vs funding share (y), both in percent.
Bubble size: late-stage pipeline share (%).
Bubble colour: quadrant (NEGLECTED / SERVED DESPITE LOW FUNDING /
                          INVESTMENT-TRANSLATION GAP / WELL-RESOURCED).
Diagonal: perfect alignment (funding share = burden share).

Reads data/processed/pipeline_cross.parquet and writes
data/processed/figures/fig1_quadrant.png and .pdf.

Run:
    PYTHONPATH=src python -m amr_gap.figure1
"""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "data" / "processed" / "pipeline_cross.parquet"
OUT = REPO_ROOT / "data" / "processed" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Quadrant -> colour (colourblind-considerate palette).
COLOR = {
    "NEGLECTED":                  "#c0392b",  # deep red: urgent
    "SERVED DESPITE LOW FUNDING": "#f39c12",  # amber: monitor
    "INVESTMENT-TRANSLATION GAP": "#8e44ad",  # purple: investigate
    "WELL-RESOURCED":             "#2c7a3d",  # green: baseline
}

# Short labels for the bubbles (italic in caption text by convention).
SHORT = {
    "klebsiella pneumoniae":    "K. pneumoniae",
    "escherichia coli":         "E. coli",
    "acinetobacter baumannii":  "A. baumannii",
    "streptococcus pneumoniae": "S. pneumoniae",
    "staphylococcus aureus":    "S. aureus",
    "pseudomonas aeruginosa":   "P. aeruginosa",
}


def main():
    df = pd.read_parquet(SRC)
    df["fund_pct"] = df["funding_share"] * 100
    df["burd_pct"] = df["burden_share"] * 100
    df["pipe_pct"] = df["pipeline_share"] * 100

    fig, ax = plt.subplots(figsize=(7.5, 6.2), dpi=200)

    # Diagonal of perfect alignment
    lo, hi = 0, max(df["fund_pct"].max(), df["burd_pct"].max()) + 5
    ax.plot([lo, hi], [lo, hi], color="#888", linestyle="--", linewidth=1,
            zorder=1, label="Perfect alignment (funding = burden)")

    # Bubbles: area ~ pipeline share (so size visually represents share)
    # Scale: 40 area-units per 1% pipeline share, with a floor for visibility.
    for _, r in df.iterrows():
        ax.scatter(r["burd_pct"], r["fund_pct"],
                   s=80 + r["pipe_pct"] * 60,
                   c=COLOR[r["quadrant"]],
                   alpha=0.78, edgecolors="white", linewidths=1.5, zorder=3)

    # Pathogen labels: offset to avoid bubble overlap
    LABEL_OFFSETS = {
        "klebsiella pneumoniae":    (1.5, -1.5),
        "escherichia coli":         (1.5,  1.0),
        "acinetobacter baumannii":  (-5.0, -3.0),
        "streptococcus pneumoniae": (1.5,  1.0),
        "staphylococcus aureus":    (1.5,  1.0),
        "pseudomonas aeruginosa":   (1.5,  2.0),
    }
    for _, r in df.iterrows():
        dx, dy = LABEL_OFFSETS.get(r["pathogen"], (1.5, 1.0))
        ax.annotate(SHORT[r["pathogen"]],
                    xy=(r["burd_pct"], r["fund_pct"]),
                    xytext=(r["burd_pct"] + dx, r["fund_pct"] + dy),
                    fontsize=10, fontstyle="italic", zorder=4)

    # Quadrant legend (colour swatches), placed bottom-right
    handles = [plt.Line2D([], [], marker='o', linestyle='', markersize=10,
                          markerfacecolor=c, markeredgecolor='white',
                          label=k) for k, c in COLOR.items()]
    leg = ax.legend(handles=handles, title="Quadrant",
                    loc="lower right", fontsize=8.5, title_fontsize=9,
                    frameon=True, framealpha=0.95)
    leg.get_frame().set_edgecolor("#bbb")

    # Bubble-size legend (separate, top-left)
    size_legend_vals = [5, 15, 30]
    size_handles = [plt.scatter([], [], s=80 + v * 60, c="#888",
                                edgecolors="white", linewidths=1.2,
                                label=f"{v}%") for v in size_legend_vals]
    size_leg = ax.legend(handles=size_handles, title="Pipeline share",
                         loc="upper left", fontsize=8.5, title_fontsize=9,
                         frameon=True, framealpha=0.95, labelspacing=1.6,
                         borderpad=1.0)
    size_leg.get_frame().set_edgecolor("#bbb")
    ax.add_artist(leg)  # keep the first legend after adding the second

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Attributable burden share (deaths, GRAM 2019), %", fontsize=11)
    ax.set_ylabel("R&D funding share (Hub base case), %", fontsize=11)
    ax.set_title("Figure 1. AMR Funding Misalignment vs Burden, with Pipeline Cross-Flag",
                 fontsize=11.5, pad=12, loc="left")
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # Annotation for the diagonal interpretation, placed clear of bubbles
    ax.text(hi * 0.02, hi * 0.55,
            "above diagonal: overfunded\nbelow diagonal: underfunded",
            ha="left", va="top", fontsize=8.5, color="#666",
            bbox=dict(facecolor="white", edgecolor="#ddd",
                      boxstyle="round,pad=0.4"))

    plt.tight_layout()
    out_png = OUT / "fig1_quadrant.png"
    out_pdf = OUT / "fig1_quadrant.pdf"
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out_png.relative_to(REPO_ROOT)}")
    print(f"wrote {out_pdf.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
