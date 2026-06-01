"""
Step 3a - Data prep for the Bayesian XDR projection.

Identifies the XDR class for each pathogen BY PROFILE (not by class number,
which is arbitrary). For each pathogen, builds a time-series table:
    (year, country, n_total_isolates, n_xdr) suitable for a Beta-binomial model.

XDR identification (programmatic):
  For each LCA class, compute P(non-susceptible | class) across panel drugs.
  XDR class = the class with the HIGHEST mean NS rate across the panel, AND
  that mean exceeds a permissive threshold (>= 0.7) to refuse to label a non-XDR
  pathogen XDR by accident.

Writes:
    data/processed/xdr_timeseries_acba.parquet
    data/processed/xdr_timeseries_klpn.parquet
    data/processed/xdr_class_profiles.csv   (for verification)

Run:
    PYTHONPATH=src python -m amr_gap.xdr_prep
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SIR = REPO_ROOT / "data" / "interim" / "core_pathogens_sir.parquet"
PROC = REPO_ROOT / "data" / "processed"

ORG_TO_SHORT = {
    "acinetobacter baumannii": "acba",
    "klebsiella pneumoniae":   "klpn",
}

PANELS = {
    "acinetobacter baumannii": [
        "amikacin", "gentamicin", "meropenem", "imipenem", "levofloxacin",
        "colistin",
    ],
    "klebsiella pneumoniae": [
        "meropenem", "imipenem", "doripenem", "ertapenem", "amikacin",
        "gentamicin", "ciprofloxacin", "levofloxacin", "aztreonam", "colistin",
    ],
}


def identify_xdr_class(lca_df: pd.DataFrame, panel: list[str]) -> tuple[int, pd.Series]:
    """Return the LCA class index whose profile is the most resistant overall.

    Computes mean NS rate across panel drugs per class; picks the class with the
    highest mean. Refuses to label any class XDR if the max mean is below 0.7.
    """
    # The LCA parquet stores binary NS per panel drug + lca_class + lca_maxprob.
    profile = lca_df.groupby("lca_class")[panel].mean()
    mean_ns = profile.mean(axis=1)
    best = mean_ns.idxmax()
    if mean_ns[best] < 0.7:
        raise RuntimeError(
            f"No class qualifies as XDR (max mean NS = {mean_ns[best]:.2f} < 0.70). "
            f"Profile:\n{profile.round(2)}")
    return int(best), profile.loc[best]


def main():
    sir = pd.read_parquet(SIR)
    profiles_out = []

    for organism, short in ORG_TO_SHORT.items():
        panel = PANELS[organism]
        lca_path = PROC / f"lca_{short}.parquet"
        lca = pd.read_parquet(lca_path)
        xdr_k, prof = identify_xdr_class(lca, panel)
        print(f"\n{organism}:")
        print(f"  panel: {panel}")
        print(f"  classes: {sorted(lca['lca_class'].unique())}")
        print(f"  -> identified XDR class = {xdr_k}")
        print(f"  XDR profile (P(NS|class), %):")
        print((prof * 100).round(1).to_string())

        # join the XDR label back to SIR data: need to attach the LCA class
        # for each isolate, then compute per-(year, country) totals and XDR counts.
        lca_min = lca[["isolate_id", "lca_class"]].copy()

        # Build per-isolate dataframe: one row per isolate with year, country, lca_class.
        # Use any pathogen-panel row to pick year/country/dataset for the isolate.
        iso_meta = (sir[sir["organism"] == organism]
                    .drop_duplicates(subset=["isolate_id"])
                    [["isolate_id", "dataset", "country", "year"]])
        iso = iso_meta.merge(lca_min, on="isolate_id", how="inner")
        iso["is_xdr"] = (iso["lca_class"] == xdr_k).astype(int)

        # Aggregate per (year, country)
        ts = (iso.groupby(["year", "country"], dropna=True)
              .agg(n_total=("isolate_id", "count"),
                   n_xdr=("is_xdr", "sum"))
              .reset_index())
        # Also a simpler year-only aggregate (the one the projection actually uses
        # first, before regional extensions)
        ts_year = (iso.groupby("year", dropna=True)
                   .agg(n_total=("isolate_id", "count"),
                        n_xdr=("is_xdr", "sum"))
                   .reset_index())

        # Write outputs
        ts_country_path = PROC / f"xdr_timeseries_{short}.parquet"
        ts.to_parquet(ts_country_path, index=False)
        ts_year_path = PROC / f"xdr_timeseries_{short}_yearonly.parquet"
        ts_year.to_parquet(ts_year_path, index=False)

        print(f"\n  Year-only XDR prevalence (the series the model will fit):")
        show = ts_year.copy()
        show["pct_xdr"] = (show["n_xdr"] / show["n_total"] * 100).round(1)
        print(show.to_string(index=False))
        print(f"\n  wrote {ts_country_path.relative_to(REPO_ROOT)} "
              f"({len(ts)} year-country rows)")
        print(f"  wrote {ts_year_path.relative_to(REPO_ROOT)} "
              f"({len(ts_year)} year rows)")

        profiles_out.append({
            "organism": organism,
            "xdr_class": xdr_k,
            "mean_ns_profile": float((prof.mean())),
            **{f"ns_{d}": float(prof[d]) for d in panel},
        })

    pd.DataFrame(profiles_out).to_csv(PROC / "xdr_class_profiles.csv", index=False)
    print(f"\nwrote {PROC / 'xdr_class_profiles.csv'} (verification)")


if __name__ == "__main__":
    main()
