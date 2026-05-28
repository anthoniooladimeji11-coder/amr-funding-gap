"""
Step 1 — Data harmonization for the AMR burden–funding gap project.

Three Vivli surveillance datasets arrive in three different shapes:

  - ATLAS     : 388 MB CSV, wide. Each drug has a value column AND an
                interpretation column (e.g. 'Meropenem' + 'Meropenem_I').
                Trailing ~25 columns are beta-lactamase GENE flags, not MICs.
  - SIDERO-WT : xlsx, wide, MIC values only (no S/I/R interpretation column).
                Antibiotic headers contain stray spaces ('Ceftazidime/ Avibactam').
  - KEYSTONE  : xlsx, wide, MIC values only. Headers contain newlines ('\\n').

This module maps all three into ONE common LONG-format schema:

    dataset, isolate_id, organism, country, region, year,
    specimen, age, gender, antibiotic, mic, interpretation

Long format = one row per (isolate x antibiotic) instead of one wide row per
isolate. This is what the LCA (Step 2) and the Bayesian model (Step 3) consume,
and it makes the three datasets line up despite their different layouts.

Beta-lactamase gene flags from ATLAS are kept SEPARATELY as an isolate-level
attribute table (they are genotype, not MIC, so they do not belong in the MIC
long table).

NOTE: this module does NOT classify MICs into S/I/R. Applying CLSI/EUCAST
breakpoints is a separate, clinically-sensitive step (the next module).
Here, 'interpretation' is only populated where the source file already provides
it (ATLAS), and is left null elsewhere.

Usage (from repo root, venv active):
    python -m amr_gap.harmonize --all
    python -m amr_gap.harmonize --dataset sidero
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------
# Paths. Adjust RAW_DIR if your folder name differs.
# --------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "raw" / "Data Challenge"
INTERIM_DIR = REPO_ROOT / "data" / "interim"

ATLAS_CSV = RAW_DIR / "ATLAS_Antibiotics" / "atlas_vivli_2004_2024.csv"
SIDERO_XLSX = (
    RAW_DIR
    / "SIDERO-WT"
    / "Updated_Shionogi Five year SIDERO-WT Surveillance data(without strain number)_Vivli_220409.xlsx"
)
KEYSTONE_XLSX = (
    RAW_DIR / "KEYSTONE" / "Omadacycline_2015 to 2025_Surveillance_data.xlsx"
)

# The common schema every dataset is mapped into.
COMMON_COLS = [
    "dataset",
    "isolate_id",
    "organism",
    "country",
    "region",
    "year",
    "specimen",
    "age",
    "gender",
    "antibiotic",
    "mic",
    "mic_numeric",
    "mic_censor",
    "interpretation",
]

# Your two core modelling pathogens (used only to write a scoped subset;
# the full harmonized data keeps ALL organisms).
CORE_ORGANISMS = ["klebsiella pneumoniae", "acinetobacter baumannii"]


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def clean_header(name: str) -> str:
    """Strip newlines / collapse repeated spaces in a column header."""
    return re.sub(r"\s+", " ", str(name).replace("\n", " ")).strip()


def norm_organism(name) -> str | float:
    """Lowercase + collapse spaces so 'Klebsiella  pneumoniae' == core list."""
    if pd.isna(name):
        return name
    return re.sub(r"\s+", " ", str(name)).strip().lower()


def parse_mic(value):
    """
    Split a raw MIC cell into (numeric_value, censor_flag).

    AMR MICs are often censored at the assay's limits, e.g.:
        '<=0.004'  -> (0.004, '<=')   at or below the lowest tested dilution
        '>64'      -> (64.0,  '>')     above the highest tested dilution
        '0.12'     -> (0.12,  'exact') a precise reading
    Returns (float|nan, str). The censor flag is kept because '<=0.004' is NOT
    the same as 0.004 for log2 transforms or breakpoint classification.
    """
    if pd.isna(value):
        return (float("nan"), "missing")
    s = str(value).strip()
    m = re.match(r"^\s*(<=|>=|<|>|=)?\s*([0-9]*\.?[0-9]+)\s*$", s)
    if not m:
        return (float("nan"), "unpar," + s[:12])  # keep a trace of odd values
    sym = m.group(1) or "exact"
    if sym == "=":
        sym = "exact"
    try:
        num = float(m.group(2))
    except ValueError:
        return (float("nan"), "unpar")
    return (num, sym)


def norm_antibiotic(name: str) -> str:
    """Standardize an antibiotic column name: clean, lowercase, unify separators."""
    n = clean_header(name).lower()
    # unify 'a/ b', 'a-b', 'a b' style combination separators to 'a-b'
    n = n.replace("/", "-").replace(" - ", "-")
    n = re.sub(r"\s*-\s*", "-", n)
    n = re.sub(r"\s+", " ", n).strip()
    return n


def _melt_wide(
    df: pd.DataFrame,
    *,
    dataset: str,
    id_map: dict,
    antibiotic_cols: list[str],
    interp_map: dict | None = None,
) -> pd.DataFrame:
    """
    Turn a wide per-isolate table into the common long format.

    id_map: maps source id/metadata columns -> common schema names, e.g.
            {'Species': 'organism', 'Year': 'year', ...}
    antibiotic_cols: the source columns holding MIC values.
    interp_map: optional {antibiotic_value_col: interpretation_col} (ATLAS only).
    """
    df = df.copy()

    # Build the isolate-level metadata frame in common-schema names.
    meta = pd.DataFrame()
    for src, common in id_map.items():
        meta[common] = df[src] if src in df.columns else pd.NA
    # Ensure every metadata column exists even if absent in this dataset.
    for c in ["isolate_id", "organism", "country", "region", "year",
              "specimen", "age", "gender"]:
        if c not in meta.columns:
            meta[c] = pd.NA
    meta["dataset"] = dataset
    meta["organism"] = meta["organism"].map(norm_organism)
    # Age/gender are encoded differently across datasets (ATLAS uses age BANDS
    # like '61+', KEYSTONE uses numeric age). Keep them as consistent strings so
    # they concatenate cleanly; we are not doing arithmetic on them in the core
    # model. A numeric band-midpoint can be derived later if ever needed.
    for c in ["age", "gender", "country", "region", "specimen", "isolate_id"]:
        meta[c] = meta[c].astype("string")
    # Year: coerce to nullable integer (handles stray non-numeric years safely).
    meta["year"] = pd.to_numeric(meta["year"], errors="coerce").astype("Int64")

    # Long-format MIC rows.
    frames = []
    for col in antibiotic_cols:
        sub = meta.copy()
        sub["antibiotic"] = norm_antibiotic(col)
        raw = df[col].values
        sub["mic"] = pd.Series(raw, index=sub.index).astype("string")
        parsed = pd.Series(raw, index=sub.index).map(parse_mic)
        sub["mic_numeric"] = parsed.map(lambda t: t[0])
        sub["mic_censor"] = parsed.map(lambda t: t[1])
        if interp_map and col in interp_map and interp_map[col] in df.columns:
            sub["interpretation"] = df[interp_map[col]].values
        else:
            sub["interpretation"] = pd.NA
        frames.append(sub)

    long = pd.concat(frames, ignore_index=True)
    # Drop rows with no usable MIC measurement (numeric parse failed/missing).
    long = long[long["mic_numeric"].notna()].reset_index(drop=True)
    return long[COMMON_COLS]


# --------------------------------------------------------------------------
# Per-dataset harmonizers
# --------------------------------------------------------------------------
def harmonize_sidero() -> pd.DataFrame:
    df = pd.read_excel(SIDERO_XLSX)
    df.columns = [clean_header(c) for c in df.columns]

    meta_cols = {
        "Organism Name": "organism",
        "Country": "country",
        "Region": "region",
        "Year Collected": "year",
        "Body Location": "specimen",
    }
    non_abx = set(meta_cols) | {"Date Collected"}
    antibiotic_cols = [c for c in df.columns if c not in non_abx]

    return _melt_wide(
        df, dataset="sidero", id_map=meta_cols, antibiotic_cols=antibiotic_cols
    )


def harmonize_keystone() -> pd.DataFrame:
    df = pd.read_excel(KEYSTONE_XLSX)
    df.columns = [clean_header(c) for c in df.columns]

    meta_cols = {
        "Collection Number": "isolate_id",
        "Organism": "organism",
        "Country": "country",
        "Continent": "region",
        "Study Year": "year",
        "Specimen Type": "specimen",
        "Age": "age",
        "Gender": "gender",
    }
    # Everything that is NOT metadata/context is treated as an antibiotic MIC col.
    context_cols = set(meta_cols) | {
        "US Census Division", "Nosocomial", "Medical Service",
        "Infection Source", "Infection Type", "Source of Bloodstream infection",
        "Ventilator-Associated Pneumonia", "Intensive Care Unit (ICU)",
        "Cystic Fibrosis (CF) Patient",
    }
    antibiotic_cols = [c for c in df.columns if c not in context_cols]

    return _melt_wide(
        df, dataset="keystone", id_map=meta_cols, antibiotic_cols=antibiotic_cols
    )


def harmonize_atlas(chunksize: int = 100_000) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    ATLAS is 388 MB — read in chunks so we never hold the whole CSV in memory.
    Returns (long_mic_table, betalactamase_gene_table).
    """
    # Peek at headers once to classify columns.
    head = pd.read_csv(ATLAS_CSV, nrows=1)
    cols = list(head.columns)

    interp_cols = [c for c in cols if c.endswith("_I")]
    value_cols = [c for c in cols if (c + "_I") in interp_cols]
    interp_map = {c: c + "_I" for c in value_cols}

    meta_cols = {
        "Isolate Id": "isolate_id",
        "Species": "organism",
        "Country": "country",
        "State": "region",
        "Year": "year",
        "Source": "specimen",
        "Age Group": "age",
        "Gender": "gender",
    }

    # Beta-lactamase gene flags = trailing columns that are neither metadata,
    # nor MIC value cols, nor interpretation cols.
    known = set(meta_cols) | set(value_cols) | set(interp_cols) | {"Study", "Speciality"}
    gene_cols = [c for c in cols if c not in known]

    long_parts: list[pd.DataFrame] = []
    gene_parts: list[pd.DataFrame] = []

    for chunk in pd.read_csv(ATLAS_CSV, chunksize=chunksize, low_memory=False):
        long_parts.append(
            _melt_wide(
                chunk,
                dataset="atlas",
                id_map=meta_cols,
                antibiotic_cols=value_cols,
                interp_map=interp_map,
            )
        )
        g = chunk[["Isolate Id"] + gene_cols].copy()
        g = g.rename(columns={"Isolate Id": "isolate_id"})
        gene_parts.append(g)

    long = pd.concat(long_parts, ignore_index=True)
    genes = pd.concat(gene_parts, ignore_index=True)
    return long, genes


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------
def write_parquet(df: pd.DataFrame, name: str) -> None:
    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    out = INTERIM_DIR / name
    df.to_parquet(out, index=False)
    print(f"  wrote {out.relative_to(REPO_ROOT)}  ({len(df):,} rows)")


def run(datasets: list[str]) -> None:
    long_all: list[pd.DataFrame] = []

    if "sidero" in datasets:
        print("Harmonizing SIDERO-WT ...")
        s = harmonize_sidero()
        write_parquet(s, "sidero_long.parquet")
        long_all.append(s)

    if "keystone" in datasets:
        print("Harmonizing KEYSTONE ...")
        k = harmonize_keystone()
        write_parquet(k, "keystone_long.parquet")
        long_all.append(k)

    if "atlas" in datasets:
        print("Harmonizing ATLAS (chunked) ...")
        a_long, a_genes = harmonize_atlas()
        write_parquet(a_long, "atlas_long.parquet")
        write_parquet(a_genes, "atlas_betalactamase_genes.parquet")
        long_all.append(a_long)

    if long_all:
        combined = pd.concat(long_all, ignore_index=True)
        write_parquet(combined, "all_long.parquet")
        scoped = combined[combined["organism"].isin(CORE_ORGANISMS)]
        write_parquet(scoped, "core_pathogens_long.parquet")
        print("\nSummary by dataset:")
        print(combined.groupby("dataset").size().to_string())
        print("\nCore-pathogen rows:", f"{len(scoped):,}")


def main() -> None:
    p = argparse.ArgumentParser(description="Harmonize Vivli AMR datasets.")
    p.add_argument("--dataset", choices=["atlas", "sidero", "keystone"])
    p.add_argument("--all", action="store_true", help="harmonize all three")
    args = p.parse_args()

    if args.all:
        run(["sidero", "keystone", "atlas"])
    elif args.dataset:
        run([args.dataset])
    else:
        p.error("pass --all or --dataset NAME")


if __name__ == "__main__":
    main()
