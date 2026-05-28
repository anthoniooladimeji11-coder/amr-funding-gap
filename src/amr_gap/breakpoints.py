"""
Step 1 (continued) — Breakpoint classification.

Turns harmonized numeric MICs into S / I / R calls using BOTH CLSI and EUCAST
breakpoints, then compares them. Designed to handle censored MICs correctly.

Breakpoint values live in data/external/breakpoints.csv and MUST be filled in
with verified values from the official tables (rows are shipped as placeholders
flagged status=VERIFY). This module will REFUSE to classify against any row
still marked VERIFY, so you cannot accidentally produce results from placeholder
numbers.

Interpretation rule (MIC in mg/L), per (organism, antibiotic, standard):
    S  if  mic <= S_max
    R  if  mic >= R_min
    I  if  S_max < mic < R_min

Censoring (mic_censor from the harmonizer) is respected:
    '<=v' (left-censored): true MIC is at or below v.
        -> if v <= S_max            => S   (certainly susceptible)
        -> else                     => UNRESOLVED (could be S or higher)
    '>v'  (right-censored): true MIC is above v.
        -> if v >= R_min            => R   (certainly resistant)
        -> else                     => UNRESOLVED
    'exact': normal rule.
This avoids the classic error of treating '<=0.5' as exactly 0.5.

Usage:
    python -m amr_gap.breakpoints            # classify core_pathogens_long
    python -m amr_gap.breakpoints --input data/interim/all_long.parquet
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
INTERIM_DIR = REPO_ROOT / "data" / "interim"
EXTERNAL_DIR = REPO_ROOT / "data" / "external"
BREAKPOINTS_CSV = EXTERNAL_DIR / "breakpoints.csv"

DEFAULT_INPUT = INTERIM_DIR / "core_pathogens_long.parquet"


def load_breakpoints() -> pd.DataFrame:
    """Load the breakpoint table, skipping comment lines, and refuse placeholders."""
    bp = pd.read_csv(BREAKPOINTS_CSV, comment="#")
    bp.columns = [c.strip() for c in bp.columns]
    # Guard: nothing still flagged VERIFY may be used.
    unverified = bp[bp["status"].str.upper() == "VERIFY"]
    verified = bp[bp["status"].str.upper() == "OK"].copy()
    if len(verified) == 0:
        raise SystemExit(
            "\nAll breakpoint rows are still placeholders (status=VERIFY).\n"
            "Fill in real values from the official EUCAST / CLSI tables in\n"
            f"  {BREAKPOINTS_CSV.relative_to(REPO_ROOT)}\n"
            "and set status=OK for verified rows before classifying.\n"
        )
    if len(unverified):
        print(f"[warn] {len(unverified)} breakpoint rows still unverified — "
              "those bug/drug/standard combos will be left unclassified.")
    for col in ["S_max", "R_min"]:
        verified[col] = pd.to_numeric(verified[col], errors="coerce")
    verified["organism"] = verified["organism"].str.strip().str.lower()
    verified["antibiotic"] = verified["antibiotic"].str.strip().str.lower()
    verified["standard"] = verified["standard"].str.strip().str.upper()
    return verified


def classify_one(mic: float, censor: str, s_max: float, r_min: float) -> str:
    """Return 'S' | 'I' | 'R' | 'U' (unresolved) for a single measurement."""
    if pd.isna(mic) or pd.isna(s_max) or pd.isna(r_min):
        return "U"
    if censor == "<=" or censor == "<":
        return "S" if mic <= s_max else "U"
    if censor == ">" or censor == ">=":
        return "R" if mic >= r_min else "U"
    # exact
    if mic <= s_max:
        return "S"
    if mic >= r_min:
        return "R"
    return "I"


def classify(df: pd.DataFrame, bp: pd.DataFrame) -> pd.DataFrame:
    """Add sir_eucast and sir_clsi columns plus an agreement flag."""
    out = df.copy()
    out["organism"] = out["organism"].astype("string").str.lower()
    out["antibiotic"] = out["antibiotic"].astype("string").str.lower()

    for standard in ["EUCAST", "CLSI"]:
        sub = bp[bp["standard"] == standard][
            ["organism", "antibiotic", "S_max", "R_min"]
        ]
        merged = out.merge(sub, on=["organism", "antibiotic"], how="left")
        col = f"sir_{standard.lower()}"
        out[col] = [
            classify_one(m, c, s, r)
            for m, c, s, r in zip(
                merged["mic_numeric"], merged["mic_censor"],
                merged["S_max"], merged["R_min"],
            )
        ]

    # Agreement only meaningful where both gave a definitive S/I/R.
    def agree(row):
        e, c = row["sir_eucast"], row["sir_clsi"]
        if e in ("U",) or c in ("U",):
            return "no_breakpoint_or_unresolved"
        return "agree" if e == c else "DISAGREE"

    out["sir_agreement"] = out.apply(agree, axis=1)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Apply CLSI/EUCAST breakpoints.")
    p.add_argument("--input", default=str(DEFAULT_INPUT))
    p.add_argument("--output", default=str(INTERIM_DIR / "core_pathogens_sir.parquet"))
    args = p.parse_args()

    bp = load_breakpoints()
    df = pd.read_parquet(args.input)
    print(f"Loaded {len(df):,} rows from {Path(args.input).name}")

    res = classify(df, bp)
    res.to_parquet(args.output, index=False)
    print(f"Wrote {args.output}")

    print("\nEUCAST calls:")
    print(res["sir_eucast"].value_counts().to_string())
    print("\nCLSI calls:")
    print(res["sir_clsi"].value_counts().to_string())
    print("\nCLSI-vs-EUCAST agreement (where both classified):")
    print(res["sir_agreement"].value_counts().to_string())


if __name__ == "__main__":
    main()
