"""
Probe the Hub 'Categories' coded hierarchy for reliable pathogen attribution.

Each Categories cell is newline-separated coded lines like:
    1100 Sector / Human / Human
    3200 Research Area / Therapeutics / Development
    1504 Infectious Agent / Bacteria / Bacteria / Klebsiella pneumoniae

Goal: find the stable numeric codes for K. pneumoniae and A. baumannii, and
the codes for the 'Gram negative' broad group, so funding can be attributed by
CODE (robust) rather than free-text matching. Also report how attribution by
code compares to the free-text counts (509 Kp / 414 Ab) found earlier.

Run:
    PYTHONPATH=src python -m amr_gap.hub_categories
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
HUB = REPO_ROOT / "data" / "raw" / "Projects.xlsx"


def main() -> None:
    df = pd.read_excel(HUB, sheet_name="data")
    cats = df["Categories"].astype("string").fillna("")

    # Collect all distinct coded lines that mention 'Infectious Agent'.
    agent_lines = Counter()
    for cell in cats:
        for line in cell.split("\n"):
            line = line.strip()
            if "Infectious Agent" in line:
                agent_lines[line] += 1

    print("=== INFECTIOUS-AGENT CODED LINES mentioning our pathogens ===")
    for line, n in agent_lines.most_common():
        low = line.lower()
        if "klebsiella" in low or "acinetobacter" in low or "gram negative" in low:
            print(f"  ({n:>5}x) {line}")
    print()

    # Extract the leading numeric code for any line naming each pathogen.
    def codes_for(term):
        cs = Counter()
        for line in agent_lines:
            if term in line.lower():
                m = re.match(r"\s*(\d+)\s", line)
                if m:
                    cs[m.group(1)] += 1
        return cs

    print("=== CANDIDATE CODES ===")
    for term in ["klebsiella pneumoniae", "klebsiella",
                 "acinetobacter baumannii", "acinetobacter", "gram negative"]:
        print(f"  '{term}': leading codes -> {dict(codes_for(term))}")
    print()

    # Attribution by code: how many projects + $ have a Categories cell whose
    # agent lines reference each pathogen (by name within Categories).
    amt = pd.to_numeric(df["Amount USD"], errors="coerce").fillna(0)
    for term in ["klebsiella", "acinetobacter", "gram negative"]:
        hit = cats.str.lower().str.contains(
            r"infectious agent.*" + term, na=False, regex=True)
        # the regex above won't span newlines well; do a line-aware check too
        hit2 = cats.apply(
            lambda c: any(("infectious agent" in ln.lower() and term in ln.lower())
                          for ln in c.split("\n")))
        print(f"  Categories names '{term}': "
              f"{int(hit2.sum()):,} projects, ${amt[hit2].sum():,.0f}")
    print()
    print("Compare to free-text earlier: klebsiella 509/$637M, "
          "acinetobacter 414/$554M.")
    print("If code/line attribution is cleaner & consistent, use it for the FAP.")


if __name__ == "__main__":
    main()
