#!/usr/bin/env python3
"""Patch draft_v2.md for PDF conversion: real figure embeds + expanded refs.

Backs up the original to draft_v2.md.bak before editing. Three edits:
  1. Insert a markdown image for Figure 1 BEFORE the existing text placeholder
     on the FIGURE 1 line (the placeholder stays as the figure caption text).
  2. Insert a markdown image and a brief Figure 2 caption immediately after
     the Section 3b text (just before Section 4 starts).
  3. Replace the compressed References block at the end with the expanded
     references that include full DOIs.

Run from the repository root:
    python report/patch_draft_for_pdf.py
"""

from __future__ import annotations
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "report" / "draft_v2.md"
BAK = REPO / "report" / "draft_v2.md.bak"

FIG1_IMG = (
    "![Figure 1. AMR funding misalignment vs burden, with pipeline cross-flag. "
    "Bubble size encodes late-stage pipeline share; colour encodes quadrant.]"
    "(../data/processed/figures/fig1_quadrant.png){ width=85% }\n\n"
)

FIG2_BLOCK = (
    "\n![Figure 2. XDR prevalence with Bayesian random-walk projection to 2030. "
    "Observed dots sized by isolate count; blue band = posterior 80% UI over the "
    "fitted period; red fan = projection 80% UI to 2030.]"
    "(../data/processed/figures/fig2_projection.png){ width=98% }\n\n"
)

EXPANDED_REFS = """## References

1. Antimicrobial Resistance Collaborators (Murray CJL, Ikuta KS, Sharara F, et al.). Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *Lancet* 2022; 399(10325): 629–655. DOI: [10.1016/S0140-6736(21)02724-0](https://doi.org/10.1016/S0140-6736(21)02724-0).

2. European Committee on Antimicrobial Susceptibility Testing (EUCAST). Breakpoint tables for interpretation of MICs and zone diameters, Version 14.0. Växjö, Sweden: EUCAST; 2024. Available at: [https://www.eucast.org/clinical_breakpoints](https://www.eucast.org/clinical_breakpoints).

3. World Health Organization. WHO Bacterial Priority Pathogens List, 2024: Bacterial pathogens of public health importance to guide research, development and strategies to prevent and control antimicrobial resistance. Geneva: WHO; 2024. Available at: [https://www.who.int/publications/i/item/9789240093461](https://www.who.int/publications/i/item/9789240093461).

4. World Health Organization. 2023 Antibacterial Agents in Clinical and Preclinical Development: An overview and analysis. Geneva: WHO; 2024. Available at: [https://www.who.int/publications/i/item/9789240094468](https://www.who.int/publications/i/item/9789240094468).

5. Magiorakos AP, Srinivasan A, Carey RB, et al. Multidrug-resistant, extensively drug-resistant and pandrug-resistant bacteria: an international expert proposal for interim standard definitions for acquired resistance. *Clin Microbiol Infect* 2012; 18(3): 268–281. DOI: [10.1111/j.1469-0691.2011.03570.x](https://doi.org/10.1111/j.1469-0691.2011.03570.x).

6. Stan Development Team. *Stan Modeling Language Users Guide and Reference Manual*, Version 2.39. 2024. Available at: [https://mc-stan.org](https://mc-stan.org).

7. Morimoto T, Vialaneix N. StepMix: A Python package for pseudo-likelihood estimation of generalized mixture models with external variables. *J Open Source Softw* 2023; 8(85): 5012. DOI: [10.21105/joss.05012](https://doi.org/10.21105/joss.05012).

8. Vivli AMR Register. The Vivli AMR Register: An open data platform for antimicrobial resistance surveillance. Available at: [https://amr.vivli.org](https://amr.vivli.org).

9. Global AMR R&D Hub. Dynamic Dashboard: Mapping the global AMR R&D landscape. Berlin: Global AMR R&D Hub; 2024. Available at: [https://globalamrhub.org](https://globalamrhub.org).

10. World Health Organization. Antibacterial agents in preclinical and clinical development reports. Geneva: WHO; published annually. Available at: [https://www.who.int/teams/control-of-neglected-tropical-diseases](https://www.who.int/teams/control-of-neglected-tropical-diseases).

11. O'Neill J. Tackling Drug-Resistant Infections Globally: Final Report and Recommendations. London: Review on Antimicrobial Resistance; 2016. Available at: [https://amr-review.org](https://amr-review.org).
"""


def main() -> None:
    if not SRC.exists():
        print(f"missing {SRC}; aborting", file=sys.stderr)
        sys.exit(1)
    text = SRC.read_text()
    BAK.write_text(text)
    print(f"backed up to {BAK.relative_to(REPO)}")

    # ---- Edit 1: insert Figure 1 image just before the existing placeholder line.
    marker_fig1 = "**[FIGURE 1."
    if marker_fig1 not in text:
        print("Figure 1 placeholder not found; not modifying figure 1.")
    else:
        text = text.replace(marker_fig1, FIG1_IMG + marker_fig1, 1)
        print("inserted Figure 1 image embed.")

    # ---- Edit 2: insert Figure 2 immediately before "## 4. Discussion".
    section4_marker = "## 4. Discussion and policy implications"
    if section4_marker not in text:
        print("Section 4 header not found; not inserting Figure 2.")
    else:
        text = text.replace(section4_marker, FIG2_BLOCK + section4_marker, 1)
        print("inserted Figure 2 block before Section 4.")

    # ---- Edit 3: replace the compressed References block.
    ref_marker = "## References (compressed)"
    if ref_marker not in text:
        # Maybe it's already been replaced; fall back to looking for "## References"
        if "## References" in text and "10.1016/S0140-6736(21)02724-0" in text:
            print("References look already expanded; skipping.")
        else:
            print(f"References marker '{ref_marker}' not found; review manually.",
                  file=sys.stderr)
            sys.exit(2)
    else:
        # Replace from "## References (compressed)" to the end of the references
        # (just before the closing italic paragraph or end of file).
        idx = text.index(ref_marker)
        # Keep everything before the marker; append expanded refs, then any
        # trailing footer (the "---" italics line if present).
        tail = ""
        # If there's a trailing "---" section after the references, preserve it.
        after = text[idx:]
        if "\n---\n" in after:
            footer_start = after.index("\n---\n")
            tail = after[footer_start:]
        text = text[:idx] + EXPANDED_REFS + tail
        print("replaced References block with expanded version.")

    SRC.write_text(text)
    print(f"wrote updated {SRC.relative_to(REPO)}")
    print("\nPreview the result with:  head -200 report/draft_v2.md")


if __name__ == "__main__":
    main()
