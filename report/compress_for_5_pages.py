#!/usr/bin/env python3
"""Compress the report body to 5 pages."""

from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "report" / "draft_v2.md"
BAK = REPO / "report" / "draft_v2.md.bak3"

PLACEHOLDER_OPENING = "**[FIGURE 1."
PLACEHOLDER_CLOSING = "]**\n"

REF_TYPO_OLD = "69789240094468"
REF_TYPO_NEW = "9789240094468"

LIMITATIONS_OLD = """**Limitations and caveats stated openly.** The funding side is at the genus level (the Hub does not encode species or drug class); we use the genus as a proxy for the priority species. For five of six pathogens this is a close fit; for *Streptococcus* it is the loosest (the genus includes group A and B streptococci that are not the burden-driving pneumoniae). The funding side is cumulative across years (1975–2026 in raw data, dominated by post-2013 activity); burden is the 2019 cross-section. Project counts proxy candidate counts in the pipeline measure: one drug can be funded by multiple projects, biasing all pathogens similarly and preserving relative comparability. The Vivli surveillance datasets carry no patient outcomes, so we rely on external GRAM attributable estimates rather than computing a clinical burden score directly from the Register; the SPIDAAR linkage dataset was requested but unavailable at submission time. The projection assumes attributable burden scales linearly with XDR prevalence and only projects two of the six pathogens; sensitivity to alternative scaling assumptions and projection scope is not explored at submission resolution. We frame these honestly: the misalignment finding is robust to plausible variation in attribution rules and metric choice, but the precision of the percentage-point gaps depends on the resolution at which the underlying data systems track pathogens."""

LIMITATIONS_NEW = """**Limitations and caveats.** Funding is attributed at genus level (the Hub does not encode species or drug class); we use genus as a proxy for the priority species. For five of six pathogens this is a close fit; for *Streptococcus* it is the loosest, since the genus includes non-pneumoniae species. Funding is cumulative (dominated by post-2013 activity); burden is the 2019 cross-section. Project counts proxy candidate counts in the pipeline measure, biasing all pathogens similarly. The Vivli surveillance datasets carry no patient outcomes, so we rely on external GRAM attributable estimates; the SPIDAAR linkage dataset was requested but unavailable. The projection assumes attributable burden scales linearly with XDR prevalence and projects only two of six pathogens. The misalignment finding is robust to plausible variation in attribution rules and metric choice, but percentage-point precision depends on the resolution of underlying data systems."""


def main():
    text = SRC.read_text()
    BAK.write_text(text)
    print(f"backed up to {BAK.relative_to(REPO)}")

    if PLACEHOLDER_OPENING in text:
        start = text.index(PLACEHOLDER_OPENING)
        end_marker = text.find(PLACEHOLDER_CLOSING, start)
        if end_marker != -1:
            end = end_marker + len(PLACEHOLDER_CLOSING)
            text = text[:start] + text[end:]
            text = text.replace("\n\n\n\n", "\n\n").replace("\n\n\n", "\n\n")
            print("removed Figure 1 placeholder")
        else:
            print("placeholder opening found but closing not; skipped")
    else:
        print("Figure 1 placeholder not found")

    if LIMITATIONS_OLD in text:
        text = text.replace(LIMITATIONS_OLD, LIMITATIONS_NEW, 1)
        print("tightened Limitations paragraph")
    else:
        print("Limitations paragraph not matched exactly; unchanged")

    if REF_TYPO_OLD in text:
        text = text.replace(REF_TYPO_OLD, REF_TYPO_NEW)
        print("fixed reference URL typo")
    else:
        print("Reference typo not found")

    SRC.write_text(text)
    print(f"wrote {SRC.relative_to(REPO)}")


if __name__ == "__main__":
    main()
