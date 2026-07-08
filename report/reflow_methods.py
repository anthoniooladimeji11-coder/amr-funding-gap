#!/usr/bin/env python3
"""Rewrite the Methods section to flow as prose."""

from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "report" / "draft_v2.md"
BAK = REPO / "report" / "draft_v2.md.bak4"

METHODS_OLD_START = "## 2. Methods"
METHODS_END_MARKER = "## 3. Results"

METHODS_NEW = """## 2. Methods

We assembled and verified seven analytical layers; full code is at the linked repository, with a complete decision log preserved as `DECISIONS.md`.

For the surveillance layer, we unified three Vivli AMR Register datasets (ATLAS, KEYSTONE, SIDERO-WT) into a long-format table of 17,477,785 isolate-drug records spanning 2004 to 2024 across 89 countries. Censored MICs (e.g. `<=0.5`, `>16`) are preserved with explicit censoring flags, and isolate identity is guaranteed unique by construction (dataset-prefixed source-row key) to enable correct per-isolate pivoting. Each MIC was then classified S/I/R against EUCAST Clinical Breakpoint Tables v14.0, verified line-by-line against the source PDF (25 verified breakpoint rows plus seven confirmed no-breakpoint pathogen-drug pairs). The classifier is censoring-aware: a value of `<=v` returns S only if v sits at or below the EUCAST S threshold, otherwise it returns Unresolved. CLSI was dropped from the analysis as values could not be verified to primary source.

For resistance phenotyping, we fitted Latent Class Analysis models on *Klebsiella pneumoniae* and *Acinetobacter baumannii* using `stepmix` with binary encoding and missingness handled by full-information ML, over the variation-bearing drugs in each pathogen's panel. *A. baumannii* resolved into 4 phenotypes (BIC-optimal): susceptible (33%), an XDR class resistant to all panel drugs except colistin (52%), a carbapenem-resistant/aminoglycoside-susceptible class (8%), and an aminoglycoside-and-fluoroquinolone-resistant/carbapenem-susceptible class (7%). *K. pneumoniae* was fixed at K=5 by clinical and parsimony judgement (BIC favoured K=6 but the sixth class was a rare 0.6% variant): susceptible, aminoglycoside-susceptible ESBL, aminoglycoside-resistant ESBL, XDR (with elevated colistin resistance, 19%), and a small carbapenem-intermediate group.

For funding attribution, we built a Funding Allocation Profile by mapping projects to pathogens via the Hub `Categories` coded field (e.g. `1504 Infectious Agent / Bacteria / Gram negative / Klebsiella spp.`), which provides structured taxonomic attribution more reliable than the free-text agent field. Multi-pathogen projects split USD equally among named target genera, and funding is recorded at project start year. Two sensitivity variants redistribute broad-spectrum funding (projects targeting "Gram negative" or "Gram positive" without naming a target genus): a *proportional* variant allocates the broad pool in proportion to base shares within each Gram class, and an *equal* variant allocates it evenly across that class's target genera. Gram-negative and Gram-positive broad pools are kept separate, so a Gram-positive project is not a candidate for redistribution onto Gram-negative pathogens.

Per-pathogen attributable deaths and disability-adjusted life-years (DALYs) for 2019 were extracted from the Lancet GRAM appendix, Table S22, row "Resistance to one or more antibiotics" (attributable columns; not associated). The six pathogens together account for 929,000 attributable deaths globally, cross-validating the extraction against the GRAM paper's headline figure. For each pathogen we then compute funding share (FAP base case across all years) and burden share (GRAM attributable deaths 2019), with per-pathogen misalignment defined as funding share minus burden share in percentage points. The Global Misalignment Index is the half-sum of absolute per-pathogen misalignments, an index in [0, 1] interpretable as the share of total funding that would need to shift to align with burden. Three funding cases (base, proportional, equal) crossed with two burden metrics (deaths, DALYs) yield six sensitivity variants. We further count late-stage therapeutic projects (Hub `Categories` containing "Therapeutics / Development" or "Therapeutics / Approval & Post approval", consistent with WHO and PEW antibacterial pipeline conventions) per pathogen, and classify pathogens into four quadrants by crossing misalignment sign with pipeline share, split at the median across the six.

For the temporal extension, we fit a Beta-binomial random-walk model on the logit of XDR prevalence to *K. pneumoniae* and *A. baumannii* (the two pathogens with locked LCA phenotypes), using historical observations through 2024 (2025 excluded as incomplete reporting). The XDR class for each pathogen is identified programmatically by profile rather than class number. The model, theta_t = theta_{t-1} + Normal(0, sigma) with weakly informative priors, propagates uncertainty into a 2030 forecast via continued-random-walk sampling. We translate projected prevalences into projected pathogen burdens under the stated assumption that attributable burden scales linearly with XDR prevalence relative to 2019. The four unprojected pathogens hold their 2019 burdens; funding shares are held constant. Uncertainty is propagated to the projected GMI by recomputing the index from each of 8,000 posterior draws.

Funding and burden are both at global, pathogen level. The funding side cannot be drug-class resolved (no such field in the Hub data); the burden side uses pathogen totals rather than pathogen-by-drug-class. Both sides match. Funding aggregates cumulative investment; burden is the 2019 cross-section. We frame the comparison as cumulative R&D effort against the burden it is intended to address.

"""


def main():
    text = SRC.read_text()
    BAK.write_text(text)
    print(f"backed up to {BAK.relative_to(REPO)}")

    if METHODS_OLD_START not in text or METHODS_END_MARKER not in text:
        print("Methods section markers not found; check draft_v2.md manually.")
        return

    start = text.index(METHODS_OLD_START)
    end = text.index(METHODS_END_MARKER)
    text = text[:start] + METHODS_NEW + text[end:]

    SRC.write_text(text)
    print(f"wrote {SRC.relative_to(REPO)}")


if __name__ == "__main__":
    main()
