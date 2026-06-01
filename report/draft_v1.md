# A Global Misalignment Index for AMR R&D Funding: Quantifying Where Investment Falls Short of Burden, and Where the Pipeline Compounds the Gap

**Team:** Anthonio Oladimeji, Babatunde Olowu
**Submission category:** AMR Visionary Award; Global AMR R&D Hub Cross-Domain Award
**Repository:** github.com/anthoniooladimeji11-coder/amr-funding-gap

---

## Abstract (≈150 words)

Global antimicrobial resistance (AMR) research and development draws an estimated USD 18.9 billion across 18,853 funded projects, yet whether this investment aligns with the disease burden it is intended to address has not been quantified at the level of individual pathogens. We construct a Global Misalignment Index (GMI) for the six pathogens responsible for ~73% of attributable AMR deaths worldwide, integrating the Vivli AMR Register (17.5 million harmonised isolate-drug records from three surveillance programmes), Global AMR R&D Hub project data, and the GRAM 2019 Lancet attributable burden estimates. We find a misalignment of 25.3% between funding share and burden share (deaths), robust across six sensitivity variants. Cross-referencing late-stage therapeutic pipeline activity reveals that *Acinetobacter baumannii* (4 late-stage projects globally) and *Klebsiella pneumoniae* are doubly neglected: underfunded **and** under-pipelined. The framework supports evidence-based reallocation toward the WHO Priority 1 carbapenem-resistant Gram-negative pathogens.

---

## 1. Background and aim

AMR is projected to become a leading global cause of death by 2050 (Murray et al., 2022). Mobilising R&D against AMR has been a sustained policy priority for a decade, with the Global AMR R&D Hub now tracking nearly USD 19 billion in committed investment across 89 funders. A central but under-examined question is whether this investment is **proportionate to where the burden actually concentrates** at the level of individual pathogens. Existing analyses describe the pipeline (WHO Antibacterial Pipeline reports) or the burden (GRAM 2019) in isolation, but do not formally measure misalignment between the two.

This work asks: **for the leading AMR pathogens, how closely does R&D funding share track attributable burden share, and where is the gap compounded by a thin late-stage pipeline?** We answer with a Global Misalignment Index (GMI) and a pipeline cross-flag that distinguishes neglected pathogens from those served despite low funding.

## 2. Methods

We assembled and verified six analytical layers (full code at the linked repository; complete decision log preserved as `DECISIONS.md`).

**Surveillance harmonisation.** We unified three Vivli AMR Register datasets (ATLAS, KEYSTONE, SIDERO-WT) into a long-format table of 17,477,785 isolate-drug records spanning 2004 to 2024 across 89 countries. Censored MICs (e.g. `<=0.5`, `>16`) are preserved with explicit censoring flags. Isolate identity is guaranteed unique by construction (dataset-prefixed source-row key), enabling correct per-isolate pivoting.

**Susceptibility classification.** Each MIC was classified S/I/R against EUCAST Clinical Breakpoint Tables v14.0 (verified line-by-line against the source PDF; 25 verified breakpoint rows plus seven confirmed no-breakpoint pathogen-drug pairs). The classifier is censoring-aware: a value of `<=v` returns S only if v sits at or below the EUCAST S threshold, otherwise it returns Unresolved. CLSI was dropped from the analysis as values could not be verified to primary source.

**Resistance phenotyping (Latent Class Analysis).** For *Klebsiella pneumoniae* and *Acinetobacter baumannii* we fitted Latent Class Analysis models (`stepmix`, binary encoding, missingness handled by full-information ML) over the variation-bearing drugs in each pathogen's panel. *A. baumannii* resolved into 4 phenotypes (BIC-optimal): susceptible (33%), an XDR class resistant to all panel drugs except colistin (52%), a carbapenem-resistant/aminoglycoside-susceptible class (8%), and an aminoglycoside-and-fluoroquinolone-resistant/carbapenem-susceptible class (7%). *K. pneumoniae* was fixed at K=5 by clinical and parsimony judgement (BIC favoured K=6 but the sixth class was a rare 0.6% variant): susceptible, aminoglycoside-susceptible ESBL, aminoglycoside-resistant ESBL, XDR (with elevated colistin resistance, 19%), and a small carbapenem-intermediate group.

**Funding Allocation Profile (FAP).** Project funding was attributed to pathogens via the Hub `Categories` coded field (e.g. `1504 Infectious Agent / Bacteria / Gram negative / Klebsiella spp.`), which provides structured taxonomic attribution and is more reliable than the free-text agent field. Multi-pathogen projects split USD equally among named target genera. Funding is recorded at project start year. Two sensitivity variants redistribute broad-spectrum funding (projects targeting "Gram negative" or "Gram positive" without naming a target genus): a *proportional* variant allocates the broad pool in proportion to base shares within each Gram class, and an *equal* variant allocates the broad pool evenly across that class's target genera. Gram-negative and Gram-positive broad pools are kept separate: a Gram-positive project is not a candidate for redistribution onto Gram-negative pathogens.

**Burden weights.** Per-pathogen attributable deaths and disability-adjusted life-years (DALYs) for 2019 were extracted from the Lancet GRAM appendix, Table S22, row "Resistance to one or more antibiotics" (attributable columns; not associated). The six pathogens together account for 929,000 attributable deaths globally, cross-validating the extraction against the GRAM paper's headline figure.

**The Misalignment Index.** For each pathogen, we compute funding share (FAP base case across all years) and burden share (GRAM attributable deaths 2019). Per-pathogen misalignment is funding share minus burden share (in percentage points). The Global Misalignment Index is the half-sum of absolute per-pathogen misalignments, an index in [0, 1] interpretable as the share of total funding that would need to shift to align with burden. Three funding cases (base, proportional, equal) crossed with two burden metrics (deaths, DALYs) yield six sensitivity variants.

**Pipeline cross-flag.** We count late-stage therapeutic projects (Hub `Categories` containing "Therapeutics / Development" or "Therapeutics / Approval & Post approval", consistent with WHO and PEW antibacterial pipeline conventions) per pathogen. Pathogens are classified into four quadrants by crossing misalignment sign with pipeline share (split at the median across the six).

**Resolution alignment.** Funding and burden are both at global, pathogen level. The funding side cannot be drug-class resolved (no such field in the Hub data); the burden side uses pathogen totals rather than pathogen-by-drug-class. Both sides match. Funding aggregates cumulative investment; burden is the 2019 cross-section. We frame the comparison as cumulative R&D effort against the burden it is intended to address.

## 3. Results

The Global Misalignment Index is **25.3%** for the six leading pathogens, base case (deaths). This implies that approximately a quarter of pathogen-attributed funding (USD 4.17 billion across the six in our base FAP) would need to redistribute to align with burden share. The index is stable across alternative metrics and funding assumptions (Table 1), ranging from 15.7% to 25.3% across all six variants; deaths versus DALYs barely moves the index.

**Table 1. GMI under sensitivity variants (% misalignment).**

| Funding case | Deaths | DALYs |
|---|---|---|
| Base (pathogen-named only) | 25.3 | 23.0 |
| Proportional broad-pool redistribution | 19.1 | 20.2 |
| Equal broad-pool redistribution | 15.7 | 16.7 |

The per-pathogen pattern (Figure 1) is clinically coherent. The three underfunded pathogens are Gram-negatives: *K. pneumoniae* (8% of funding against 21% of burden; gap of –12.4 percentage points), *E. coli* (–6.8 pp), and *A. baumannii* (–6.1 pp). The three overfunded pathogens are *P. aeruginosa* (+11.7 pp), *S. aureus* (+7.6 pp), and *S. pneumoniae* (+5.9 pp).

**The pipeline cross-flag sharpens this into an actionable map (Figure 1, Table 2).** *A. baumannii* and *K. pneumoniae* both fall into the **NEGLECTED** quadrant: underfunded relative to burden **and** under-represented in the late-stage pipeline. *A. baumannii* is the starkest case in the dataset, with only 4 late-stage therapeutic projects identified worldwide (5.3% of the 251 late-stage projects across our six pathogens) against 14% of attributable AMR deaths. *E. coli* falls into **SERVED DESPITE LOW FUNDING**: its funding share lags its burden, but its late-stage pipeline (25%) is healthy, consistent with substantial broad-spectrum and urinary-tract-infection-targeted development not always captured in pathogen-named funding. *S. pneumoniae* sits in **INVESTMENT-TRANSLATION GAP** (overfunded with thin pipeline at 2.6%), an artefact partly of vaccine-heavy research that does not produce novel antibiotic candidates, and partly of the genus-to-species attribution caveat noted below. *S. aureus* and *P. aeruginosa* are **WELL-RESOURCED** (overfunded with healthy pipelines).

**Table 2. GMI with pipeline cross-flag (base funding, attributable deaths 2019).**

| Pathogen | Funding % | Burden % | f-b (pp) | Pipeline % | Quadrant |
|---|---|---|---|---|---|
| *K. pneumoniae* | 8.4 | 20.8 | –12.4 | 11.8 | NEGLECTED |
| *E. coli* | 16.8 | 23.6 | –6.8 | 25.0 | SERVED DESPITE LOW FUNDING |
| *A. baumannii* | 8.1 | 14.2 | –6.1 | 5.3 | NEGLECTED |
| *S. pneumoniae* | 19.1 | 13.1 | +5.9 | 2.6 | INVESTMENT-TRANSLATION GAP |
| *S. aureus* | 26.8 | 19.2 | +7.6 | 30.3 | WELL-RESOURCED |
| *P. aeruginosa* | 20.9 | 9.1 | +11.7 | 25.0 | WELL-RESOURCED |

**[FIGURE 1 placeholder.** Burden share (x-axis) vs funding share (y-axis), per pathogen; bubble size encodes late-stage pipeline share; bubble colour encodes quadrant. Diagonal line shows perfect alignment. Pathogens above the diagonal are overfunded; below are underfunded. *A. baumannii* and *K. pneumoniae* sit in the bottom-left underfunded zone with the smallest bubbles. **End placeholder.]**

**[FIGURE 2 placeholder.** Per-pathogen tri-bar chart: side-by-side bars for funding share, burden share, and pipeline share for each of the six pathogens, sorted by misalignment. Visually conveys the gap structure for readers who do not parse Figure 1. **End placeholder.]**

## 4. Discussion and policy implications

The doubly-neglected finding is the central message: **two of the six leading AMR pathogens (*A. baumannii* and *K. pneumoniae*) are receiving disproportionately little funding relative to the deaths they cause, and that gap is compounded by a thin late-stage pipeline.** These two pathogens are precisely the WHO Priority 1 critical-priority carbapenem-resistant Gram-negatives. The convergence of three independent signals (burden share, funding share, pipeline share) on the same conclusion strengthens it materially.

**Implications for R&D allocation.** A reallocation on the order of 10 to 15 percentage points of within-the-six funding share toward *A. baumannii* and *K. pneumoniae* would close the misalignment for these pathogens. In absolute terms, applied to the USD 4.17 billion of pathogen-named funding in our base FAP, that is roughly USD 400 to 600 million of cumulative R&D effort. Even the most generous redistribution assumption (the equal broad-pool sensitivity, which credits broad Gram-negative funding fully to specific Gram-negative pathogens) leaves a GMI of 15.7%, so the qualitative case for reallocation does not depend on a particular attribution rule.

**Limitations and caveats stated openly.** The funding side is at the genus level (the Hub does not encode species or drug class); we use the genus as a proxy for the priority species. For five of six pathogens this is a close fit; for *Streptococcus* it is the loosest (the genus includes group A and B streptococci that are not the burden-driving pneumoniae). The funding side is cumulative across years (1975-2026 in raw data, dominated by post-2013 activity); burden is the 2019 cross-section. Project counts proxy candidate counts in the pipeline measure: one drug can be funded by multiple projects, biasing all pathogens similarly and preserving relative comparability. The Vivli surveillance datasets carry no patient outcomes, so we rely on external GRAM attributable estimates rather than computing a clinical burden score directly from the Register; the SPIDAAR linkage dataset was requested but unavailable at submission time. We frame these honestly: the misalignment finding is robust to plausible variation in attribution rules and metric choice, but the precision of the percentage-point gaps depends on the resolution at which the underlying data systems track pathogens.

**Future extensions.** A Bayesian projection of resistance prevalence to 2030 (Stan model infrastructure is in place in the repository) would yield a projected GMI showing whether the gap widens. Regional disaggregation requires regional funding attribution, which the current Hub data does not robustly support. Drug-class resolution requires extending the Hub taxonomy or linking to a separate pipeline-tracking dataset such as the WHO antibacterial agents preclinical and clinical pipelines.

**Policy ask.** The Global AMR R&D Hub and funders coordinating through it should consider an explicit alignment audit at the pathogen level on a recurring basis, with portfolio rebalancing toward the doubly-neglected pathogens identified here. The framework presented is reproducible (full code and decision log open), extensible to additional pathogens, and updatable as new GRAM and Hub data are released. The Vivli AMR Register, by providing the granular surveillance backbone, makes this kind of alignment audit possible for the first time.

## References (compressed)

Murray CJL et al. Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *Lancet* 399:629-655 (2022).
European Committee on Antimicrobial Susceptibility Testing. EUCAST Clinical Breakpoint Tables v14.0 (2024). www.eucast.org
World Health Organization. WHO Bacterial Priority Pathogens List, 2024. Geneva: WHO (2024).
World Health Organization. 2023 Antibacterial Agents in Clinical and Preclinical Development. Geneva: WHO (2024).
Vivli AMR Register. amr.vivli.org
Global AMR R&D Hub. Dynamic Dashboard. globalamrhub.org

---

*Repository (code, decision log, reproducibility): github.com/anthoniooladimeji11-coder/amr-funding-gap*

*All analyses use open data. Surveillance data accessed via the Vivli AMR Register under the 2026 Data Challenge; funding data from the Global AMR R&D Hub; burden estimates from the GRAM 2019 Lancet appendix Table S22.*
