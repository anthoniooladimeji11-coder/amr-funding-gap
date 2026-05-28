# Project Decision Log — AMR Burden–Funding Gap (GFAO / VALOR)

A running record of decisions, findings, and rationale. Newest entries at the
top of each section. Every entry should answer: *what* we decided, and *why*.
This log exists so we can defend choices to the judging panel, stay reproducible,
and keep both teammates in sync.

---

## Project at a glance

- **Goal:** Quantify and map the *misalignment* between where antimicrobial
  resistance (AMR) burden concentrates and where global R&D funding flows;
  provide a counterfactual reallocation tool (VALOR dashboard).
- **Why not causal:** A naive "funding suppresses resistance" causal claim is
  unidentifiable here (temporal mismatch, reverse causality, global/product-
  directed funding that doesn't map to country outcomes). We therefore do a
  **descriptive + counterfactual optimization** design, not causal inference.
- **Award target:** AMR Visionary Award + Global AMR R&D Hub Cross-Domain Award.
- **Team:** 2 people. One codes (lead), one supports (clinical/writing).
- **Repo:** github.com/anthoniooladimeji11-coder/amr-funding-gap (public, MIT).

---

## Pipeline (5 steps)

1. **Harmonize + calibrate** surveillance data; reconstruct representative
   baselines (Deville–Särndal calibration raked to WHO GLASS, hierarchical
   fallback for sparse cells).
2. **Latent Class Analysis** (stepmix, Python) to group isolates into
   multidrug-resistance phenotypes (e.g. XDR-CRAB) as single metrics.
3. **Hierarchical Bayesian Multivariate Probit** (Stan via cmdstanpy), scoped to
   *K. pneumoniae* and *A. baumannii* vs. cefiderocol / omadacycline /
   carbapenems; project resistance to 2030, pooling toward sparse LMIC settings.
4. **Global Misalignment Index (GMI):** compare region-specific, burden-weighted
   resistance against pathogen(–drug-class) funding shares; cross-flag cells
   against pipeline candidates to separate genuine neglect from rational
   deprioritization.
5. **VALOR dashboard** (Streamlit): gap map + reallocation counterfactual engine,
   with modeled core clearly demarcated from illustrative layers.

---

## Key decisions

### D-008 — Censored MICs preserved with numeric + flag (2026-05-28)
**Decision:** Store each MIC three ways: `mic` (raw string), `mic_numeric`
(parsed number), `mic_censor` (`<=`, `<`, `>=`, `>`, `exact`, `missing`).
**Why:** `<=0.004` is left-censored (true value at/below assay floor), not
literally 0.004. Keeping the flag lets the breakpoint step and log2 transform
treat censoring correctly. Surfaced as a pyarrow write error on KEYSTONE.

### D-007 — Language: Python-only (2026-05-28)
**Decision:** Drop the R/poLCA part of the original plan; do everything in
Python. LCA via `stepmix` (categorical mixture model, sklearn-compatible);
Bayesian via Stan/`cmdstanpy` (PyMC/numpyro as fallback).
**Why:** Solo coder + 8-week clock. Two-language pipelines add a fragile R↔Python
hand-off and double the debugging surface. stepmix reproduces poLCA's method.
**Write-up impact:** "LCA via stepmix" instead of "poLCA".

### D-006 — Clinical Burden Score re-grounded on external estimates (2026-05-28)
**Decision:** Drop the SPIDAAR-based mortality × length-of-stay CBS. Re-ground
burden weighting on external open data (IHME/GRAM AMR burden, WHO GLASS), keyed
by pathogen × region.
**Why:** Data inspection confirmed NONE of the three available datasets carry
outcomes (mortality/LOS) fields. SPIDAAR was requested but not received/available.
**Upside:** External burden estimates are more globally credible than 4-country
SPIDAAR data and sidestep the transportability problem.
**Open item:** confirm exact external source + vintage; document its cells.

### D-005 — GMI granularity: pathogen-level on the funding side (PENDING) (2026-05-28)
**Issue:** Hub funding data has NO clean "drug class" column (only Infectious
Agent / Individual Infectious Agent, Research Area/Subcategory, Sector, Product
Name). The GMI was formulated per pathogen × drug class.
**Leaning:** Compute GMI at **pathogen level** on the funding side (burden side
can stay drug-class-resolved from MICs), stating the asymmetry explicitly.
**Status:** NOT finalized — need to inspect actual values in Hub columns
(Individual Infectious Agent, Research Subcategory) before locking. See TODO.

### D-004 — Keep full data + scoped subset (2026-05-28)
**Decision:** Harmonize ALL pathogens/drugs (chunked, → Parquet), and also emit
a `core_pathogens_long.parquet` filtered to K. pneumoniae + A. baumannii.
**Why:** LCA needs the full antibiotic panel to call XDR/PDR; keeping all data
preserves the option to add pathogens later. Parquet = ~10x smaller/faster than
the 388 MB ATLAS CSV.

### D-003 — Separate ATLAS β-lactamase genes from MIC table (2026-05-28)
**Decision:** ATLAS's trailing gene-flag columns (KPC, NDM, OXA, VIM, CTX-M, …)
go into a separate isolate-level table, not the MIC long table.
**Why:** They are genotype, not MIC measurements. Useful as covariates / for
mechanism-aware phenotyping (appeals to Bradford). Also note: this confirms
ATLAS *does* carry KPC/OXA genotype — relevant to the separate "JARS" team.

### D-002 — Common long-format schema (2026-05-28)
**Decision:** Map all datasets to one long schema:
`dataset, isolate_id, organism, country, region, year, specimen, age, gender,
antibiotic, mic, mic_numeric, mic_censor, interpretation`.
**Why:** Three datasets arrive in three different wide shapes with inconsistent
names (ATLAS "Species/Year/Source" + `_I` interp cols; SIDERO "Organism Name/
Year Collected/Body Location"; KEYSTONE "Organism/Study Year/Specimen Type" with
newline-laden headers). Long format is what LCA + Bayesian models consume.

### D-001 — Descriptive + counterfactual, not causal (pre-session)
**Decision:** Frame the project as misalignment mapping + allocation
counterfactuals, not a causal effect of funding on resistance.
**Why:** Causal identification fails here (see "Why not causal" above). This was
the central pivot that made the proposal defensible.

---

## Data inventory (as received)

| Dataset | File | Size | Rows×Cols | Range | Notes |
|---|---|---|---|---|---|
| ATLAS | atlas_vivli_2004_2024.csv | 388 MB | many × 127 | 2004–2024 | wide; drug+`_I` pairs; β-lactamase gene flags |
| KEYSTONE | Omadacycline_...xlsx | 19 MB | × 46 | 2015–2025 | omadacycline; rich patient context; newline headers |
| SIDERO-WT | Updated_Shionogi...SIDERO-WT...xlsx | 5.7 MB | × 20 | 2013–2019 | cefiderocol; MIC only (no interp col) |
| Hub funding | Projects.xlsx | 60 MB | 18,854 × 26 | 2017–? | funder/agent/sector/amount; NO drug-class col |

**Confirmed absent across all surveillance data:** mortality, length-of-stay,
any clinical outcome. (Drove decision D-006.)

**SIDERO content check (verified clean):** 388,018 long rows; organisms &
antibiotic names normalize correctly; core pathogens = 95,352 rows; MICs are
real, some censored.

---

## Environment

- macOS (Apple Silicon, arm64), Python 3.13.9 in `.venv`.
- Core: pandas, numpy, scipy, scikit-learn, matplotlib, openpyxl, pyarrow.
- LCA: stepmix. Bayesian: cmdstanpy + CmdStan 2.39.0 (compiled OK).
- Dashboard (later): streamlit.
- Pinned in `requirements.txt`.
- **Data NEVER committed** — `data/` and all spreadsheet formats gitignored.

---

## Open TODOs / risks to resolve

- [ ] **Inspect Hub funding column VALUES** (Individual Infectious Agent,
      Research Subcategory, Sector) to finalize D-005 (GMI granularity).
- [ ] Decide funding-year attribution rule (Start Year vs spread vs End Year).
- [ ] Choose + document external burden source for CBS (IHME/GRAM vs GLASS).
- [ ] Pull WHO GLASS reference margins for calibration; note coverage gaps.
- [ ] Build breakpoint module (CLSI/EUCAST) → derive S/I/R, handle censoring.
- [ ] Confirm with Vivli whether SPIDAAR is available (low priority).
- [ ] Verify official rules: EOI word count basis, Cross-Domain Award value.
- [ ] 5-page final report limit (incl. figures) — design figures early.
