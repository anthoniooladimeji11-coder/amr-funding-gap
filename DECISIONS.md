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

### D-009 — Isolate identity assigned before melting (2026-05-28)
**Decision:** In the harmonizer, every source row gets a guaranteed-unique
`isolate_id` BEFORE the wide→long melt: real id where present
(`atlas:<id>`, `keystone:<id>`), else source row position (`sidero:row<N>`),
always dataset-prefixed.
**Why:** After melting, isolate identity cannot be recovered. SIDERO has no
native id column; without this fix, grouping its drug-rows back into isolates
was impossible (a metadata-tuple approximation would wrongly merge distinct
isolates). This fix makes per-isolate pivoting correct for all three datasets.
**Impact:** Step 1 re-run; row counts unchanged (17.48M total, 2.74M core).

### D-010 — Age string / year Int64 across datasets (2026-05-28)
**Decision:** `age`, `gender`, `country`, `region`, `specimen`, `isolate_id`
forced to string; `year` coerced to nullable Int64.
**Why:** Datasets encode age incompatibly (ATLAS uses bands like '61+';
KEYSTONE uses numeric age). Mixed types broke the Parquet concat. Age is a
covariate, not used in arithmetic, so string is lossless and safe.

### D-011 — Breakpoints: EUCAST v14.0 only, verified against source PDF (2026-05-28)
**Decision:** S/I/R classification uses EUCAST Clinical Breakpoint Tables
v14.0 (2024) ONLY. CLSI dropped entirely.
**Why:** EUCAST tables are free/official and were verified row-by-row against
the actual PDF. CLSI M100 could not be accessed for verification; rather than
enter unverifiable (and partly wrong) search-derived CLSI values, CLSI is
omitted and documented as a limitation. EUCAST is the sole primary standard.
**Corrections caught vs. search values:** Klebsiella ciprofloxacin was NOT
0.001 (real ~S≤0.25/R>0.5); Acinetobacter ciprofloxacin IS 0.001 (genuine);
Acinetobacter amikacin/gentamicin DO have (bracketed) breakpoints (search
wrongly said none); Acinetobacter cefiderocol has NO clinical breakpoint
(IE / PK-PD cutoff only); doripenem has real published values (not "derived").
**Table:** 25 verified EUCAST rows + 7 confirmed no-breakpoint, committed to
repo at data/external/breakpoints.csv (reference data, force-added past
.gitignore). Panel widened with tobramycin + aztreonam for XDR category coverage.

### D-012 — Censored-MIC-aware classification rule (2026-05-28)
**Decision:** Classification respects censoring:
'<=v' → S only if v ≤ S_max else Unresolved(U); '>v' → R only if v ≥ R_min
else U; exact → standard rule. R_min set to next doubling-dilution above the
EUCAST "R>" value so `mic >= R_min` reproduces EUCAST "R if mic > Rval".
**Why:** '<=0.5' is not literally 0.5; forcing a call on a censored value that
straddles the breakpoint would be wrong. U is honest uncertainty, not error.
**Result on core pathogens:** of ~997k panel-drug calls, ~95% definitive
(665k S, 263k R, 25k I); 4.4% U, concentrated in newer agents (cef-avi,
cef-taz, cefiderocol) where censoring/edge values legitimately can't resolve.

### D-013 — LCA design (Step 2) (2026-05-28)
**Decision:** Latent Class Analysis via stepmix `binary_nan`, run SEPARATELY
per pathogen, binary encoding (S=0; I/R=1; U=missing), missingness handled by
full-information ML, fixed random_state=42, n_init=10. Clustering panels chosen
by VARIATION (drugs near 0% or 100% NS, or very sparse, excluded):
  - A. baumannii: amikacin, gentamicin, meropenem, imipenem, levofloxacin,
    colistin (6 drugs); min 3 drugs tested per isolate → 45,012 isolates.
  - K. pneumoniae: meropenem, imipenem, doripenem, ertapenem, amikacin,
    gentamicin, ciprofloxacin, levofloxacin, aztreonam, colistin (10 drugs);
    min 4 drugs → 73,824 isolates.
**Why separate panels/pathogens:** different breakpoint availability and
different resistance biology; pooling would just separate species. Thresholds
chosen at the "cliff" in the isolate-survival curve (keeps large, well-
characterized samples).
**Pre-registered expectation (from co-resistance diagnostic):** Ab ≈ 2–3
phenotypes incl. an XDR-but-colistin-sparing class; Kp ≈ 4–5 incl. susceptible
/ FQ-ESBL / carbapenem-R / XDR. (Used to validate LCA output, not bias it.)

### D-014 — LCA results & class-count choice (2026-05-28)
**A. baumannii (K=4 by BIC, mean max-prob 0.92 — confident):**
  - C0 (52%): XDR-like — ~96–100% NS to amik/gent/mero/imip/levo, colistin 4.7%
  - C1 (33%): susceptible (wild-type)
  - C2 (8%): carbapenem-R, aminoglycoside-S
  - C3 (7%): aminoglycoside+FQ-R, carbapenem-S
  Colistin correctly isolated as the spared last-line agent. Matches prediction.
**K. pneumoniae (BIC favored K=6):** classes interpretable (susceptible;
several ESBL/FQ-resistant variants; carbapenem-R; XDR with colistin 19%; a tiny
425-isolate carbapenem-R/FQ-S group). **OPEN:** whether K=6 is too granular vs
K=4/5 for parsimony/defensibility — BIC kept improving, but two classes are
small variants. To be decided with clinical teammate (looking at K=4/5 profiles).

### Files added this session
- src/amr_gap/harmonize.py (updated: unique ids, censoring, type fixes)
- src/amr_gap/breakpoints.py (CLSI/EUCAST engine, censoring-aware, VERIFY-guard)
- data/external/breakpoints.csv (EUCAST v14.0 verified table — committed)
- src/amr_gap/lca_diagnostic.py, lca_threshold.py, coresistance.py (planning)
- src/amr_gap/lca.py (Step 2 LCA)
- data/processed/lca_acba.parquet, lca_klpn.parquet (labeled isolates — gitignored)

### Updated open TODOs
- [ ] DECIDE Klebsiella K (4/5/6) — see D-014.
- [ ] Inspect Hub funding column VALUES → finalize GMI granularity (D-005).
- [ ] Choose/document external burden source for CBS (IHME/GRAM vs GLASS).
- [ ] Funding-year attribution rule (start vs spread vs end year).
- [ ] Calibration weighting (Deville–Särndal raked to GLASS) — Step 1 finish.
- [ ] Bayesian projection (Step 3); GMI (Step 4); VALOR dashboard (Step 5).
- [ ] Name LCA phenotypes clinically (with teammate) for the write-up.
- [ ] Verify EOI/prize/rules details on official Vivli page.

### D-005 RESOLVED — GMI funding side is pathogen(genus)-level via Categories codes (2026-05-28)
**Evidence (hub_inspect + hub_categories):** Hub total $18.87B / 18,853 projects.
No drug-class field exists (Research Area = type: Therapeutics/Diagnostics/etc;
Research Subcategory = stage: Discovery/Development/Approval; Product Name 97%
missing). Therefore GMI funding side CANNOT be drug-class-resolved.
**Decision:** GMI is **pathogen × year** on the funding side. Drug-class
resolution dropped on funding (kept on burden side from MICs); asymmetry stated
openly.
**Attribution source:** the structured `Categories` field (coded taxonomy,
e.g. "1504 Infectious Agent / Bacteria / Gram negative / Klebsiella spp.") is
PRIMARY — more reliable than the free-text 'Individual Infectious Agent' field
(57% missing). Cross-validated against free-text: Categories 515 Kp / 469 Ab
vs free-text 509 / 414 — consistent, Categories catches slightly more.
**Genus-vs-species caveat (must state in report):** Hub attributes to GENUS
("Klebsiella spp.", "Acinetobacter spp."), surveillance is SPECIES
(K. pneumoniae, A. baumannii). K. pneumoniae / A. baumannii dominate clinical
isolates of their genera, so genus funding ≈ species proxy — declared as an
assumption, not hidden.
**Funding figures (Categories, Amount USD):**
  - Klebsiella spp.: 515 projects, ~$593M
  - Acinetobacter spp.: 469 projects, ~$610M
  - Broad "Gram negative" (unattributable pool): 5,130 projects, ~$4.90B
**Attribution rule (pre-registered, both built):**
  - BASE CASE: pathogen-named funding only (the ~$593M / ~$610M).
  - SENSITIVITY: add a proportional share of the $4.90B broad Gram-negative
    pool allocated down to each genus.

### D-015 — Funding-year attribution (OPEN, leaning Start Year) (2026-05-28)
Hub has Start Year (1975–2026) and End Year (to 2036). Leaning: attribute a
project's USD to its Start Year for the FAP (simplest, defensible), with a
sensitivity option to spread evenly across Start..End. To finalize when building
the FAP. Note funding predates surveillance (pre-2004) — will restrict FAP years
to overlap the burden window.

### D-014 RESOLVED — Klebsiella K = 5 (clinical + parsimony choice) (2026-05-28)
**Decision:** Klebsiella pneumoniae LCA fixed at **K=5 classes**, overriding the
BIC-optimal K=6. Acinetobacter remains K=4 (BIC-optimal and clean).
**Why:** Comparing K=4/5/6 profiles with the clinical teammate:
  - 4->5 is justified: it splits the ESBL/FQ-resistant class into
    aminoglycoside-RESISTANT vs aminoglycoside-SUSCEPTIBLE variants — a real,
    treatment-relevant distinction (aminoglycosides are a viable option for the
    susceptible group).
  - 5->6 only adds a rare (~0.6%, 425-isolate) carbapenem-R / FQ-S class. Genuine
    but tiny; on a 5-page submission it invites "why is this its own phenotype?"
    and isn't worth defending vs the parsimony cost.
**K=5 Klebsiella phenotypes (P(NS|class)):**
  - Susceptible (~53%): all ~0%.
  - Aminoglycoside-S ESBL (~15%): cipro 100%, levo 70%, aztreonam 61%, gent 0%, amik 1%, carbapenems 0%.
  - Aminoglycoside-R ESBL (~18%): cipro 100%, gent 98%, aztreonam 95%, levo 89%, carbapenems 0%.
  - XDR (~13%): carbapenems 96–100%, FQ 98–100%, aztreonam 96%, amik 63%, colistin 19%.
  - Small carbapenem-intermediate / aztreonam group (~1%): partial carbapenem, aztreonam 94%, FQ low.
  Assignment clean (mean max-prob 0.92).
**Implementation:** `FORCED_K = {"klebsiella pneumoniae": 5}` in lca.py; the
module still fits & reports K=2..6 BIC for transparency, then selects the forced K.
**Caveat to note in report:** the ~1% class is the least clean phenotype (may
reflect incomplete testing or a true intermediate mechanism) — do not over-interpret.
