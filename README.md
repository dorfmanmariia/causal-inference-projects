[README.md](https://github.com/user-attachments/files/28160060/README.md)
# causal-inference-projects

Stata code for two applied microeconometrics projects. Each project uses real-world data and implements a complete causal inference pipeline — from data cleaning and validity checks through to robustness analysis and visualization.

---

## Project 1: Female Board Representation and Firm Performance

**File:** `project1_female_board_roa.do`

### Research Question
Does appointing a female director causally improve firm profitability (ROA)?

### Data
Cross-sectional panel of publicly listed firms, 2010–2019. Source: Refinitiv (board composition, financials, TRBC industry classification). Treatment: firm appoints at least one female director after 2015 corporate governance recommendations. Outcome: average ROA over 2015–2019.

### Methods
The observational setting creates selection bias — firms that appoint women may differ systematically from those that do not. The code implements three complementary strategies to address this:

| Estimator | Key assumption |
|-----------|---------------|
| Propensity Score Matching (PSM, 1:1 nearest-neighbor) | Conditional independence given observable baseline characteristics |
| Coarsened Exact Matching (CEM) | Same |
| IPWRA (Doubly Robust) | Consistent if *either* the outcome model or the propensity score model is correctly specified |

Full diagnostics included: propensity score overlap, covariate balance (standardized bias before/after matching), placebo test on pre-treatment ROA, heterogeneity by firm size, and IPW-weighted dynamic treatment effects by year.

### Pipeline
```
Step 0  Data cleaning, winsorization, panel → cross-section conversion
Step 1  Summary statistics and pre-treatment balance t-tests
Step 2  Propensity score logit with industry fixed effects
Step 3  Overlap check (kernel density plot)
Step 4  PSM + covariate balance (pstest)
Step 5  Treatment effect estimation: OLS / PSM / CEM / IPWRA
Step 6  Robustness: alternative outcome (ΔROA), placebo, size heterogeneity
Step 7  Figures: PS overlap, ROA density after matching, dynamic effects
```

### Required packages
```stata
ssc install psmatch2
ssc install cem
ssc install winsor2
ssc install estout
ssc install coefplot
```

---

## Project 2: Credit Scores, Interest Rates, and Mortgage Default

**File:** `project2_rdd_iv_mortgage.do`

### Research Question
Does a higher mortgage interest rate causally increase the probability of delinquency?

### Design
Lenders use FICO = 620 as an internal pricing threshold: borrowers just above receive lower rates than observably similar borrowers just below. This creates a **fuzzy RDD** — the threshold predicts rates but does not deterministically assign them, making it a valid instrument for 2SLS estimation.

```
First stage:   FICO ≥ 620  →  lower interest rate   (F > 100)
Reduced form:  FICO ≥ 620  →  lower delinquency
2SLS (LATE):   interest rate  →  delinquency, instrumented by FICO ≥ 620
```

### Methods

**Validity checks**
- McCrary density test (no bunching at threshold)
- Covariate balance at cutoff (9 covariates)
- Placebo cutoffs at 590, 600, 610, 630, 640, 650
- Donut RDD (excluding observations exactly at FICO = 620)
- Discontinuous missingness analysis for covariates

**Main estimation**
- Reduced form and 2SLS across three control specifications
- Polynomial orders p = 1, 2, 3 × bandwidths from ±20 to ±80
- MSE-optimal bandwidth via `rdrobust` (uniform and triangular kernels)
- Probit marginal effect as robustness to functional form

**Heterogeneity**
- Bank vs. broker-originated loans
- Racial heterogeneity (white vs. non-white borrowers)
- Refinance vs. purchase mortgages
- Formal interaction tests

**Bandwidth sensitivity**
- 2SLS coefficient plotted across bandwidths 25–60 with 95% CIs
- Local randomization tests in narrow windows (±5, ±10)

### Pipeline
```
Step 0  Data loading, diagnostics, logical consistency checks
Step 1  Visual diagnostics: binned scatter, FICO histogram, first stage table
Step 2  RDD validity: McCrary, covariate balance, placebo cutoffs, donut RDD
Step 3  Main results: reduced form + 2SLS (Table 4)
Step 4  Heterogeneity analysis (Table 5)
Step 5  Bandwidth sensitivity graph + local randomization tests
```

### Required packages
```stata
ssc install estout
ssc install rdrobust
ssc install rddensity
ssc install lpdensity
net install rdlocrand, from(https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata)
```

---

## Data

The datasets are not included in this repository as they were provided for coursework. The `.do` files expect:
- `data_ps2.dta` — firm-level panel data (Project 1)
- `data_ps4.dta` — mortgage-level cross-section (Project 2)

---
