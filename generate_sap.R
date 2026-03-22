# =============================================================================
# Statistical Analysis Plan (SAP)
# Phase III RCT — Overall Survival and PFS Analysis
# NICE HTA Submission Format
# =============================================================================

sap <- "
# Statistical Analysis Plan
## Phase III RCT: Treatment X vs Best Supportive Care
### Version 1.0 | Pre-database lock

---

## 1. Study Overview

**Design:** Phase III, two-arm, open-label, randomised controlled trial

**Population:** Adults with advanced solid tumour, ECOG PS 0-2,
1-3 prior lines of systemic therapy

**Randomisation:** 1:1 stratified by ECOG PS (0 vs 1-2) and
prior lines (1 vs 2-3)

**Primary endpoint:** Overall Survival (OS)

**Secondary endpoints:** Progression-Free Survival (PFS),
Overall Response Rate (ORR), Duration of Response (DoR)

---

## 2. Primary Analysis: Overall Survival

### 2.1 Kaplan-Meier Analysis
- KM survival curves with 95% CI (Greenwood formula)
- Median OS with 95% CI (Brookmeyer-Crowley method)
- OS rates at 6, 12, 18, 24 months with 95% CI
- Number at risk tables

### 2.2 Primary Comparison
- Two-sided stratified log-rank test, alpha = 0.05
- Stratified by randomisation factors

### 2.3 Cox Proportional Hazards
- Univariable: treatment arm only
- Multivariable: arm + age + ECOG PS + prior lines + PD-L1 + biomarker
- Both report HR with 95% CI and p-value

### 2.4 PH Assumption Testing
- Schoenfeld residuals (cox.zph)
- Log(-log(S(t))) plots
- If violated: time-varying coefficients or RMST as sensitivity analysis

### 2.5 Parametric Models (for HTA extrapolation)
Six distributions fitted: exponential, Weibull, log-normal,
log-logistic, Gompertz, generalised gamma.
Selection via AIC/BIC. Extrapolation plots produced for NICE submission.

### 2.6 RMST Analysis
Restricted Mean Survival Time at tau = 24 months.
Provides summary measure robust to PH violation.

---

## 3. Subgroup Analyses (pre-specified, exploratory)

| Subgroup | Levels |
|---|---|
| Age | < 65 vs >= 65 |
| ECOG PS | 0 vs 1-2 |
| PD-L1 | Positive vs Negative |
| Prior lines | 1 vs 2+ |

Results in forest plot. No formal interaction testing.

---

## 4. Sensitivity Analyses

| Analysis | Rationale |
|---|---|
| Per-Protocol Set | Assess protocol deviation impact |
| Landmark at 6 months | Remove guarantee-time bias |
| RMST | Robust to PH violation |
| Censoring at crossover | Address post-progression treatment effect |

---

## 5. Missing Data

- < 5% missingness: complete case analysis
- >= 5% missingness: multiple imputation (m=20, PMM)
- Sensitivity: best/worst case imputation

---

## 6. Software

R (>= 4.3.0): survival, flexsurv, survminer, survRM2, ggplot2
All code version-controlled. Random seeds documented.

---

## 7. Reporting Standards

- CONSORT 2010
- NICE DSU Technical Support Documents
- ISPOR Good Practices for Outcomes Research

---

*SAP v1.0 | Author: Atrija Haldar | Finalised prior to database lock*
"

writeLines(sap, file.path(output_dir, "statistical_analysis_plan.md"))
cat("  Saved: output/statistical_analysis_plan.md\n\n")
