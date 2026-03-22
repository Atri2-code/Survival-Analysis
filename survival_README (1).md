# Advanced Clinical Trial Survival Analysis Pipeline

A comprehensive R pipeline implementing survival analysis methods used in
NICE health technology assessment (HTA) submissions and medical communications —
from raw patient-level RCT data through to publication-ready outputs and
a pre-specified statistical analysis plan.

---

## Methods implemented

| Method | Purpose |
|---|---|
| Kaplan-Meier | Non-parametric survival estimation with risk tables |
| Log-rank test | Between-arm comparison, stratified |
| Cox PH (univariable) | Unadjusted hazard ratio |
| Cox PH (multivariable) | Adjusted HR (age, ECOG PS, prior lines, PD-L1, biomarker) |
| Schoenfeld residuals | Proportional hazards assumption testing |
| Parametric models (×6) | Exponential, Weibull, log-normal, log-logistic, Gompertz, generalised gamma |
| AIC/BIC model selection | Best-fit parametric model for extrapolation |
| Extrapolation plots | Long-term survival projection for HTA cost-effectiveness models |
| RMST | Restricted Mean Survival Time — robust to PH violation |
| Landmark analysis | Removes guarantee-time bias at 6-month landmark |
| Forest plot | Subgroup HRs across 9 pre-specified subgroups |
| SAP document | Pre-specified Statistical Analysis Plan in NICE format |

---

## Endpoints analysed

- **Overall Survival (OS)** — primary endpoint
- **Progression-Free Survival (PFS)** — secondary endpoint

---

## Simulated trial characteristics

- N = 500 patients (1:1 randomisation)
- Treatment HR (OS) = 0.65 — 35% reduction in hazard
- Treatment HR (PFS) = 0.55 — 45% reduction in hazard
- Median OS (control) = 14 months
- Administrative censoring at 36 months
- Covariates: age, ECOG PS, prior lines, PD-L1 status, biomarker

---

## Outputs

| File | Description |
|---|---|
| `output/km_os.png` | OS KM curves with risk table and p-value |
| `output/km_pfs.png` | PFS KM curves with risk table |
| `output/ph_test.png` | Schoenfeld residuals plots |
| `output/survival_extrapolation.png` | Parametric extrapolation comparison |
| `output/landmark_analysis.png` | Landmark KM from 6-month timepoint |
| `output/forest_plot.png` | Subgroup HR forest plot |
| `output/parametric_model_fit.csv` | AIC/BIC comparison across distributions |
| `output/cox_results.csv` | Multivariable Cox regression results |
| `output/subgroup_analysis.csv` | Subgroup HR estimates |
| `output/baseline_table.csv` | Baseline characteristics by arm |
| `output/statistical_analysis_plan.md` | SAP in NICE submission format |

---

## How to run

```r
# Install packages
install.packages(c("survival", "survminer", "ggplot2", "dplyr",
                   "flexsurv", "broom", "gridExtra", "survRM2"))

# Run full pipeline
source("analysis.R")
```

---

## Relevance to NICE submissions

Parametric modelling and extrapolation are required in NICE Single Technology
Appraisals (STAs) to estimate long-term survival beyond trial follow-up for
cost-effectiveness modelling. This pipeline implements all six distributions
recommended in NICE DSU Technical Support Document 14.

The SAP format follows NICE DSU conventions and CONSORT 2010 reporting
standards — the same framework used in commercial HTA submissions.

---

## Dependencies

```r
survival    # Core survival analysis
survminer   # KM plots with risk tables
flexsurv    # Parametric survival models
survRM2     # Restricted Mean Survival Time
ggplot2     # Visualisation
dplyr       # Data manipulation
broom       # Tidy model outputs
```

---

## Author

Atrija Haldar
[LinkedIn](https://www.linkedin.com/in/atrija-haldar-196a3b221/)
MSc Engineering, Technology and Business Management — University of Leeds
