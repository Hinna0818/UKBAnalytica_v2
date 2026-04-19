# UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for UK Biobank RAP Data

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/Hinna0818/UKBAnalytica_v2?style=flat)](https://github.com/Hinna0818/UKBAnalytica_v2/stargazers)
[![GitHub last commit](https://img.shields.io/github/last-commit/Hinna0818/UKBAnalytica_v2)](https://github.com/Hinna0818/UKBAnalytica_v2/commits/main)
[![Visits](https://hits.sh/github.com/Hinna0818/UKBAnalytica_v2.svg)](https://hits.sh/github.com/Hinna0818/UKBAnalytica_v2/)
<!-- badges: end -->

<img src="man/figures/logo.png" align="right" height="139" alt="UKBAnalytica logo" />

**UKBAnalytica** is a high-performance R package for processing UK Biobank 
Research Analysis Platform (RAP) data exports. It focuses on standardized
phenotyping, survival-ready datasets, scalable preprocessing, and downstream analysis.

**For details, please visit**: [Full documentation for UKBAnalytica](https://hinna0818.github.io/UKBAnalytica_v2/)

![](docs/image.png)


## Installation
You can install the development version of `UKBAnalytica` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Hinna0818/UKBAnalytica_v2")
```

Sometimes due to the network problem, it is not easy to use `devtools` to install, so you can install in this way:
```{r}
# install.packages("pak")
pak::pkg_install("Hinna0818/UKBAnalytica_v2")
```

Or just clone this repo and intall it locally:
```{bash}
git clone https://github.com/Hinna0818/UKBAnalytica_v2.git
cd UKBAnalytica
R CMD INSTALL .
```

## Quick start

```r
library(UKBAnalytica)
library(data.table)

ukb_data <- fread("population.csv")

diseases <- get_predefined_diseases()[
  c("AA", "Hypertension", "Diabetes")
]

analysis_dt <- build_survival_dataset(
  dt = ukb_data,
  disease_definitions = diseases,
  prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
  outcome_sources = c("ICD10", "ICD9", "Death"),
  primary_disease = "AA",
  show_flow = TRUE,
  dt_threads = 8
)

head(analysis_dt[, .(
  eid,
  AA_history,
  Hypertension_history,
  Diabetes_history,
  outcome_status,
  outcome_surv_time
)])

# Optional: retrieve participant flow table printed by show_flow
flow_dt <- attr(analysis_dt, "participant_flow")
if (!is.null(flow_dt)) print(flow_dt)
```

## Recent survival updates (v0.6.2)

- Added optional `show_flow` in `build_survival_dataset()` to print participant attrition in terminal and attach a reusable flow table via `attr(result, "participant_flow")`.
- Added optional `dt_threads` in `build_survival_dataset()` to temporarily control `data.table` threads for large runs, with automatic restoration on exit.
- Added algorithm-source column compatibility for both `p{field}_i0` and `p{field}` naming styles.
- Improved date robustness in ICD/self-report/death parsing to prevent malformed date values from interrupting full-pipeline execution.

## What this package covers

### Core Functionality
- RAP data download helpers (Python scripts)
- Baseline preprocessing with standardized mappings
- Multi-source disease definitions (ICD-10, ICD-9, self-report, death)
- Survival analysis datasets with prevalent/incident classification
- Baseline Table 1 summaries and multiple imputation

### Advanced Analysis Modules (v0.5.0+ supports)
- **Subgroup Analysis**: Stratified analysis with interaction p-values
- **Propensity Score Methods**: PSM matching and IPTW weighting
- **Mediation Analysis**: Causal mediation using regmedint backend
- **MI Pooling**: Multiple imputation result combining (Rubin's Rules)
- **Sensitivity Analysis Preprocessing**: Exclude early events or rows with missing covariates before regression
- **Machine Learning**: Unified ML interface with SHAP interpretation
- **Visualization**: Forest plots, K-M curves, balance plots, SHAP plots

## Machine Learning Example

```r
library(UKBAnalytica)

# Train a Random Forest classifier
ml_rf <- ukb_ml_model(
  diabetes ~ age + bmi + sbp + smoking,
  data = ukb_data,
  model = "rf",
  task = "classification",
  seed = 42
)

# Model evaluation
print(ml_rf)  # AUC, accuracy, etc.
ukb_ml_metrics(ml_rf, ci = TRUE)

# ROC curve
roc <- ukb_ml_roc(ml_rf)
plot(roc)

# SHAP interpretation (requires fastshap)
shap <- ukb_shap(ml_rf, sample_n = 1000)
plot_shap_summary(shap)
plot_shap_dependence(shap, feature = "age")

# Model comparison
ml_xgb <- ukb_ml_model(diabetes ~ ., data, model = "xgboost")
comparison <- ukb_ml_compare(ml_rf, ml_xgb)
plot(comparison)

# Survival ML
surv_rf <- ukb_ml_survival(
  Surv(time, event) ~ age + sex + bmi,
  data = ukb_data,
  model = "rsf"
)
print(surv_rf)  # C-index
```

## Advanced Analysis Example

```r
# Subgroup analysis
results <- run_subgroup_analysis(
  data = dt, exposure = "treatment", outcome = "event",
  subgroup_var = "age_group", model_type = "cox",
  endpoint = c("time", "status")
)
plot_forest(results)

# Multiple imputation pooling
pooled <- pool_mi_models(
  datasets = mi_datasets,
  formula = Surv(time, status) ~ treatment + age + sex,
  model_type = "cox"
)
summary(pooled)

# Mediation analysis
med <- run_mediation(
  data = dt, exposure = "treatment", mediator = "biomarker",
  outcome = "event", outcome_type = "cox"
)
plot_mediation(med, type = "effects")
```

## Sensitivity Analysis Example

```r
# Remove events occurring in the first 2 years of follow-up
dt_sens1 <- sensitivity_exclude_early_events(
  data = analysis_dt,
  endpoint = c("outcome_surv_time", "outcome_status"),
  n_years = 2
)

# Remove rows with any missing adjustment covariate
dt_sens2 <- sensitivity_exclude_missing_covariates(
  data = dt_sens1,
  covariates = c("age", "sex", "bmi", "smoking")
)

# Pass directly to the standard regression interface
cox_sens <- runmulti_cox(
  data = dt_sens2,
  main_var = c("bmi", "sbp"),
  covariates = c("age", "sex", "bmi", "smoking"),
  endpoint = c("outcome_surv_time", "outcome_status")
)
```

## Supplementary Materials
Here we provide some learning materials for UK Biobank in which you may be interested:
- [UK Biobank database browser](https://biobank.ndph.ox.ac.uk/ukb/index.cgi)
- [UK Biobank RAP platform](https://ukbiobank.dnanexus.com/landing)
- [UK Biobank learning guides supported by our team](https://hinna0818.github.io/Bioinfo-SMU/Epidemiology/UK_Biobank/) 

