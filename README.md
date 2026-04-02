# UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for UK Biobank RAP Data

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/Hinna0818/UKBAnalytica?style=flat)](https://github.com/Hinna0818/UKBAnalytica/stargazers)
[![GitHub last commit](https://img.shields.io/github/last-commit/Hinna0818/UKBAnalytica)](https://github.com/Hinna0818/UKBAnalytica/commits/main)
[![Visits](https://hits.sh/github.com/Hinna0818/UKBAnalytica.svg)](https://hits.sh/github.com/Hinna0818/UKBAnalytica/)
<!-- badges: end -->

<img src="man/figures/logo.png" align="right" height="139" alt="UKBAnalytica logo" />

**UKBAnalytica** is a high-performance R package for processing UK Biobank (UKB)
Research Analysis Platform (RAP) data exports. It focuses on standardized
phenotyping, survival-ready datasets, and scalable preprocessing.

**For details, please visit**: [Full documentation for UKBAnalytica](https://hinna0818.github.io/UKBAnalytica/)

![](docs/image.png)


## Installation
You can install the development version of `UKBAnalytica` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Hinna0818/UKBAnalytica")
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
  primary_disease = "AA"
)

head(analysis_dt[, .(
  eid,
  AA_history,
  Hypertension_history,
  Diabetes_history,
  outcome_status,
  outcome_surv_time
)])
```

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

## Supplementary Materials
Here we provide some learning materials for UK Biobank in which you may be interested:
- [UK Biobank database browser](https://biobank.ndph.ox.ac.uk/ukb/index.cgi)
- [UK Biobank RAP platform](https://ukbiobank.dnanexus.com/landing)
- [UK Biobank learning guides supported by our team](https://hinna0818.github.io/Bioinfo-SMU/Epidemiology/UK_Biobank/) 


