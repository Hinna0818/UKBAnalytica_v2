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

## Installation

```r
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

### Advanced Analysis Modules (v0.5.0+)
- **Subgroup Analysis**: Stratified analysis with interaction p-values
- **Propensity Score Methods**: PSM matching and IPTW weighting
- **Mediation Analysis**: Causal mediation using regmedint backend
- **MI Pooling**: Multiple imputation result combining (Rubin's Rules)
- **Visualization**: Forest plots, K-M curves, balance plots, diagnostic plots

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

## License

MIT License © 2024 Nan He

