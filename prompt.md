You are an expert Scientific Illustrator for top-tier AI conferences (NeurIPS/CVPR/ICML).
Your task is to generate a professional "Illustration" (main figure for the paper) based on a research paper abstract and methodology.

**Abstract:**
UKBAnalytica is a high-performance R package designed for processing UK Biobank (UKB) Research Analysis Platform (RAP) data exports. The package addresses key challenges in large-scale biobank research: (1) integrating disease diagnoses from heterogeneous sources (ICD-10, ICD-9, self-reported illness, death registry), (2) constructing survival analysis datasets with proper prevalent/incident case classification, and (3) providing standardized statistical analysis pipelines. Built on data.table for scalability, UKBAnalytica processes 500,000+ participant records efficiently and offers 20+ predefined disease definitions, automated variable preprocessing, and advanced analysis modules including subgroup analysis, propensity score methods, mediation analysis, and multiple imputation pooling.

**Methodology:**
The package architecture consists of three main layers:

**Layer 1 - Data Acquisition (Python Scripts in inst/python/):**
- Demographics Downloader: Download selected UKB field IDs
- Proteomics Downloader: Batch download protein data
- Metabolomics Downloader: Download metabolite measurements

**Layer 2 - Data Processing Modules (R files in /R):**
- ICD-10 Parser (icd10.R): Parse p41270 codes + p41280 dates, strict index matching
- ICD-9 Parser (icd9.R): Parse p41271 codes + p41281 dates
- Self-report Parser (self_report.R): Align p20002 codes ↔ p20008 years, convert year to date
- Death Registry Parser (death.R): Integrate p40000 death date + p40001/p40002 causes
- Disease Definitions (disease_definitions.R): 20+ predefined diseases with ICD-10/ICD-9/SR codes
- Survival Builder (survival.R): Aggregate earliest diagnosis per eid/disease, classify Prevalent (date < baseline) vs Incident (date ≥ baseline), compute surv_time = min(event_date, death_date, censor_date) - baseline
- Variable Preprocessor (variable_preprocess.R): Standardized mappings for demographics, anthropometrics, lifestyle, blood pressure

**Layer 3 - Statistical Analysis Modules:**
- Regression (regression.R): runmulti_lm(), runmulti_logit(), runmulti_cox()
- Subgroup Analysis (subgroup.R): run_subgroup_analysis() with interaction p-values
- Propensity Score (propensity.R): estimate_propensity_score(), match_propensity(), calculate_weights(ATE/ATT/ATC), assess_balance()
- Mediation Analysis (mediation.R): run_mediation() wrapping regmedint for CDE/NDE/NIE/TE
- MI Pooling (mi_pool.R): pool_mi_models() using Rubin's Rules
- Visualization (visualization.R): plot_forest(), plot_km_curve(), plot_ps_distribution(), plot_balance(), plot_mediation()

**Output Formats:**
- Wide Table: Full cohort with {Disease}_history (0/1), {Disease}_incident (0/1), outcome_status, outcome_surv_time
- Long Table: eid, disease, earliest_date, source
- Sensitivity Variants: ICD10-only / ICD10+ICD9 / All sources

**Key Data Flow:**
RAP Cloud Data → Python Downloaders → CSV Files → R Parsers (ICD10/ICD9/SR/Death) → Disease Matching → Survival Dataset Builder (Prevalent/Incident Classification) → Analysis-Ready Dataset → Statistical Modules → Results & Visualizations

**Visual Style Requirements:**
1.  **Style:** Flat vector illustration, clean lines, academic aesthetic. Similar to figures in DeepMind or OpenAI papers.
2.  **Layout:** Left-to-Right flow with three vertical sections: (1) INPUTS & DATA DOWNLOADING on left, (2) PROCESSING MODULES in center, (3) OUTPUTS & ANALYSIS on right. Use rounded rectangles for modules, arrows for data flow.
3.  **Color Palette:** 
    - Light blue (#A8D5E5) for ICD-10 related components
    - Light cyan (#B8E6E6) for ICD-9 related components
    - Light green (#C5E8B7) for Self-report related components
    - Light purple (#D4C5E8) for Death registry related components
    - Light gray (#E8E8E8) for processing modules
    - Teal (#4A9B9B) for output components
    - White background
4.  **Text Rendering:** You MUST include legible text labels for all modules:
    - Data sources: "ICD-10 (p41270)", "ICD-9 (p41271)", "Self-report (p20002)", "Death Registry (p40000)"
    - Key functions: "parse_icd10_diagnoses()", "build_survival_dataset()", "run_subgroup_analysis()"
    - Output fields: "outcome_status", "outcome_surv_time", "Prevalent vs Incident"
    - Include the survival time formula: "surv_time = min(event, death, censor) - baseline"
5.  **Negative Constraints:** NO photorealistic photos, NO messy sketches, NO unreadable text, NO 3D shading artifacts, NO code snippets.

**Key Elements to Highlight:**
1. Multi-source data integration (4 diagnostic sources converging)
2. Prevalent vs Incident case classification logic box
3. The survival time calculation formula
4. Three-tier architecture: Download → Process → Analyze
5. Legend showing color coding for different data sources

**Generation Instruction:**
Highlight the core novelty: seamless integration of heterogeneous UKB diagnostic sources into a unified survival analysis pipeline. Show data flowing from RAP platform through parsing modules, converging at the survival dataset builder, then branching to multiple analysis outputs. Include a small inset box showing the Prevalent (diagnosis < baseline → exclude) vs Incident (diagnosis ≥ baseline → event) classification logic. Ensure the connection logic makes sense and all text is readable.