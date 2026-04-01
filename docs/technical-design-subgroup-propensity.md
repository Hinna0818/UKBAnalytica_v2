# UKBAnalytica 新模块技术设计文档

## 亚组分析、倾向性评分、中介分析与可视化模块

**版本**: 0.4.0 (规划)  
**作者**: UKBAnalytica 开发团队  
**日期**: 2026-04-01

---

## 目录

1. [概述](#概述)
2. [模块一：亚组分析 (subgroup.R)](#模块一亚组分析-subgroupr)
3. [模块二：倾向性评分 (propensity.R)](#模块二倾向性评分-propensityr)
4. [模块三：中介分析 (mediation.R)](#模块三中介分析-mediationr)
5. [模块四：可视化 (visualization.R)](#模块四可视化-visualizationr)
6. [模块五：多重插补结果合并 (mi_pool.R)](#模块五多重插补结果合并-mi_poolr)
7. [依赖管理](#依赖管理)
8. [函数依赖关系](#函数依赖关系)
9. [使用示例](#使用示例)
10. [开发计划](#开发计划)

---

## 概述

本文档描述五个新模块的技术设计：

| 模块 | 文件 | 功能 |
|------|------|------|
| 亚组分析 | `subgroup.R` | 分层分析与交互作用检验 |
| 倾向性评分 | `propensity.R` | PSM匹配与IPTW加权 |
| 中介分析 | `mediation.R` | 因果中介分析 (基于regmedint) |
| 可视化 | `visualization.R` | 森林图、K-M曲线、中介效应图等 |
| 多重插补合并 | `mi_pool.R` | 多重插补结果合并 (Rubin's Rules) |

这些模块将与现有的 `regression.R`、`survival.R` 等模块紧密集成，形成完整的流行病学分析工作流。

---

## 模块一：亚组分析 (subgroup.R)

### 1.1 核心函数

#### `run_subgroup_analysis()`

主函数：在不同亚组中运行回归模型，计算交互作用P值。

**函数签名**:

```r
run_subgroup_analysis(
  data,

data,
  exposure,
  outcome,
  subgroup_var,
  covariates = NULL,
  model_type = c("cox", "logistic", "linear"),
  endpoint = NULL,
  ref_level = NULL
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `data` | data.frame/data.table | 包含所有变量的数据集 |
| `exposure` | character | 暴露变量名 |
| `outcome` | character | 结局变量名 (线性/逻辑回归) |
| `subgroup_var` | character | 亚组变量名 |
| `covariates` | character vector | 协变量名称向量 |
| `model_type` | character | 模型类型: "cox", "logistic", "linear" |
| `endpoint` | character(2) | Cox模型的生存时间和状态列 c("time", "status") |
| `ref_level` | character | 参考组水平 (可选) |

**返回值**: data.frame

| 列名 | 说明 |
|------|------|
| `subgroup` | 亚组名称 |
| `n` | 样本量 |
| `n_event` | 事件数 (Cox/logistic) |
| `estimate` | HR/OR/Beta |
| `lower95` | 95%CI下限 |
| `upper95` | 95%CI上限 |
| `pvalue` | P值 |
| `p_interaction` | 交互作用P值 |

---

#### `run_multi_subgroup()`

批量运行多个亚组变量的分析。

**函数签名**:

```r
run_multi_subgroup(
  data,
  exposure,
  outcome,
  subgroup_vars,
  covariates = NULL,
  model_type = "cox",
  endpoint = NULL
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `subgroup_vars` | character vector | 多个亚组变量名 |
| 其他参数 | - | 同 `run_subgroup_analysis()` |

**返回值**: 合并的 data.frame，包含所有亚组的分析结果

---

#### `.calculate_interaction_pvalue()` (内部函数)

计算暴露与亚组变量的交互作用P值。

**函数签名**:

```r
.calculate_interaction_pvalue(
  data,
  exposure,
  outcome,
  subgroup_var,
  covariates,
  model_type,
  endpoint
)
```

---

### 1.2 技术实现要点

#### 模型类型支持

| 模型类型 | 公式 |
|----------|------|
| Cox | `Surv(time, status) ~ exposure * subgroup + covariates` |
| Logistic | `outcome ~ exposure * subgroup + covariates` (family = binomial) |
| Linear | `outcome ~ exposure * subgroup + covariates` |

#### 交互作用检验

1. 构建含交互项的完整模型
2. 提取交互项系数的 Wald 检验 P 值
3. 连续型亚组变量需先进行分组处理

#### 缺失值处理

- 默认使用完整案例分析 (complete case analysis)
- 可选参数控制缺失值行为

---

## 模块二：倾向性评分 (propensity.R)

### 2.1 核心函数

#### `estimate_propensity_score()`

计算倾向性评分。

**函数签名**:

```r
estimate_propensity_score(
  data,
  treatment,
  covariates,
  method = c("logistic", "gbm"),
  formula = NULL
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `data` | data.frame/data.table | 数据集 |
| `treatment` | character | 处理组变量名 (0/1) |
| `covariates` | character vector | 用于估计PS的协变量 |
| `method` | character | 估计方法: "logistic" 或 "gbm" |
| `formula` | formula | 可选自定义公式 |

**返回值**: data.table

| 列名 | 说明 |
|------|------|
| `eid` | 参与者ID |
| `ps` | 倾向性评分 |
| `treatment` | 原始分组 |

---

#### `match_propensity()`

倾向性评分匹配。

**函数签名**:

```r
match_propensity(
  data,
  ps_col = "ps",
  treatment,
  ratio = 1,
  caliper = 0.2,
  method = c("nearest", "optimal"),
  replace = FALSE,
  exact_match = NULL
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `data` | data.table | 含PS的数据集 |
| `ps_col` | character | PS列名 |
| `treatment` | character | 处理组变量名 |
| `ratio` | numeric | 匹配比例 (1:n) |
| `caliper` | numeric | 卡尺宽度 (PS标准差的倍数) |
| `method` | character | 匹配方法: "nearest" 或 "optimal" |
| `replace` | logical | 是否有放回匹配 |
| `exact_match` | character vector | 精确匹配变量 |

**返回值**: data.table (匹配后数据集，新增 `match_id`, `match_distance` 列)

---

#### `calculate_weights()`

计算 IPTW 权重。

**函数签名**:

```r
calculate_weights(
  data,
  ps_col = "ps",
  treatment,
  weight_type = c("ATE", "ATT", "ATC"),
  stabilized = TRUE,
  truncate = c(0.01, 0.99)
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `weight_type` | character | 权重类型: ATE/ATT/ATC |
| `stabilized` | logical | 是否使用稳定化权重 |
| `truncate` | numeric(2) | 权重截断分位数 |

**返回值**: data.table (新增 `weight` 列)

**权重计算公式**:

```r
# ATE权重 (Average Treatment Effect)
w = treatment / ps + (1 - treatment) / (1 - ps)

# ATT权重 (Average Treatment Effect on Treated)
w = treatment + (1 - treatment) * ps / (1 - ps)

# ATC权重 (Average Treatment Effect on Controls)
w = treatment * (1 - ps) / ps + (1 - treatment)

# 稳定化权重
w_stabilized = w * P(treatment) / ps
```

---

#### `assess_balance()`

评估匹配/加权后的平衡性。

**函数签名**:

```r
assess_balance(
  data,
  treatment,
  covariates,
  method = c("unmatched", "matched", "weighted"),
  weight_col = NULL,
  threshold = 0.1
)
```

**返回值**: data.frame

| 列名 | 说明 |
|------|------|
| `variable` | 变量名 |
| `mean_treated` | 处理组均值 |
| `mean_control` | 对照组均值 |
| `smd` | 标准化均值差 |
| `variance_ratio` | 方差比 |
| `balanced` | 是否平衡 (SMD < threshold) |

**平衡性指标**:

| 指标 | 公式 | 目标值 |
|------|------|--------|
| SMD | `(mean_t - mean_c) / sqrt((var_t + var_c) / 2)` | < 0.1 |
| 方差比 | `var_t / var_c` | 0.5 - 2.0 |

---

#### `run_weighted_analysis()`

运行加权回归分析。

**函数签名**:

```r
run_weighted_analysis(
  data,
  exposure,
  outcome,
  covariates = NULL,
  weight_col = "weight",
  model_type = c("cox", "logistic", "linear"),
  endpoint = NULL,
  robust_se = TRUE
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `weight_col` | character | 权重列名 |
| `robust_se` | logical | 是否使用稳健标准误 |

---

### 2.2 技术实现要点

#### 匹配算法

| 方法 | 说明 | 依赖包 |
|------|------|--------|
| Nearest neighbor | 最近邻匹配 | MatchIt |
| Optimal matching | 最优匹配 | optmatch |

- 支持卡尺约束 (caliper)
- 支持精确匹配 (exact matching)
- 支持有/无放回匹配

#### 稳健标准误

- 使用 `sandwich::vcovHC()` 估计器
- 聚类稳健SE (cluster-robust) 可选

---

## 模块三：中介分析 (mediation.R)

### 3.1 概述

中介分析模块基于 **regmedint** 包实现，支持 Valeri & VanderWeele (2013, 2015) 提出的回归基因果中介分析方法。该模块封装了 regmedint 的核心功能，提供更友好的接口和与 UKBAnalytica 其他模块的集成。

**支持的模型组合**:

| 结局模型 (yreg) | 中介模型 (mreg) |
|-----------------|-----------------|
| linear | linear, logistic |
| logistic | linear, logistic |
| poisson | linear, logistic |
| negbin | linear, logistic |
| survCox | linear, logistic |
| survAFT_exp | linear, logistic |
| survAFT_weibull | linear, logistic |

### 3.2 核心函数

#### `run_mediation()`

主函数：运行因果中介分析。

**函数签名**:

```r
run_mediation(
  data,
  exposure,
  mediator,
  outcome,
  covariates = NULL,
  exposure_levels = c(0, 1),
  mediator_value = 0,
  covariate_values = NULL,
  mediator_type = c("continuous", "binary"),
  outcome_type = c("linear", "logistic", "cox"),
  endpoint = NULL,
  interaction = TRUE,
  boot = FALSE,
  boot_n = 1000,
  conf_level = 0.95
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `data` | data.frame/data.table | 数据集 |
| `exposure` | character | 暴露变量名 (处理变量) |
| `mediator` | character | 中介变量名 |
| `outcome` | character | 结局变量名 (Cox模型为时间变量) |
| `covariates` | character vector | 协变量名称向量 |
| `exposure_levels` | numeric(2) | 暴露水平 c(a0, a1)，a0为参考值，a1为比较值 |
| `mediator_value` | numeric | CDE评估时的中介变量固定值 |
| `covariate_values` | numeric vector | 协变量条件值 (用于计算条件效应) |
| `mediator_type` | character | 中介变量类型: "continuous" 或 "binary" |
| `outcome_type` | character | 结局模型类型: "linear", "logistic", "cox" |
| `endpoint` | character(2) | Cox模型需要 c("time_col", "status_col") |
| `interaction` | logical | 是否包含暴露-中介交互项 |
| `boot` | logical | 是否进行Bootstrap置信区间估计 |
| `boot_n` | integer | Bootstrap重抽样次数 |
| `conf_level` | numeric | 置信水平 (默认0.95) |

**返回值**: `mediation_result` 对象 (S3类)，包含:

| 组件 | 说明 |
|------|------|
| `effects` | data.frame，包含各效应的估计值、SE、CI、P值 |
| `mediator_model` | 中介模型对象 |
| `outcome_model` | 结局模型对象 |
| `regmedint_obj` | 原始regmedint对象 |
| `call` | 函数调用 |

**效应指标说明**:

| 效应 | 全称 | 解释 |
|------|------|------|
| `cde` | Controlled Direct Effect | 控制直接效应：将中介变量固定在某值时，暴露的效应 |
| `pnde` | Pure Natural Direct Effect | 纯自然直接效应：不通过中介的暴露效应 (≈传统NDE) |
| `tnde` | Total Natural Direct Effect | 总自然直接效应 |
| `pnie` | Pure Natural Indirect Effect | 纯自然间接效应 |
| `tnie` | Total Natural Indirect Effect | 总自然间接效应：通过中介的暴露效应 (≈传统NIE) |
| `te` | Total Effect | 总效应 = NDE + NIE |
| `pm` | Proportion Mediated | 中介比例 = NIE / TE |

---

#### `run_multi_mediator()`

多中介变量分析：逐一检验多个潜在中介。

**函数签名**:

```r
run_multi_mediator(
  data,
  exposure,
  mediators,
  outcome,
  covariates = NULL,
  mediator_type = "continuous",
  outcome_type = "linear",
  endpoint = NULL,
  ...
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `mediators` | character vector | 多个中介变量名 |
| `...` | - | 传递给 `run_mediation()` 的其他参数 |

**返回值**: data.frame，每行为一个中介变量的分析结果

---

#### `run_sensitivity_mediation()`

中介分析敏感性分析：评估未观测混杂的影响。

**函数签名**:

```r
run_sensitivity_mediation(
  mediation_result,
  rho_values = seq(-0.9, 0.9, by = 0.1)
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `mediation_result` | mediation_result | `run_mediation()` 的返回对象 |
| `rho_values` | numeric vector | 敏感性参数值序列 (未观测混杂的相关系数) |

**返回值**: data.frame，包含不同rho值下的效应估计

---

### 3.3 技术实现要点

#### 模型类型映射

| UKBAnalytica参数 | regmedint参数 |
|------------------|---------------|
| `mediator_type = "continuous"` | `mreg = "linear"` |
| `mediator_type = "binary"` | `mreg = "logistic"` |
| `outcome_type = "linear"` | `yreg = "linear"` |
| `outcome_type = "logistic"` | `yreg = "logistic"` |
| `outcome_type = "cox"` | `yreg = "survCox"` |

#### Bootstrap置信区间

当 `boot = TRUE` 时，使用非参数Bootstrap估计置信区间：

```r
# 伪代码
boot_results <- replicate(boot_n, {
  boot_data <- data[sample(nrow(data), replace = TRUE), ]
  fit_mediation(boot_data)
})
ci <- quantile(boot_results, c((1-conf_level)/2, (1+conf_level)/2))
```

#### 协变量条件值

- 若 `covariate_values = NULL`，默认使用协变量的均值 (连续) 或众数 (分类)
- 用户可指定具体值以计算条件效应

#### 罕见事件假设

对于 logistic 和 Cox 结局模型，效应估计依赖于罕见事件假设。当事件率较高时，结果可能有偏差。建议：
- 事件率 < 10% 时结果可靠
- 事件率 > 10% 时考虑使用 log-binomial 模型或报告警告

---

### 3.4 辅助函数

#### `summary.mediation_result()`

打印中介分析结果摘要。

```r
summary(mediation_obj, exponentiate = FALSE)
```

#### `coef.mediation_result()`
提取效应估计。

```r
coef(mediation_obj)
```

#### `confint.mediation_result()`

提取置信区间。

```r
confint(mediation_obj, level = 0.95)
```

---

## 模块四：可视化 (visualization.R)

### 3.1 核心函数

#### `plot_forest()`

森林图绘制。

**函数签名**:

```r
plot_forest(
  results,
  estimate_col = "HR",
  lower_col = "lower95",
  upper_col = "upper95",
  label_col = "subgroup",
  pvalue_col = "pvalue",
  p_interaction_col = NULL,
  null_value = 1,
  log_scale = TRUE,
  colors = NULL,
  title = "Subgroup Analysis",
  xlab = "Hazard Ratio (95% CI)",
  show_n = TRUE,
  show_events = TRUE
)
```

**返回值**: ggplot2 对象

**森林图布局示例**:

```
|亚组标签 |  N   |事件|    图形区域     | HR(95%CI)    | P值   |P-int|
|---------|------|----|-----------------| -------------|-------|-----|
| 整体    | 1000 | 150|    ●────        | 1.5(1.2-1.8) | 0.001 |  -  |
| 年龄    |      |    |                 |              |       | 0.12|
|   <65   | 500  |  60|   ●───          | 1.3(1.0-1.7) | 0.03  |     |
|   ≥65   | 500  |  90|     ●─────      | 1.7(1.3-2.2) |<0.001 |     |
| 性别    |      |    |                 |              |       | 0.45|
|   男    | 600  | 100|    ●────        | 1.5(1.1-2.0) | 0.008 |     |
|   女    | 400  |  50|   ●──           | 1.4(1.0-1.9) | 0.04  |     |
```

---

#### `plot_km_curve()`

Kaplan-Meier 生存曲线。

**函数签名**:

```r
plot_km_curve(
  data,
  time_col,
  status_col,
  group_col = NULL,
  conf_int = TRUE,
  risk_table = TRUE,
  censor_marks = TRUE,
  palette = "jco",
  title = NULL,
  xlab = "Time (years)",
  ylab = "Survival Probability",
  legend_title = "Group",
  median_line = TRUE,
  pvalue = TRUE,
  xlim = NULL,
  break_time = NULL
)
```

**K-M曲线组件**:

- 生存曲线 (按组着色)
- 置信区间带
- 删失标记 (censoring marks)
- 风险表 (number at risk)
- Log-rank P值
- 中位生存时间参考线

---

#### `plot_ps_distribution()`

倾向性评分分布图。

**函数签名**:

```r
plot_ps_distribution(
  data,
  ps_col = "ps",
  treatment,
  type = c("histogram", "density", "mirror"),
  matched = FALSE,
  match_col = NULL
)
```

**图形类型**:

| 类型 | 说明 |
|------|------|
| histogram | 直方图 |
| density | 密度曲线 |
| mirror | 镜像直方图 (处理/对照组分上下) |

---

#### `plot_balance()`

平衡性诊断图 (Love plot)。

**函数签名**:

```r
plot_balance(
  balance_before,
  balance_after,
  threshold = 0.1,
  title = "Covariate Balance",
  xlab = "Standardized Mean Difference"
)
```

**图形特点**:

- X轴: 标准化均值差 (SMD)
- Y轴: 协变量名称
- 点: 匹配前 (空心) vs 匹配后 (实心)
- 参考线: SMD = ±0.1 阈值

---

#### `plot_calibration()`

校准曲线。

**函数签名**:

```r
plot_calibration(
  data,
  predicted,
  observed,
  n_bins = 10,
  smooth = TRUE,
  conf_int = TRUE
)
```

---

#### `plot_mediation()`

中介分析路径图/效应图。

**函数签名**:

```r
plot_mediation(
  mediation_result,
  type = c("path", "effects", "decomposition"),
  show_ci = TRUE,
  show_pvalue = TRUE,
  exponentiate = FALSE,
  title = NULL,
  colors = NULL
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `mediation_result` | mediation_result | `run_mediation()` 返回对象 |
| `type` | character | 图形类型: "path"(路径图), "effects"(效应条形图), "decomposition"(效应分解饼图) |
| `show_ci` | logical | 是否显示置信区间 |
| `show_pvalue` | logical | 是否显示P值 |
| `exponentiate` | logical | 是否指数化 (用于HR/OR) |
| `title` | character | 图标题 |
| `colors` | character vector | 自定义颜色 |

**图形类型说明**:

1. **路径图 (path)**: 展示暴露→中介→结局的因果路径

```
                    ┌─────────┐
         a路径      │         │      b路径
    ┌──────────────►│ Mediator │──────────────┐
    │               │         │               │
    │               └─────────┘               ▼
┌───┴────┐                               ┌────────┐
│Exposure│───────────────────────────────►│Outcome │
└────────┘         c'路径 (直接效应)       └────────┘
```

2. **效应条形图 (effects)**: 展示各效应估计及CI

3. **效应分解图 (decomposition)**: 展示TE分解为NDE+NIE的比例

---

#### `plot_mediation_forest()`

多中介变量结果的森林图。

**函数签名**:

```r
plot_mediation_forest(
  multi_mediation_result,
  effect_type = c("tnie", "pnde", "te", "pm"),
  exponentiate = FALSE,
  null_value = 0,
  title = "Mediation Analysis: Multiple Mediators"
)
```

---

### 4.2 可视化设计规范

#### 配色方案

默认使用 `ggsci` 包的期刊配色方案:

| 方案 | 说明 |
|------|------|
| `jco` | Journal of Clinical Oncology |
| `nejm` | New England Journal of Medicine |
| `lancet` | The Lancet |
| `npg` | Nature Publishing Group |

支持自定义调色板。

#### 输出格式

所有可视化函数返回 `ggplot2` 对象，可以:

- 进一步自定义样式
- 使用 `ggsave()` 保存为多种格式 (PNG, PDF, SVG等)
- 组合多图 (使用 `patchwork` 或 `cowplot`)

---

## 依赖管理

### 新增 Imports

```
scales        # 颜色/轴缩放工具
```

### 新增 Suggests

```
MatchIt       # PSM匹配
sandwich      # 稳健标准误
lmtest        # 稳健SE提取
regmedint     # 中介分析核心包 (重要!)
survminer     # K-M曲线增强
cobalt        # 平衡诊断
gbm           # 梯度提升PS估计
ggsci         # 期刊配色方案
```

### 更新 DESCRIPTION 示例

```
Imports: 
    data.table,
    stringi,
    ggplot2,
    scales
Suggests:
    testthat,
    survival,
    tableone,
    mice,
    MatchIt,
    sandwich,
    lmtest,
    regmedint,
    survminer,
    cobalt,
    gbm,
    ggsci
```

---

## 函数依赖关系

```
subgroup.R
├── run_subgroup_analysis()
│   ├── .validate_regression_inputs()    ← [from regression.R]
│   ├── .calculate_interaction_pvalue()  ← [internal]
│   └── runmulti_cox/lm/logit()          ← [reuse existing]
└── run_multi_subgroup()
    └── run_subgroup_analysis()

propensity.R
├── estimate_propensity_score()
│   └── stats::glm() / gbm::gbm()
├── match_propensity()
│   └── MatchIt::matchit()
├── calculate_weights()
├── assess_balance()
└── run_weighted_analysis()
    └── sandwich::vcovHC()

mediation.R
├── run_mediation()
│   └── regmedint::regmedint()
├── run_multi_mediator()
│   └── run_mediation()
├── run_sensitivity_mediation()
├── summary.mediation_result()
├── coef.mediation_result()
└── confint.mediation_result()

visualization.R
├── plot_forest()
│   └── ggplot2
├── plot_km_curve()
│   └── survival::survfit() + ggplot2
├── plot_ps_distribution()
│   └── ggplot2
├── plot_balance()
│   └── ggplot2
├── plot_calibration()
│   └── ggplot2
├── plot_mediation()              ← [新增]
│   └── ggplot2
└── plot_mediation_forest()       ← [新增]
    └── ggplot2
```

---

## 使用示例

### 完整工作流示例

```r
library(UKBAnalytica)

# =============================================================================
# 步骤 1: 准备数据
# =============================================================================
surv_data <- build_survival_dataset(
  ukb_data, 
  diseases, 
  primary_disease = "CVD"
)

# =============================================================================
# 步骤 2: 亚组分析
# =============================================================================
subgroup_results <- run_multi_subgroup(
  data = surv_data,
  exposure = "smoking",
  outcome = NULL,
  subgroup_vars = c("sex", "age_group", "bmi_cat"),
  covariates = c("education", "income"),
  model_type = "cox",
  endpoint = c("outcome_surv_time", "outcome_status")
)

# 查看结果
print(subgroup_results)

# =============================================================================
# 步骤 3: 森林图可视化
# =============================================================================
forest_plot <- plot_forest(
  subgroup_results, 
  title = "Smoking and CVD Risk by Subgroups"
)
print(forest_plot)

# 保存图片
ggsave("forest_plot.png", forest_plot, width = 10, height = 8, dpi = 300)

# =============================================================================
# 步骤 4: 倾向性评分估计
# =============================================================================
ps_data <- estimate_propensity_score(
  data = surv_data,
  treatment = "smoking_binary",
  covariates = c("age", "sex", "bmi", "education", "income")
)

# 检查PS分布
ps_dist_plot <- plot_ps_distribution(
  ps_data, 
  treatment = "smoking_binary",
  type = "mirror"
)
print(ps_dist_plot)

# =============================================================================
# 步骤 5: 倾向性评分匹配
# =============================================================================
matched_data <- match_propensity(
  data = ps_data,
  treatment = "smoking_binary",
  ratio = 1,
  caliper = 0.2,
  method = "nearest"
)

cat("匹配前样本量:", nrow(ps_data), "\n")
cat("匹配后样本量:", nrow(matched_data), "\n")

# =============================================================================
# 步骤 6: 平衡性诊断
# =============================================================================
covars <- c("age", "sex", "bmi", "education", "income")

balance_before <- assess_balance(
  ps_data, 
  treatment = "smoking_binary",
  covariates = covars,
  method = "unmatched"
)

balance_after <- assess_balance(
  matched_data, 
  treatment = "smoking_binary",
  covariates = covars,
  method = "matched"
)

# Love plot
balance_plot <- plot_balance(balance_before, balance_after)
print(balance_plot)

# =============================================================================
# 步骤 7: K-M生存曲线
# =============================================================================
km_plot <- plot_km_curve(
  matched_data, 
  time_col = "outcome_surv_time",
  status_col = "outcome_status",
  group_col = "smoking_binary",
  risk_table = TRUE,
  pvalue = TRUE,
  palette = "jco"
)
print(km_plot)

# =============================================================================
# 步骤 8: 匹配后效应估计
# =============================================================================
matched_results <- runmulti_cox(
  data = matched_data,
  main_var = "smoking_binary",
  covariates = covars,
  endpoint = c("outcome_surv_time", "outcome_status")
)

print(matched_results)

# =============================================================================
# 可选: IPTW加权分析
# =============================================================================
weighted_data <- calculate_weights(
  data = ps_data,
  treatment = "smoking_binary",
  weight_type = "ATE",
  stabilized = TRUE
)

weighted_results <- run_weighted_analysis(
  data = weighted_data,
  exposure = "smoking_binary",
  outcome = NULL,
  covariates = covars,
  weight_col = "weight",
  model_type = "cox",
  endpoint = c("outcome_surv_time", "outcome_status"),
  robust_se = TRUE
)

print(weighted_results)
```

---

### 中介分析工作流示例

```r
library(UKBAnalytica)

# =============================================================================
# 示例: 探索BMI是否中介吸烟与CVD的关系
# =============================================================================

# 准备数据
analysis_data <- surv_data[!is.na(bmi) & !is.na(smoking_binary), ]

# =============================================================================
# 步骤 1: 单中介变量分析 (连续型中介)
# =============================================================================
mediation_bmi <- run_mediation(
  data = analysis_data,
  exposure = "smoking_binary",
  mediator = "bmi",
  outcome = "outcome_surv_time",
  covariates = c("age", "sex", "education"),
  exposure_levels = c(0, 1),       # 非吸烟 vs 吸烟
  mediator_value = 25,             # CDE评估时BMI固定在25
  covariate_values = c(55, 1, 2),  # 年龄55, 男性, 中等教育
  mediator_type = "continuous",
  outcome_type = "cox",
  endpoint = c("outcome_surv_time", "outcome_status"),
  interaction = TRUE
)

# 查看结果
summary(mediation_bmi)

# 提取效应估计
coef(mediation_bmi)
#>       effect       est        se      lower     upper        p
#> 1       cde  0.4521    0.1832    0.0931    0.8112   0.0134
#> 2      pnde  0.4102    0.1654    0.0860    0.7344   0.0131
#> 3      tnie  0.0523    0.0287    0.0040    0.1006   0.0341
#> 4      tnde  0.4156    0.1672    0.0879    0.7433   0.0128
#> 5      pnie  0.0469    0.0265   -0.0050    0.0988   0.0764
#> 6        te  0.4625    0.1698    0.1297    0.7953   0.0065
#> 7        pm  0.1131    0.0612    0.0086    0.2176   0.0341

# 中介比例约11.3%，P=0.034

# =============================================================================
# 步骤 2: 可视化中介效应
# =============================================================================

# 路径图
path_plot <- plot_mediation(mediation_bmi, type = "path")
print(path_plot)

# 效应分解图
decomp_plot <- plot_mediation(mediation_bmi, type = "decomposition")
print(decomp_plot)

# =============================================================================
# 步骤 3: 多中介变量分析
# =============================================================================
multi_med_results <- run_multi_mediator(
  data = analysis_data,
  exposure = "smoking_binary",
  mediators = c("bmi", "blood_pressure", "ldl", "hba1c"),
  outcome = "outcome_surv_time",
  covariates = c("age", "sex", "education"),
  mediator_type = "continuous",
  outcome_type = "cox",
  endpoint = c("outcome_surv_time", "outcome_status")
)

print(multi_med_results)
#>      mediator    tnie     tnie_se    tnie_p      pm     pm_se
#> 1         bmi  0.0523    0.0287     0.0341   0.113    0.061
#> 2 blood_press  0.0312    0.0198     0.1152   0.067    0.043
#> 3         ldl  0.0089    0.0145     0.5392   0.019    0.031
#> 4       hba1c  0.0156    0.0112     0.1632   0.034    0.024

# 森林图展示多中介变量结果
forest_med <- plot_mediation_forest(
  multi_med_results,
  effect_type = "tnie",
  title = "Indirect Effects via Different Mediators"
)
print(forest_med)

# =============================================================================
# 步骤 4: 二分类中介变量分析 (糖尿病作为中介)
# =============================================================================
mediation_dm <- run_mediation(
  data = analysis_data,
  exposure = "smoking_binary",
  mediator = "diabetes_history",  # 0/1 二分类
  outcome = "outcome_surv_time",
  covariates = c("age", "sex", "education"),
  exposure_levels = c(0, 1),
  mediator_value = 0,             # CDE评估时假设无糖尿病
  mediator_type = "binary",       # 重要: 指定为二分类
  outcome_type = "cox",
  endpoint = c("outcome_surv_time", "outcome_status"),
  interaction = TRUE
)

summary(mediation_dm)

# =============================================================================
# 步骤 5: Bootstrap置信区间 (更稳健)
# =============================================================================
mediation_boot <- run_mediation(
  data = analysis_data,
  exposure = "smoking_binary",
  mediator = "bmi",
  outcome = "outcome_surv_time",
  covariates = c("age", "sex"),
  mediator_type = "continuous",
  outcome_type = "cox",
  endpoint = c("outcome_surv_time", "outcome_status"),
  boot = TRUE,          # 启用Bootstrap
  boot_n = 1000,        # 1000次重抽样
  conf_level = 0.95
)

# Bootstrap CI 更适合小样本或非正态分布
confint(mediation_boot)

# =============================================================================
# 步骤 6: 线性结局模型 (连续型结局)
# =============================================================================
mediation_linear <- run_mediation(
  data = analysis_data,
  exposure = "smoking_binary",
  mediator = "bmi",
  outcome = "systolic_bp",        # 连续型结局
  covariates = c("age", "sex"),
  mediator_type = "continuous",
  outcome_type = "linear",        # 线性模型
  interaction = TRUE
)

summary(mediation_linear)

# =============================================================================
# 步骤 7: Logistic结局模型 (二分类结局)
# =============================================================================
mediation_logit <- run_mediation(
  data = analysis_data,
  exposure = "smoking_binary",
  mediator = "bmi",
  outcome = "hypertension",       # 0/1 二分类结局
  covariates = c("age", "sex"),
  mediator_type = "continuous",
  outcome_type = "logistic",      # Logistic模型
  interaction = TRUE
)

# 指数化以获得OR
summary(mediation_logit, exponentiate = TRUE)
```

---

## 开发计划

| 阶段 | 内容 | 预计输出 |
|------|------|----------|
| **Phase 1** | `subgroup.R` 核心函数实现 | 亚组分析功能 ✅ |
| **Phase 2** | `propensity.R` 核心函数实现 | PSM/IPTW功能 ✅ |
| **Phase 3** | `visualization.R` 可视化函数 | 森林图、K-M曲线等 ✅ |
| **Phase 4** | `mediation.R` 中介分析实现 | 因果中介分析功能 |
| **Phase 5** | 单元测试与文档 | testthat测试、roxygen2文档 |
| **Phase 6** | 集成测试与示例 | vignette、使用示例 |

---

## 附录

### A. 参考文献

**倾向性评分:**

1. Austin PC. An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. *Multivariate Behav Res*. 2011;46(3):399-424.

2. Rosenbaum PR, Rubin DB. The central role of the propensity score in observational studies for causal effects. *Biometrika*. 1983;70(1):41-55.

3. Wang X, et al. Propensity score methods in health technology assessment: principles, extended applications, and recent advances. *Front Pharmacol*. 2022.

**因果中介分析:**

4. Valeri L, VanderWeele TJ. Mediation analysis allowing for exposure-mediator interactions and causal interpretation: theoretical assumptions and implementation with SAS and SPSS macros. *Psychological Methods*. 2013;18(2):137-150.

5. Valeri L, VanderWeele TJ. SAS macro for causal mediation analysis with survival data. *Epidemiology*. 2015;26(2):e23-e24.

6. VanderWeele TJ. *Explanation in Causal Inference: Methods for Mediation and Interaction*. Oxford University Press. 2015.

7. Yoshida K, Li S, Mathur MB. regmedint: Regression-Based Causal Mediation Analysis with Interaction and Effect Modification Terms. R package. 2021. https://kaz-yos.github.io/regmedint/

### B. 相关R包

| 包名 | 功能 |
|------|------|
| MatchIt | 倾向性评分匹配 |
| cobalt | 平衡性诊断 |
| WeightIt | 权重估计 |
| survminer | 生存曲线可视化 |
| forestplot | 森林图 |
| **regmedint** | **因果中介分析 (本模块核心依赖)** |
| mediation | 传统中介分析 (Baron & Kenny方法) |
| medflex | 自然效应模型中介分析 |

### C. 中介效应术语对照表

| 英文缩写 | 英文全称 | 中文翻译 |
|----------|----------|----------|
| CDE | Controlled Direct Effect | 控制直接效应 |
| NDE | Natural Direct Effect | 自然直接效应 |
| NIE | Natural Indirect Effect | 自然间接效应 |
| PNDE | Pure Natural Direct Effect | 纯自然直接效应 |
| TNDE | Total Natural Direct Effect | 总自然直接效应 |
| PNIE | Pure Natural Indirect Effect | 纯自然间接效应 |
| TNIE | Total Natural Indirect Effect | 总自然间接效应 |
| TE | Total Effect | 总效应 |
| PM | Proportion Mediated | 中介比例 |

### D. 中介分析因果假设

使用本模块进行因果中介分析时，需满足以下假设：

1. **无未测量的暴露-结局混杂** (given covariates)
2. **无未测量的中介-结局混杂** (given covariates and exposure)
3. **无未测量的暴露-中介混杂** (given covariates)
4. **无中介-结局混杂受暴露影响** (cross-world independence)

若假设4不成立，只有CDE和TE可识别；NDE和NIE不可识别。

---

*文档结束*

---

## 模块五：多重插补结果合并 (mi_pool.R)

### 5.1 概述

多重插补 (Multiple Imputation, MI) 是处理缺失数据的标准方法。在完成多重插补后，需要在每个插补数据集上分别运行分析，然后使用 Rubin's Rules 合并结果。本模块提供了便捷的接口，支持多种回归模型的结果合并。

**核心依赖**: `mitools` 包

**支持的模型类型**:
| 模型类型 | R函数 | 估计量 |
|---------|-------|--------|
| 线性回归 | `lm()` | Beta系数 |
| 逻辑回归 | `glm(family=binomial)` | Odds Ratio |
| Poisson回归 | `glm(family=poisson)` | Rate Ratio |
| Cox回归 | `survival::coxph()` | Hazard Ratio |
| 负二项回归 | `MASS::glm.nb()` | Rate Ratio |

---

### 5.2 核心函数

#### `pool_mi_models()`

主函数：合并多个插补数据集上的回归模型结果。

**函数签名**:

```r
pool_mi_models(
  models = NULL,
  datasets = NULL,
  formula = NULL,
  model_type = c("lm", "logistic", "poisson", "cox", "negbin"),
  family = NULL,
  df.complete = Inf,
  conf.level = 0.95,
  exponentiate = NULL
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `models` | list | 已拟合模型的列表 (每个插补数据集一个模型) |
| `datasets` | list/imputationList | 插补数据集列表 (当models=NULL时使用) |
| `formula` | formula | 回归公式 (当datasets提供时需要) |
| `model_type` | character | 模型类型: "lm", "logistic", "poisson", "cox", "negbin" |
| `family` | family | glm的family参数 (可选，自动根据model_type推断) |
| `df.complete` | numeric | 完整数据自由度 (用于小样本校正) |
| `conf.level` | numeric | 置信区间水平，默认0.95 |
| `exponentiate` | logical | 是否对系数取指数 (NULL=自动判断) |

**返回值**: `mi_pooled_result` 对象，包含:

| 组件 | 说明 |
|------|------|
| `pooled` | 合并后的系数表 (data.frame) |
| `mi_result` | mitools::MIcombine 返回的原始对象 |
| `n_imputations` | 插补数据集数量 |
| `model_type` | 模型类型 |
| `formula` | 使用的公式 |
| `call` | 函数调用 |

**pooled data.frame 列**:

| 列名 | 说明 |
|------|------|
| `term` | 变量名 |
| `estimate` | 合并估计值 (或exp后的HR/OR) |
| `std.error` | 合并标准误 |
| `statistic` | Wald统计量 |
| `p.value` | P值 |
| `conf.low` | 置信区间下限 |
| `conf.high` | 置信区间上限 |
| `fmi` | 缺失信息比例 (Fraction of Missing Information) |
| `df` | Rubin's 调整自由度 |

---

#### `fit_mi_models()`

辅助函数：在多个插补数据集上拟合回归模型。

**函数签名**:

```r
fit_mi_models(
  datasets,
  formula,
  model_type = c("lm", "logistic", "poisson", "cox", "negbin"),
  family = NULL,
  ...
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `datasets` | list/imputationList | 插补数据集列表 |
| `formula` | formula | 回归公式 |
| `model_type` | character | 模型类型 |
| `family` | family | glm的family参数 |
| `...` | | 传递给模型函数的额外参数 |

**返回值**: 模型对象列表

---

#### `create_imputation_list()`

辅助函数：将数据集列表转换为 `imputationList` 对象。

**函数签名**:

```r
create_imputation_list(
  datasets,
  validate = TRUE
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `datasets` | list | data.frame列表 |
| `validate` | logical | 是否验证数据结构一致性 |

**返回值**: `imputationList` 对象

---

#### `pool_custom_estimates()`

高级函数：合并自定义估计量（不限于回归系数）。

**函数签名**:

```r
pool_custom_estimates(
  estimates,
  variances,
  df.complete = Inf,
  conf.level = 0.95,
  labels = NULL
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `estimates` | list | 估计量向量的列表 |
| `variances` | list | 方差-协方差矩阵的列表 |
| `df.complete` | numeric | 完整数据自由度 |
| `conf.level` | numeric | 置信区间水平 |
| `labels` | character | 估计量标签 |

**返回值**: `mi_pooled_result` 对象

---

### 5.3 S3方法

```r
# 打印方法
print.mi_pooled_result(x, digits = 3, ...)

# 摘要方法
summary.mi_pooled_result(object, exponentiate = NULL, ...)

# 提取系数
coef.mi_pooled_result(object, ...)

# 提取置信区间
confint.mi_pooled_result(object, parm, level = 0.95, ...)

# 提取方差-协方差矩阵
vcov.mi_pooled_result(object, ...)

# 整洁输出 (兼容broom包风格)
tidy.mi_pooled_result(x, conf.int = TRUE, conf.level = 0.95, 
                       exponentiate = FALSE, ...)
```

---

### 5.4 可视化函数

#### `plot_mi_pooled()`

绘制多重插补合并结果的森林图。

**函数签名**:

```r
plot_mi_pooled(
  mi_result,
  terms = NULL,
  exponentiate = NULL,
  null_value = NULL,
  title = "Pooled Estimates (Multiple Imputation)",
  colors = NULL,
  show_fmi = TRUE
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `mi_result` | mi_pooled_result | pool_mi_models的返回对象 |
| `terms` | character | 要显示的变量 (NULL=全部) |
| `exponentiate` | logical | 是否取指数 |
| `null_value` | numeric | 无效值参考线 (0或1) |
| `title` | character | 图表标题 |
| `colors` | character | 颜色向量 |
| `show_fmi` | logical | 是否显示FMI |

**返回值**: ggplot2对象

---

#### `plot_mi_diagnostics()`

绘制多重插补诊断图。

**函数签名**:

```r
plot_mi_diagnostics(
  mi_result,
  type = c("fmi", "variance_ratio", "df"),
  title = NULL
)
```

**参数说明**:

| 参数 | 类型 | 说明 |
|------|------|------|
| `mi_result` | mi_pooled_result | pool_mi_models的返回对象 |
| `type` | character | 诊断图类型 |
| `title` | character | 图表标题 |

**返回值**: ggplot2对象

---

### 5.5 与现有模块的集成

本模块设计为独立模块，不修改现有函数。用户可以通过以下方式集成：

#### 与亚组分析集成

```r
# 在每个插补数据集上运行亚组分析
subgroup_results <- lapply(mi_datasets, function(d) {
  run_subgroup_analysis(
    data = d,
    exposure = "treatment",
    outcome = "event",
    subgroup_var = "age_group",
    model_type = "cox",
    endpoint = c("time", "status")
  )
})

# 手动合并特定亚组的估计
```

#### 与倾向性评分集成

```r
# 在每个插补数据集上计算PS并运行加权分析
ps_results <- lapply(mi_datasets, function(d) {
  ps <- estimate_propensity_score(d, "treatment", covariates)
  d$weight <- calculate_weights(ps$ps, d$treatment, "ate")
  run_weighted_analysis(d, "treatment", "outcome", "weight", "logistic")
})

# 合并加权分析结果
pooled <- pool_mi_models(ps_results)
```

#### 与中介分析集成

```r
# regmedint 本身支持 bootstrap，可在每个插补数据集上运行
med_results <- lapply(mi_datasets, function(d) {
  run_mediation(d, "exposure", "mediator", "outcome",
                covariates = covs, outcome_type = "cox")
})

# 合并中介效应 (使用 pool_custom_estimates)
nie_estimates <- lapply(med_results, function(r) r$effects$est[r$effects$effect == "tnie"])
nie_se <- lapply(med_results, function(r) r$effects$se[r$effects$effect == "tnie"])
pooled_nie <- pool_custom_estimates(
  estimates = nie_estimates,
  variances = lapply(nie_se, function(x) matrix(x^2)),
  labels = "TNIE"
)
```

---

### 5.6 使用示例

#### 基础用法：合并已拟合模型

```r
library(UKBAnalytica)
library(mice)
library(survival)

# 假设已有多重插补数据
# mi_data <- mice(original_data, m = 5)
# mi_datasets <- complete(mi_data, action = "all")

# 方法1: 先拟合模型，再合并
models <- lapply(mi_datasets, function(d) {
  coxph(Surv(time, status) ~ treatment + age + sex, data = d)
})
pooled <- pool_mi_models(models, model_type = "cox")
summary(pooled)

# 方法2: 一步完成拟合和合并
pooled <- pool_mi_models(
  datasets = mi_datasets,
  formula = Surv(time, status) ~ treatment + age + sex,
  model_type = "cox"
)
print(pooled)
```

#### 使用 imputationList

```r
# 创建 imputationList
mi_list <- create_imputation_list(mi_datasets)

# 使用 with() 语法
models <- with(mi_list, coxph(Surv(time, status) ~ treatment + age + sex))

# 合并结果
pooled <- pool_mi_models(models, model_type = "cox")
```

#### 逻辑回归示例

```r
# 二元结局的逻辑回归
pooled_logit <- pool_mi_models(
  datasets = mi_datasets,
  formula = outcome ~ exposure + age + sex + bmi,
  model_type = "logistic",
  exponentiate = TRUE  # 输出OR
)

# 查看合并结果
summary(pooled_logit)
#>              term estimate std.error statistic  p.value conf.low conf.high    fmi
#> 1     (Intercept)    0.123     0.234    -2.345    0.019    0.078     0.194  0.052
#> 2        exposure    1.456     0.156     2.345    0.019    1.072     1.978  0.078
#> 3             age    1.023     0.008     2.875    0.004    1.007     1.039  0.045
#> ...
```

#### 可视化

```r
# 森林图
plot_mi_pooled(pooled, exponentiate = TRUE)

# 诊断图：缺失信息比例
plot_mi_diagnostics(pooled, type = "fmi")
```

#### 自定义估计量合并

```r
# 假设要合并每个插补数据集的均值差
mean_diffs <- lapply(mi_datasets, function(d) {
  mean(d$outcome[d$treatment == 1]) - mean(d$outcome[d$treatment == 0])
})

# 估计方差 (简化示例)
vars <- lapply(mi_datasets, function(d) {
  n1 <- sum(d$treatment == 1)
  n0 <- sum(d$treatment == 0)
  var(d$outcome[d$treatment == 1])/n1 + var(d$outcome[d$treatment == 0])/n0
})

pooled_diff <- pool_custom_estimates(
  estimates = mean_diffs,
  variances = lapply(vars, function(v) matrix(v)),
  labels = "Mean Difference"
)
summary(pooled_diff)
```

---

### 5.7 Rubin's Rules 原理

多重插补合并使用 Rubin's Rules：

**合并点估计**:
$$\bar{Q} = \frac{1}{m} \sum_{i=1}^{m} \hat{Q}_i$$

**合并方差** (总方差 = 组内方差 + 组间方差):
$$T = \bar{U} + (1 + \frac{1}{m}) B$$

其中：
- $\bar{U} = \frac{1}{m} \sum_{i=1}^{m} U_i$ (组内方差的平均)
- $B = \frac{1}{m-1} \sum_{i=1}^{m} (\hat{Q}_i - \bar{Q})^2$ (组间方差)

**缺失信息比例 (FMI)**:
$$\lambda = \frac{(1 + m^{-1}) B}{T}$$

**调整自由度** (Barnard-Rubin):
$$\nu = (m-1)(1 + \frac{1}{\lambda (1 + m^{-1})})^2$$

---

### 5.8 注意事项

1. **插补数量**: 建议 m ≥ 20 个插补数据集，以获得稳定的估计
2. **模型一致性**: 所有插补数据集上的模型规格必须相同
3. **收敛检查**: 确保所有模型都成功收敛
4. **Cox模型**: 需确保 `survival` 包已安装
5. **FMI解释**: FMI > 0.5 表示缺失数据对估计有较大影响

---

### 5.9 依赖包

| 包 | 类型 | 用途 |
|----|------|------|
| mitools | Suggests | MI结果合并核心功能 |
| survival | Suggests | Cox回归模型 |
| MASS | Suggests | 负二项回归 |

---
