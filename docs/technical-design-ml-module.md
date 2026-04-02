# UKBAnalytica 机器学习模块技术设计文档

## 1. 概述

### 1.1 模块目标

为 UK Biobank 数据分析提供统一的机器学习接口，支持：
- 常见分类/回归模型的训练与预测
- 生存分析机器学习模型
- SHAP 可解释性分析
- 模型评估与比较
- 与现有包模块的无缝集成

### 1.2 设计原则

1. **统一接口**: 所有模型使用相同的 API 风格
2. **惰性依赖**: 机器学习包作为 Suggests，按需加载
3. **UKB 优化**: 针对大规模数据优化，支持采样策略
4. **可解释性优先**: 内置 SHAP 解释，适合医学研究发表
5. **与现有模块集成**: 可与倾向性评分、多重插补等模块配合

---

## 2. 依赖包选择

### 2.1 核心依赖 (Suggests)

| 包名 | 用途 | 选择理由 |
|------|------|----------|
| `ranger` | 随机森林 | 快速、内存高效、支持生存分析 |
| `xgboost` | 梯度提升树 | 高性能、广泛使用 |
| `glmnet` | 正则化回归 | LASSO/Ridge/Elastic Net |
| `e1071` | SVM | 经典支持向量机 |
| `nnet` | 神经网络 | 单层神经网络，轻量 |
| `fastshap` | SHAP 值 | 快速 SHAP 计算 |
| `shapviz` | SHAP 可视化 | 高质量 SHAP 图表 |
| `pROC` | ROC 分析 | AUC 计算与置信区间 |
| `caret` | 模型框架 | 交叉验证、预处理 |
| `mlr3` | 模型框架 | 现代 ML 框架（可选） |

### 2.2 生存分析 ML

| 包名 | 用途 |
|------|------|
| `randomForestSRC` | 随机生存森林 |
| `gbm` | 梯度提升生存模型 |
| `survival` | Cox 模型基础 |

---

## 3. API 设计

### 3.1 核心函数概览

```
ml_module.R
├── ukb_ml_model()           # 统一模型训练接口
├── ukb_ml_predict()         # 统一预测接口
├── ukb_ml_cv()              # 交叉验证
├── ukb_ml_tune()            # 超参数调优
├── ukb_ml_compare()         # 模型比较
├── ukb_ml_importance()      # 变量重要性
│
ml_shap.R
├── ukb_shap()               # SHAP 值计算
├── ukb_shap_summary()       # SHAP 汇总图
├── ukb_shap_dependence()    # SHAP 依赖图
├── ukb_shap_force()         # SHAP 力图（单样本）
├── ukb_shap_interaction()   # SHAP 交互分析
│
ml_survival.R
├── ukb_ml_survival()        # 生存分析 ML 模型
├── ukb_ml_survival_predict()# 生存预测
├── ukb_ml_survival_shap()   # 生存模型 SHAP
│
ml_evaluate.R
├── ukb_ml_metrics()         # 评估指标
├── ukb_ml_calibration()     # 校准曲线
├── ukb_ml_roc()             # ROC 曲线
├── ukb_ml_confusion()       # 混淆矩阵
```

### 3.2 详细 API 设计

#### 3.2.1 ukb_ml_model() - 统一模型训练

```r
ukb_ml_model(
  formula,
  data,
  model = c("rf", "xgboost", "glmnet", "svm", "nnet", "logistic"),
  task = c("classification", "regression"),
  

  # 数据处理
  split_ratio = 0.8,
  stratify = TRUE,
  seed = NULL,
  sample_n = NULL,              # 大数据采样
  
  # 模型参数
  params = list(),              # 模型特定参数
  
  # 交叉验证
  cv = FALSE,
  cv_folds = 5,
  
  # 控制
  verbose = TRUE,
  ...
)
```

**返回值**: `ukb_ml` S3 对象

```r
# 返回结构
list(
  model = <fitted_model>,       # 原始模型对象
  model_type = "rf",            # 模型类型
  task = "classification",      # 任务类型
  formula = formula,            # 公式
  predictors = c(...),          # 预测变量名
  outcome = "...",              # 结局变量名
  train_data = <data>,          # 训练数据

  test_data = <data>,           # 测试数据
  train_idx = c(...),           # 训练集索引
  test_idx = c(...),            # 测试集索引
  metrics = list(...),          # 评估指标
  cv_results = <if cv=TRUE>,    # 交叉验证结果
  call = <call>                 # 调用信息
)
```

**示例**:

```r
# 随机森林分类
ml_rf <- ukb_ml_model(
  diabetes ~ age + bmi + sbp + smoking + alcohol,
  data = ukb_data,
  model = "rf",
  task = "classification",
  split_ratio = 0.8,
  params = list(num.trees = 500, mtry = 3)
)

# XGBoost 回归
ml_xgb <- ukb_ml_model(
  bmi ~ age + sex + income + education,
  data = ukb_data,
  model = "xgboost",
  task = "regression",
  cv = TRUE,
  cv_folds = 5
)

# LASSO 逻辑回归
ml_lasso <- ukb_ml_model(
  cvd ~ .,
  data = ukb_data,
  model = "glmnet",
  task = "classification",
  params = list(alpha = 1)  # LASSO
)
```

#### 3.2.2 ukb_ml_predict() - 统一预测

```r
ukb_ml_predict(
  object,                       # ukb_ml 对象

  newdata = NULL,               # 新数据，NULL 则用测试集
  type = c("response", "prob", "class", "link"),
  ...
)
```

#### 3.2.3 ukb_ml_cv() - 交叉验证

```r
ukb_ml_cv(
  formula,
  data,
  model = "rf",
  task = "classification",
  folds = 5,
  repeats = 1,
  stratify = TRUE,
  metrics = c("auc", "accuracy", "sensitivity", "specificity"),
  parallel = FALSE,
  seed = NULL,
  ...
)
```

**返回值**: `ukb_ml_cv` S3 对象

```r
list(
  cv_metrics = data.frame(...),  # 每折指标
  mean_metrics = c(...),         # 平均指标
  sd_metrics = c(...),           # 标准差
  fold_models = list(...),       # 每折模型（可选保留）
  fold_predictions = list(...)   # 每折预测
)
```

#### 3.2.4 ukb_ml_tune() - 超参数调优

```r
ukb_ml_tune(
  formula,
  data,
  model = "rf",
  task = "classification",
  
  # 参数网格
  param_grid = list(
    num.trees = c(100, 500, 1000),
    mtry = c(2, 4, 6),
    min.node.size = c(1, 5, 10)
  ),
  
  # 搜索策略
  search = c("grid", "random"),
  n_iter = 20,                   # 随机搜索迭代数
  
  # 验证
  cv_folds = 5,
  metric = "auc",                # 优化目标
  
  parallel = FALSE,
  seed = NULL,
  verbose = TRUE
)
```

**返回值**: `ukb_ml_tune` S3 对象

```r
list(
  best_params = list(...),       # 最优参数
  best_score = 0.85,             # 最优得分
  results = data.frame(...),     # 所有参数组合结果
  best_model = <ukb_ml>          # 最优模型
)
```

#### 3.2.5 ukb_ml_compare() - 模型比较

```r
ukb_ml_compare(
  ...,                           # 多个 ukb_ml 对象
  models = list(),               # 或以列表传入
  metrics = c("auc", "accuracy", "brier"),
  test_data = NULL,              # 统一测试集
  plot = TRUE
)
```

**返回值**: `ukb_ml_compare` S3 对象

```r
list(
  comparison = data.frame(
    model = c("rf", "xgboost", "glmnet"),
    auc = c(0.85, 0.87, 0.82),
    accuracy = c(0.80, 0.82, 0.78),
    ...
  ),
  statistical_tests = list(...), # 模型间统计检验
  plot = <ggplot>
)
```

---

### 3.3 SHAP 可解释性 API

#### 3.3.1 ukb_shap() - SHAP 值计算

```r
ukb_shap(

  object,                        # ukb_ml 对象
  data = NULL,                   # 解释数据，NULL 用测试集
  nsim = 100,                    # 蒙特卡洛采样数
  sample_n = NULL,               # 采样观测数（大数据）
  seed = NULL,
  ...
)
```

**返回值**: `ukb_shap` S3 对象

```r
list(
  shap_values = <matrix>,        # SHAP 值矩阵 (n x p)
  baseline = <numeric>,          # 基线值
  feature_names = c(...),        # 特征名
  feature_values = <data.frame>, # 原始特征值
  model = <ukb_ml>               # 原模型引用
)
```

#### 3.3.2 ukb_shap_summary() - SHAP 汇总图

```r
ukb_shap_summary(
  shap_object,                   # ukb_shap 对象
  max_features = 20,             # 最多显示特征数
  type = c("beeswarm", "bar"),   # 图类型
  color_palette = "viridis",
  ...
)
```

**输出**: ggplot2 图形

#### 3.3.3 ukb_shap_dependence() - SHAP 依赖图

```r
ukb_shap_dependence(
  shap_object,
  feature,                       # 主特征
  color_feature = NULL,          # 颜色特征（交互）
  alpha = 0.5,
  smooth = TRUE,
  ...
)
```

#### 3.3.4 ukb_shap_force() - 力图（单样本解释）

```r
ukb_shap_force(
  shap_object,
  row_id = 1,                    # 样本索引
  max_features = 10,
  ...
)
```

#### 3.3.5 ukb_shap_interaction() - 交互分析

```r
ukb_shap_interaction(
  object,                        # ukb_ml 对象
  data = NULL,
  features = NULL,               # 指定特征对
  ...
)
```

---

### 3.4 生存分析 ML API

#### 3.4.1 ukb_ml_survival() - 生存 ML 模型

```r
ukb_ml_survival(
  formula,                       # Surv(time, event) ~ x1 + x2 + ...
  data,
  model = c("rsf", "gbm_surv", "coxnet"),
  
  # RSF 参数
  ntree = 500,
  mtry = NULL,
  nodesize = 15,
  
  # GBM 参数
  n.trees = 500,
  interaction.depth = 3,
  shrinkage = 0.01,
  
  # CoxNet 参数
  alpha = 1,                     # 1=LASSO, 0=Ridge
  lambda = NULL,
  
  split_ratio = 0.8,
  seed = NULL,
  ...
)
```

**返回值**: `ukb_ml_surv` S3 对象

```r
list(
  model = <survival_model>,
  model_type = "rsf",
  formula = formula,
  train_data = <data>,
  test_data = <data>,
  c_index = 0.75,                # Harrell's C-index
  brier_score = <time_dependent>,
  call = <call>
)
```

#### 3.4.2 ukb_ml_survival_predict()

```r
ukb_ml_survival_predict(
  object,                        # ukb_ml_surv 对象
  newdata = NULL,
  times = c(1, 3, 5, 10),        # 预测时间点（年）
  type = c("risk", "survival", "chf")
)
```

#### 3.4.3 ukb_ml_survival_shap()

```r
ukb_ml_survival_shap(
  object,                        # ukb_ml_surv 对象
  data = NULL,
  time_point = 5,                # 特定时间点的 SHAP
  nsim = 50,
  ...
)
```

---

### 3.5 模型评估 API

#### 3.5.1 ukb_ml_metrics() - 评估指标

```r
ukb_ml_metrics(
  object,                        # ukb_ml 对象
  newdata = NULL,
  metrics = NULL                 # NULL = 自动选择
)
```

**分类任务返回**:
- AUC, Accuracy, Sensitivity, Specificity
- PPV, NPV, F1 Score
- Brier Score
- 95% CI (bootstrap)

**回归任务返回**:
- RMSE, MAE, R²
- 95% CI

#### 3.5.2 ukb_ml_roc() - ROC 曲线

```r
ukb_ml_roc(
  object,                        # ukb_ml 对象或列表
  newdata = NULL,
  plot = TRUE,
  ci = TRUE,
  ci_method = "delong"
)
```

#### 3.5.3 ukb_ml_calibration() - 校准曲线

```r
ukb_ml_calibration(
  object,
  newdata = NULL,
  n_bins = 10,
  method = c("loess", "isotonic"),
  plot = TRUE
)
```

#### 3.5.4 ukb_ml_confusion() - 混淆矩阵

```r
ukb_ml_confusion(
  object,
  newdata = NULL,
  threshold = 0.5,
  plot = TRUE
)
```

---

## 4. S3 方法

### 4.1 ukb_ml 类

```r
# 打印
print.ukb_ml(x, ...)

# 摘要
summary.ukb_ml(object, ...)

# 预测
predict.ukb_ml(object, newdata = NULL, type = "response", ...)

# 系数 (glmnet)
coef.ukb_ml(object, ...)

# 绘图
plot.ukb_ml(x, type = c("importance", "roc", "calibration"), ...)
```

### 4.2 ukb_shap 类

```r
print.ukb_shap(x, ...)
summary.ukb_shap(object, n = 10, ...)
plot.ukb_shap(x, type = c("summary", "bar"), ...)
```

### 4.3 ukb_ml_surv 类

```r
print.ukb_ml_surv(x, ...)
summary.ukb_ml_surv(object, ...)
predict.ukb_ml_surv(object, newdata, times, ...)
plot.ukb_ml_surv(x, type = c("importance", "survival"), ...)
```

---

## 5. 与现有模块集成

### 5.1 与倾向性评分集成

```r
# 使用 ML 估计倾向性评分
ps_ml <- ukb_ml_model(
  treatment ~ age + sex + bmi + smoking,
  data = ukb_data,
  model = "xgboost",
  task = "classification"
)

# 提取预测概率作为倾向性评分
ukb_data$ps_ml <- ukb_ml_predict(ps_ml, type = "prob")[, 2]

# 用于 IPTW
ps_result <- estimate_propensity(
  data = ukb_data,
  treatment = "treatment",
  covariates = c("age", "sex", "bmi", "smoking"),
  method = "custom",
  custom_ps = ukb_data$ps_ml
)
```

### 5.2 与多重插补集成

```r
# 对每个插补数据集训练模型
mi_models <- lapply(imputed_list, function(imp_data) {
  ukb_ml_model(
    outcome ~ .,
    data = imp_data,
    model = "rf"
  )
})

# 合并预测
mi_predictions <- pool_mi_predictions(mi_models, newdata = test_data)
```

### 5.3 与亚组分析集成

```r
# 在亚组内进行 ML 分析
subgroup_ml <- function(subgroup_data) {
  ukb_ml_model(
    outcome ~ predictors,
    data = subgroup_data,
    model = "rf"
  )
}

# 应用到所有亚组
results <- lapply(subgroups, subgroup_ml)
```

---

## 6. 大数据优化策略

### 6.1 采样策略

```r
# 训练时采样
ml_model <- ukb_ml_model(
  ...,
  sample_n = 50000,              # 从大数据中采样
  stratify = TRUE                # 分层采样保持类别比例
)

# SHAP 计算采样
shap_values <- ukb_shap(
  ml_model,
  sample_n = 1000                # SHAP 计算用 1000 样本
)
```

### 6.2 并行计算

```r
# 交叉验证并行
cv_result <- ukb_ml_cv(
  ...,
  parallel = TRUE,
  n_cores = 4
)

# 超参数调优并行
tune_result <- ukb_ml_tune(
  ...,
  parallel = TRUE
)
```

### 6.3 内存管理

```r
# 不保留训练数据
ml_model <- ukb_ml_model(
  ...,
  keep_data = FALSE
)

# 流式预测（分批）
predictions <- ukb_ml_predict_batch(
  model,
  newdata = large_data,
  batch_size = 10000
)
```

---

## 7. 可视化函数

### 7.1 plot_ml_importance() - 变量重要性

```r
plot_ml_importance(
  object,                        # ukb_ml 对象
  n_features = 20,
  type = c("bar", "dot"),
  color = "#3182BD"
)
```

### 7.2 plot_ml_roc() - ROC 曲线

```r
plot_ml_roc(
  ...,                           # 一个或多个 ukb_ml 对象
  models = list(),
  ci = TRUE,
  ci_alpha = 0.2
)
```

### 7.3 plot_ml_calibration() - 校准曲线

```r
plot_ml_calibration(
  object,
  n_bins = 10,
  smooth = TRUE,
  rug = TRUE
)
```

### 7.4 plot_ml_confusion() - 混淆矩阵热图

```r
plot_ml_confusion(
  object,
  threshold = 0.5,
  normalize = TRUE,
  colors = c("white", "#E34A33")
)
```

### 7.5 plot_ml_compare() - 模型比较

```r
plot_ml_compare(
  compare_object,                # ukb_ml_compare 结果
  metric = "auc",
  type = c("bar", "dot", "radar")
)
```

---

## 8. 执行计划

### Phase 1: 核心框架 (Week 1)

**文件**: `R/ml_model.R`

| 任务 | 函数 | 优先级 |
|------|------|--------|
| 统一模型接口 | `ukb_ml_model()` | P0 |
| 预测接口 | `ukb_ml_predict()` | P0 |
| S3 方法 | `print`, `summary`, `predict` | P0 |
| 变量重要性 | `ukb_ml_importance()` | P1 |

### Phase 2: 模型评估 (Week 1-2)

**文件**: `R/ml_evaluate.R`

| 任务 | 函数 | 优先级 |
|------|------|--------|
| 评估指标 | `ukb_ml_metrics()` | P0 |
| ROC 曲线 | `ukb_ml_roc()` | P0 |
| 校准曲线 | `ukb_ml_calibration()` | P1 |
| 混淆矩阵 | `ukb_ml_confusion()` | P1 |
| 模型比较 | `ukb_ml_compare()` | P1 |

### Phase 3: SHAP 解释 (Week 2)

**文件**: `R/ml_shap.R`

| 任务 | 函数 | 优先级 |
|------|------|--------|
| SHAP 计算 | `ukb_shap()` | P0 |
| SHAP 汇总图 | `ukb_shap_summary()` | P0 |
| SHAP 依赖图 | `ukb_shap_dependence()` | P1 |
| 力图 | `ukb_shap_force()` | P2 |

### Phase 4: 生存分析 ML (Week 3)

**文件**: `R/ml_survival.R`

| 任务 | 函数 | 优先级 |
|------|------|--------|
| 生存 ML | `ukb_ml_survival()` | P0 |
| 生存预测 | `ukb_ml_survival_predict()` | P0 |
| 生存 SHAP | `ukb_ml_survival_shap()` | P1 |

### Phase 5: 高级功能 (Week 3-4)

| 任务 | 函数 | 优先级 |
|------|------|--------|
| 交叉验证 | `ukb_ml_cv()` | P1 |
| 超参数调优 | `ukb_ml_tune()` | P2 |
| 可视化函数 | `plot_ml_*()` | P1 |

---

## 9. 测试计划

### 9.1 单元测试

```r
# tests/testthat/test-ml_model.R
test_that("ukb_ml_model works for classification", {
  # 测试各模型类型
})

test_that("ukb_ml_model works for regression", {
  # 测试回归任务
})

test_that("ukb_ml_predict returns correct format", {
  # 测试预测输出
})
```

### 9.2 集成测试

```r
# 测试完整工作流
test_that("complete ML workflow works", {
  # 训练 -> 预测 -> 评估 -> SHAP
})
```

---

## 10. 示例代码

### 10.1 基本使用

```r
library(UKBAnalytica)

# 准备数据
data <- ukb_data %>%
  select(cvd, age, sex, bmi, sbp, dbp, smoking, diabetes)

# 训练随机森林
ml_rf <- ukb_ml_model(
  cvd ~ .,
  data = data,
  model = "rf",
  task = "classification",
  split_ratio = 0.8,
  seed = 42
)

# 查看结果
print(ml_rf)
summary(ml_rf)

# 评估
metrics <- ukb_ml_metrics(ml_rf)
print(metrics)

# ROC 曲线
plot_ml_roc(ml_rf)
```

### 10.2 SHAP 分析

```r
# 计算 SHAP 值
shap <- ukb_shap(ml_rf, sample_n = 1000)

# SHAP 汇总图
ukb_shap_summary(shap)

# 单个变量的 SHAP 依赖
ukb_shap_dependence(shap, feature = "age", color_feature = "sex")
```

### 10.3 模型比较

```r
# 训练多个模型
ml_rf <- ukb_ml_model(cvd ~ ., data, model = "rf")
ml_xgb <- ukb_ml_model(cvd ~ ., data, model = "xgboost")
ml_lasso <- ukb_ml_model(cvd ~ ., data, model = "glmnet")

# 比较
comparison <- ukb_ml_compare(ml_rf, ml_xgb, ml_lasso)
print(comparison)
plot(comparison)
```

### 10.4 生存分析 ML

```r
# 随机生存森林
surv_rf <- ukb_ml_survival(
  Surv(time, event) ~ age + sex + bmi + smoking + diabetes,
  data = ukb_data,
  model = "rsf",
  ntree = 500
)

# 5年生存预测
pred <- ukb_ml_survival_predict(surv_rf, times = 5)

# 变量重要性
plot(surv_rf, type = "importance")
```

---

## 11. 注意事项

### 11.1 包检查兼容性

- 所有 ML 包放入 `Suggests`
- 使用时检查包是否安装
- 提供友好的安装提示

```r
.check_ml_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Package '%s' is required. Install with: install.packages('%s')",
      pkg, pkg
    ))
  }
}
```

### 11.2 大数据考虑

- 默认启用采样策略
- 提供内存使用估算
- 支持分批预测

### 11.3 可重复性

- 所有随机操作支持 `seed` 参数
- 保存完整的模型信息
- 支持模型导出/导入

---

## 12. 文件结构

```
R/
├── ml_model.R         # 核心模型函数
├── ml_evaluate.R      # 模型评估函数
├── ml_shap.R          # SHAP 解释函数
├── ml_survival.R      # 生存分析 ML
└── ml_utils.R         # 辅助函数

tests/testthat/
├── test-ml_model.R
├── test-ml_evaluate.R
├── test-ml_shap.R
└── test-ml_survival.R

man/
├── ukb_ml_model.Rd
├── ukb_shap.Rd
└── ...

vignettes/
└── machine-learning.Rmd
```

---

## 13. 版本规划

- **v0.6.0**: 核心 ML 模块 (Phase 1-2)
- **v0.7.0**: SHAP 解释 (Phase 3)
- **v0.8.0**: 生存分析 ML + 高级功能 (Phase 4-5)
