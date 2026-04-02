#' Survival Machine Learning Module
#'
#' @description
#' Machine learning models for survival analysis including
#' random survival forests and gradient boosting survival models.
#'
#' @name ml_survival
#' @keywords internal
NULL

# Main Survival ML Function

#' Train Survival Machine Learning Model
#'
#' @description
#' Unified interface for training machine learning models for survival analysis.
#' Supports random survival forests (RSF), gradient boosting survival, and
#' regularized Cox models.
#'
#' @param formula Survival formula (e.g., Surv(time, event) ~ x1 + x2)
#' @param data Data frame
#' @param model Model type: "rsf" (random survival forest), "gbm_surv" (gradient boosting),
#'   "coxnet" (regularized Cox)
#' @param split_ratio Train/test split ratio (default 0.8)
#' @param seed Random seed
#' @param params List of model-specific parameters
#' @param verbose Print progress
#' @param ... Additional arguments
#'
#' @return A ukb_ml_surv object containing:
#' \itemize{
#'   \item model: Fitted survival model
#'   \item c_index: Harrell's C-index on test data
#'   \item train_data, test_data: Split datasets
#' }
#'
#' @examples
#' \dontrun{
#' # Random Survival Forest
#' surv_rf <- ukb_ml_survival(
#'   Surv(time, event) ~ age + sex + bmi + smoking,
#'   data = ukb_data,
#'   model = "rsf"
#' )
#'
#' # Prediction
#' pred <- ukb_ml_survival_predict(surv_rf, times = c(1, 3, 5))
#' }
#'
#' @export
ukb_ml_survival <- function(formula,
                            data,
                            model = c("rsf", "gbm_surv", "coxnet"),
                            split_ratio = 0.8,
                            seed = NULL,
                            params = list(),
                            verbose = TRUE,
                            ...) {
  
  model <- match.arg(model)
  
  .check_ml_package("survival")
  
  if (!is.null(seed)) set.seed(seed)
  
  # Parse formula
  terms_obj <- terms(formula, data = data)
  response_vars <- all.vars(formula[[2]])  # Surv(time, event)
  predictor_vars <- attr(terms_obj, "term.labels")
  
  # Handle "." notation
  if (length(predictor_vars) == 1 && predictor_vars == ".") {
    predictor_vars <- setdiff(names(data), response_vars)
  }
  
  # Get time and event variables
  time_var <- response_vars[1]
  event_var <- response_vars[2]
  
  # Remove incomplete cases
  all_vars <- c(time_var, event_var, predictor_vars)
  complete_idx <- complete.cases(data[, all_vars, drop = FALSE])
  clean_data <- data[complete_idx, , drop = FALSE]
  
  if (nrow(clean_data) < nrow(data)) {
    message(sprintf("Removed %d rows with missing values", nrow(data) - nrow(clean_data)))
  }
  
  # Split data
  n <- nrow(clean_data)
  train_idx <- sample(n, round(n * split_ratio))
  test_idx <- setdiff(seq_len(n), train_idx)
  
  train_data <- clean_data[train_idx, ]
  test_data <- clean_data[test_idx, ]
  
  if (verbose) {
    message(sprintf("Training %s survival model", model))
    message(sprintf("Train: %d, Test: %d observations", nrow(train_data), nrow(test_data)))
  }
  
  # Fit model
  fitted <- switch(model,
    rsf = .fit_rsf(formula, train_data, params, verbose, ...),
    gbm_surv = .fit_gbm_surv(formula, train_data, time_var, event_var, predictor_vars, params, verbose, ...),
    coxnet = .fit_coxnet(formula, train_data, time_var, event_var, predictor_vars, params, verbose, ...)
  )
  
  # Create result object
  result <- list(
    model = fitted,
    model_type = model,
    formula = formula,
    time_var = time_var,
    event_var = event_var,
    predictors = predictor_vars,
    train_data = train_data,
    test_data = test_data,
    train_idx = train_idx,
    test_idx = test_idx,
    seed = seed,
    call = match.call()
  )
  
  class(result) <- "ukb_ml_surv"
  
  # Calculate C-index
  result$c_index <- .calculate_c_index(result, test_data)
  
  if (verbose) {
    message(sprintf("Test C-index: %.3f", result$c_index))
  }
  
  result
}

# Model Fitting Functions

#' Fit Random Survival Forest
#' @keywords internal
.fit_rsf <- function(formula, data, params, verbose, ...) {
  .check_ml_package("randomForestSRC")
  
  default_params <- list(
    ntree = 500,
    mtry = NULL,
    nodesize = 15,
    importance = TRUE
  )
  
  model_params <- modifyList(default_params, params)
  
  do.call(randomForestSRC::rfsrc, c(
    list(formula = formula, data = data),
    model_params
  ))
}

#' Fit GBM Survival
#' @keywords internal
.fit_gbm_surv <- function(formula, data, time_var, event_var, predictor_vars, params, verbose, ...) {
  .check_ml_package("gbm")
  
  default_params <- list(
    distribution = "coxph",
    n.trees = 500,
    interaction.depth = 3,
    shrinkage = 0.01,
    n.minobsinnode = 10,
    cv.folds = 0
  )
  
  model_params <- modifyList(default_params, params)
  
  # GBM requires Surv object as response
  gbm_formula <- as.formula(paste0(
    "survival::Surv(", time_var, ", ", event_var, ") ~ ",
    paste(predictor_vars, collapse = " + ")
  ))
  
  do.call(gbm::gbm, c(
    list(formula = gbm_formula, data = data),
    model_params
  ))
}

#' Fit Cox with Elastic Net
#' @keywords internal
.fit_coxnet <- function(formula, data, time_var, event_var, predictor_vars, params, verbose, ...) {
  .check_ml_package("glmnet")
  
  default_params <- list(
    family = "cox",
    alpha = 1  # LASSO
  )
  
  model_params <- modifyList(default_params, params)
  
  # Prepare X matrix and Surv object
  X <- model.matrix(~ . - 1, data = data[, predictor_vars, drop = FALSE])
  y <- survival::Surv(data[[time_var]], data[[event_var]])
  
  do.call(glmnet::cv.glmnet, c(
    list(x = X, y = y),
    model_params
  ))
}

# Prediction Function

#' Predict from Survival ML Model
#'
#' @description
#' Generate predictions from a survival ML model.
#'
#' @param object A ukb_ml_surv object
#' @param newdata Optional new data
#' @param times Time points for survival prediction
#' @param type Prediction type: "risk", "survival", "chf" (cumulative hazard)
#' @param ... Additional arguments
#'
#' @return Matrix of predictions (observations x time points)
#'
#' @export
ukb_ml_survival_predict <- function(object,
                                    newdata = NULL,
                                    times = c(1, 3, 5, 10),
                                    type = c("survival", "risk", "chf"),
                                    ...) {
  
  type <- match.arg(type)
  
  if (is.null(newdata)) {
    newdata <- object$test_data
  }
  
  model <- object$model
  model_type <- object$model_type
  
  pred <- switch(model_type,
    rsf = .predict_rsf(model, newdata, times, type),
    gbm_surv = .predict_gbm_surv(model, newdata, times, type, object),
    coxnet = .predict_coxnet(model, newdata, times, type, object)
  )
  
  pred
}

#' @keywords internal
.predict_rsf <- function(model, newdata, times, type) {
  pred <- predict(model, newdata = newdata)
  
  # Get survival probabilities at specified times
  time_idx <- sapply(times, function(t) {
    which.min(abs(pred$time.interest - t))
  })
  
  if (type == "survival") {
    pred$survival[, time_idx, drop = FALSE]
  } else if (type == "chf") {
    pred$chf[, time_idx, drop = FALSE]
  } else {
    # Risk scores (higher = worse prognosis)
    pred$predicted
  }
}

#' @keywords internal
.predict_gbm_surv <- function(model, newdata, times, type, object) {
  .check_ml_package("gbm")
  
  # GBM returns relative risk
  risk <- predict(model, newdata = newdata, n.trees = model$n.trees, type = "link")
  
  if (type == "risk") {
    return(risk)
  }
  
  # For survival probabilities, need baseline hazard estimation
  # Using Cox model baseline
  .check_ml_package("survival")
  
  train_data <- object$train_data
  train_risk <- predict(model, newdata = train_data, n.trees = model$n.trees, type = "link")
  
  # Fit baseline hazard
  surv_obj <- survival::Surv(train_data[[object$time_var]], train_data[[object$event_var]])
  basehaz_fit <- survival::basehaz(survival::coxph(surv_obj ~ offset(train_risk), data = train_data))
  
  # Interpolate baseline hazard at requested times
  H0 <- approx(basehaz_fit$time, basehaz_fit$hazard, times, rule = 2)$y
  
  # Calculate survival or CHF
  surv_mat <- matrix(NA, nrow = length(risk), ncol = length(times))
  
  for (i in seq_along(times)) {
    if (type == "survival") {
      surv_mat[, i] <- exp(-H0[i] * exp(risk))
    } else {
      surv_mat[, i] <- H0[i] * exp(risk)
    }
  }
  
  colnames(surv_mat) <- paste0("t", times)
  surv_mat
}

#' @keywords internal
.predict_coxnet <- function(model, newdata, times, type, object) {
  .check_ml_package("glmnet")
  .check_ml_package("survival")
  
  X <- model.matrix(~ . - 1, data = newdata[, object$predictors, drop = FALSE])
  risk <- as.numeric(predict(model, newx = X, s = "lambda.min", type = "link"))
  
  if (type == "risk") {
    return(risk)
  }
  
  # Estimate baseline hazard from training data
  train_X <- model.matrix(~ . - 1, data = object$train_data[, object$predictors, drop = FALSE])
  train_risk <- as.numeric(predict(model, newx = train_X, s = "lambda.min", type = "link"))
  
  surv_obj <- survival::Surv(object$train_data[[object$time_var]], object$train_data[[object$event_var]])
  basehaz_fit <- survival::basehaz(survival::coxph(surv_obj ~ offset(train_risk), data = object$train_data))
  
  H0 <- approx(basehaz_fit$time, basehaz_fit$hazard, times, rule = 2)$y
  
  surv_mat <- matrix(NA, nrow = length(risk), ncol = length(times))
  
  for (i in seq_along(times)) {
    if (type == "survival") {
      surv_mat[, i] <- exp(-H0[i] * exp(risk))
    } else {
      surv_mat[, i] <- H0[i] * exp(risk)
    }
  }
  
  colnames(surv_mat) <- paste0("t", times)
  surv_mat
}

# C-index Calculation

#' @keywords internal
.calculate_c_index <- function(object, test_data) {
  .check_ml_package("survival")
  
  # Get risk scores
  risk <- switch(object$model_type,
    rsf = predict(object$model, newdata = test_data)$predicted,
    gbm_surv = predict(object$model, newdata = test_data, n.trees = object$model$n.trees, type = "link"),
    coxnet = {
      X <- model.matrix(~ . - 1, data = test_data[, object$predictors, drop = FALSE])
      as.numeric(predict(object$model, newx = X, s = "lambda.min", type = "link"))
    }
  )
  
  surv_obj <- survival::Surv(test_data[[object$time_var]], test_data[[object$event_var]])
  
  # Calculate concordance
  conc <- survival::concordance(surv_obj ~ risk)
  conc$concordance
}

# Variable Importance for Survival

#' Get Variable Importance for Survival Model
#'
#' @param object A ukb_ml_surv object
#' @param ... Additional arguments
#'
#' @return Data frame with variable importance
#'
#' @export
ukb_ml_survival_importance <- function(object, ...) {
  
  model <- object$model
  model_type <- object$model_type
  
  importance <- switch(model_type,
    rsf = {
      imp <- model$importance
      data.frame(
        variable = names(imp),
        importance = as.numeric(imp),
        stringsAsFactors = FALSE
      )
    },
    gbm_surv = {
      .check_ml_package("gbm")
      imp <- summary(model, plotit = FALSE)
      data.frame(
        variable = imp$var,
        importance = imp$rel.inf,
        stringsAsFactors = FALSE
      )
    },
    coxnet = {
      coefs <- as.matrix(coef(model, s = "lambda.min"))
      data.frame(
        variable = rownames(coefs),
        importance = abs(coefs[, 1]),
        stringsAsFactors = FALSE
      )
    }
  )
  
  importance <- importance[order(importance$importance, decreasing = TRUE), ]
  rownames(importance) <- NULL
  
  importance
}

# SHAP for Survival Models

#' SHAP Values for Survival Models
#'
#' @description
#' Compute SHAP values for survival ML models at a specific time point.
#'
#' @param object A ukb_ml_surv object
#' @param data Data for SHAP computation
#' @param time_point Time point for SHAP calculation
#' @param nsim Number of Monte Carlo samples
#' @param sample_n Subsample size
#' @param seed Random seed
#' @param verbose Print progress
#' @param ... Additional arguments
#'
#' @return A ukb_shap object
#'
#' @export
ukb_ml_survival_shap <- function(object,
                                 data = NULL,
                                 time_point = 5,
                                 nsim = 50,
                                 sample_n = NULL,
                                 seed = NULL,
                                 verbose = TRUE,
                                 ...) {
  
  .check_ml_package("fastshap")
  
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(data)) {
    data <- object$test_data
  }
  
  X <- data[, object$predictors, drop = FALSE]
  
  if (!is.null(sample_n) && sample_n < nrow(X)) {
    idx <- sample(nrow(X), sample_n)
    X <- X[idx, , drop = FALSE]
    data <- data[idx, , drop = FALSE]
    if (verbose) message(sprintf("Using %d sampled observations for SHAP", sample_n))
  }
  
  if (verbose) message(sprintf("Computing SHAP values at time = %g...", time_point))
  
  # Create prediction wrapper for survival probability at time_point
  pred_wrapper <- function(model, newdata) {
    # This is a simplified approach - for RSF
    if (object$model_type == "rsf") {
      pred <- predict(model, newdata = newdata)
      time_idx <- which.min(abs(pred$time.interest - time_point))
      pred$survival[, time_idx]
    } else {
      # For other models, compute survival at time_point
      full_data <- cbind(newdata, data[, c(object$time_var, object$event_var)])
      ukb_ml_survival_predict(object, newdata = full_data, times = time_point, type = "survival")[, 1]
    }
  }
  
  # Compute SHAP
  shap_values <- fastshap::explain(
    object = object$model,
    feature_names = object$predictors,
    X = X,
    pred_wrapper = pred_wrapper,
    nsim = nsim,
    ...
  )
  
  result <- list(
    shap_values = as.matrix(shap_values),
    baseline = NA,  # Baseline is complex for survival
    feature_names = object$predictors,
    feature_values = X,
    model_type = object$model_type,
    task = "survival",
    time_point = time_point
  )
  
  class(result) <- "ukb_shap"
  
  if (verbose) message("SHAP computation complete")
  
  result
}

# S3 Methods

#' @export
print.ukb_ml_surv <- function(x, ...) {
  cat("\n")
  cat("UKB Survival Machine Learning Model\n")
  cat("\n")
  cat(sprintf("Model: %s\n", x$model_type))
  cat(sprintf("Time variable: %s\n", x$time_var))
  cat(sprintf("Event variable: %s\n", x$event_var))
  cat(sprintf("Predictors: %d variables\n", length(x$predictors)))
  cat(sprintf("Train size: %d\n", nrow(x$train_data)))
  cat(sprintf("Test size: %d\n", nrow(x$test_data)))
  cat(sprintf("\nTest C-index: %.3f\n", x$c_index))
  invisible(x)
}

#' @export
summary.ukb_ml_surv <- function(object, ...) {
  print(object)
  
  cat("\nVariable Importance (Top 10):\n")
  imp <- ukb_ml_survival_importance(object)
  print(head(imp, 10))
  
  invisible(object)
}

#' @export
predict.ukb_ml_surv <- function(object, newdata = NULL, times = c(1, 3, 5), type = "survival", ...) {
  ukb_ml_survival_predict(object, newdata = newdata, times = times, type = type, ...)
}
