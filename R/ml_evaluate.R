#' Machine Learning Model Evaluation
#'
#' @description
#' Functions for evaluating and visualizing ML model performance.
#'
#' @name ml_evaluate
#' @keywords internal
NULL

# Metrics Function

#' Calculate Model Performance Metrics
#'
#' @description
#' Compute performance metrics for a trained ML model.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data for evaluation
#' @param metrics Specific metrics to compute (NULL for defaults)
#' @param ci Logical; compute confidence intervals (default FALSE)
#' @param ci_method Method for CI: "bootstrap" or "delong" (for AUC)
#' @param n_boot Number of bootstrap samples
#' @param verbose Print results
#' @param ... Additional arguments
#'
#' @return Named vector or list with metrics and optional CIs
#'
#' @export
ukb_ml_metrics <- function(object,
                           newdata = NULL,
                           metrics = NULL,
                           ci = FALSE,
                           ci_method = c("bootstrap", "delong"),
                           n_boot = 1000,
                           verbose = TRUE,
                           ...) {
  
  ci_method <- match.arg(ci_method)
  task <- object$task
  
  # Get predictions and true values
  if (is.null(newdata)) {
    y_true <- object$y_test
    pred <- ukb_ml_predict(object, type = "prob")
  } else {
    y_true <- newdata[[object$outcome]]
    if (task == "classification") y_true <- as.factor(y_true)
    pred <- ukb_ml_predict(object, newdata = newdata, type = "prob")
  }
  
  # Default metrics
  if (is.null(metrics)) {
    metrics <- if (task == "classification") {
      c("auc", "accuracy", "sensitivity", "specificity", "ppv", "npv", "f1", "brier")
    } else {
      c("rmse", "mae", "rsquared")
    }
  }
  
  # Calculate metrics
  result <- .calculate_all_metrics(y_true, pred, task, metrics)
  
  # Add confidence intervals if requested
  if (ci) {
    if (task == "classification" && ci_method == "delong" && "auc" %in% metrics) {
      .check_ml_package("pROC")
      prob <- if (is.matrix(pred)) pred[, 2] else pred
      roc_obj <- pROC::roc(y_true, prob, quiet = TRUE)
      auc_ci <- as.numeric(pROC::ci.auc(roc_obj))
      result <- list(
        metrics = result,
        ci = list(auc = c(lower = auc_ci[1], estimate = auc_ci[2], upper = auc_ci[3]))
      )
    } else {
      # Bootstrap CI
      boot_results <- .bootstrap_metrics(y_true, pred, task, metrics, n_boot)
      result <- list(
        metrics = result,
        ci = boot_results
      )
    }
  }
  
  if (verbose && !ci) {
    cat("\nModel Performance Metrics:\n")
    for (m in names(result)) {
      cat(sprintf("  %s: %.4f\n", m, result[m]))
    }
  }
  
  result
}

#' @keywords internal
.calculate_all_metrics <- function(y_true, y_pred, task, metrics) {
  result <- numeric()
  
  if (task == "classification") {
    .check_ml_package("pROC")
    
    # Get probability for positive class
    if (is.matrix(y_pred)) {
      prob <- y_pred[, 2]
    } else {
      prob <- y_pred
    }
    
    # Binary predictions
    pred_class <- ifelse(prob > 0.5, levels(y_true)[2], levels(y_true)[1])
    pred_class <- factor(pred_class, levels = levels(y_true))
    
    # Confusion matrix elements
    tp <- sum(pred_class == levels(y_true)[2] & y_true == levels(y_true)[2])
    tn <- sum(pred_class == levels(y_true)[1] & y_true == levels(y_true)[1])
    fp <- sum(pred_class == levels(y_true)[2] & y_true == levels(y_true)[1])
    fn <- sum(pred_class == levels(y_true)[1] & y_true == levels(y_true)[2])
    
    for (m in metrics) {
      val <- switch(m,
        auc = as.numeric(pROC::auc(pROC::roc(y_true, prob, quiet = TRUE))),
        accuracy = (tp + tn) / (tp + tn + fp + fn),
        sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA,
        specificity = if (tn + fp > 0) tn / (tn + fp) else NA,
        ppv = if (tp + fp > 0) tp / (tp + fp) else NA,
        npv = if (tn + fn > 0) tn / (tn + fn) else NA,
        f1 = if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA,
        brier = mean((prob - as.numeric(y_true == levels(y_true)[2]))^2),
        NA
      )
      result[m] <- val
    }
  } else {
    y_true <- as.numeric(y_true)
    y_pred <- as.numeric(y_pred)
    
    for (m in metrics) {
      val <- switch(m,
        rmse = sqrt(mean((y_true - y_pred)^2)),
        mae = mean(abs(y_true - y_pred)),
        rsquared = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2),
        mape = mean(abs((y_true - y_pred) / y_true)) * 100,
        NA
      )
      result[m] <- val
    }
  }
  
  result
}

#' @keywords internal
.bootstrap_metrics <- function(y_true, y_pred, task, metrics, n_boot) {
  n <- length(y_true)
  boot_mat <- matrix(NA, nrow = n_boot, ncol = length(metrics))
  colnames(boot_mat) <- metrics
  
  for (b in seq_len(n_boot)) {
    idx <- sample(n, replace = TRUE)
    y_b <- y_true[idx]
    pred_b <- if (is.matrix(y_pred)) y_pred[idx, ] else y_pred[idx]
    boot_mat[b, ] <- .calculate_all_metrics(y_b, pred_b, task, metrics)
  }
  
  # Calculate CIs
  ci_list <- list()
  for (m in metrics) {
    ci_list[[m]] <- c(
      lower = quantile(boot_mat[, m], 0.025, na.rm = TRUE),
      estimate = mean(boot_mat[, m], na.rm = TRUE),
      upper = quantile(boot_mat[, m], 0.975, na.rm = TRUE)
    )
  }
  
  ci_list
}

# ROC Curve

#' ROC Curve Analysis
#'
#' @description
#' Generate ROC curve and calculate AUC with optional confidence intervals.
#'
#' @param object A ukb_ml object or list of objects
#' @param newdata Optional new data
#' @param plot Whether to create ROC plot (default TRUE)
#' @param ci Compute confidence interval for AUC
#' @param ci_method Method: "delong" (default) or "bootstrap"
#' @param ... Additional arguments
#'
#' @return ukb_ml_roc object with ROC curve data
#'
#' @export
ukb_ml_roc <- function(object,
                       newdata = NULL,
                       plot = TRUE,
                       ci = TRUE,
                       ci_method = c("delong", "bootstrap"),
                       ...) {
  
  .check_ml_package("pROC")
  ci_method <- match.arg(ci_method)
  
  # Handle single model or list
  if (inherits(object, "ukb_ml")) {
    objects <- list(object)
    names(objects) <- .get_model_label(object$model_type)
  } else if (is.list(object)) {
    objects <- object
    if (is.null(names(objects))) {
      names(objects) <- sapply(objects, function(o) .get_model_label(o$model_type))
    }
  }
  
  roc_list <- list()
  auc_df <- data.frame(
    model = character(),
    auc = numeric(),
    lower = numeric(),
    upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(objects)) {
    obj <- objects[[name]]
    
    if (is.null(newdata)) {
      y_true <- obj$y_test
      pred <- ukb_ml_predict(obj, type = "prob")
    } else {
      y_true <- as.factor(newdata[[obj$outcome]])
      pred <- ukb_ml_predict(obj, newdata = newdata, type = "prob")
    }
    
    prob <- if (is.matrix(pred)) pred[, 2] else pred
    
    roc_obj <- pROC::roc(y_true, prob, quiet = TRUE)
    roc_list[[name]] <- roc_obj
    
    # AUC with CI
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    if (ci) {
      if (ci_method == "delong") {
        auc_ci <- as.numeric(pROC::ci.auc(roc_obj))
      } else {
        auc_ci <- as.numeric(pROC::ci.auc(roc_obj, method = "bootstrap"))
      }
      auc_df <- rbind(auc_df, data.frame(
        model = name,
        auc = auc_val,
        lower = auc_ci[1],
        upper = auc_ci[3],
        stringsAsFactors = FALSE
      ))
    } else {
      auc_df <- rbind(auc_df, data.frame(
        model = name,
        auc = auc_val,
        lower = NA,
        upper = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  result <- list(
    roc = roc_list,
    auc = auc_df
  )
  
  class(result) <- "ukb_ml_roc"
  
  if (plot) {
    result$plot <- plot_ml_roc(result)
  }
  
  result
}

#' @export
print.ukb_ml_roc <- function(x, ...) {
  cat("\nROC Analysis Results\n")
  cat("\n")
  print(x$auc, row.names = FALSE)
  invisible(x)
}

# Calibration Curve

#' Calibration Curve Analysis
#'
#' @description
#' Generate calibration curve to assess prediction reliability.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data
#' @param n_bins Number of bins for calibration (default 10)
#' @param method Smoothing method: "loess", "isotonic", or "none"
#' @param plot Whether to create calibration plot
#' @param ... Additional arguments
#'
#' @return ukb_ml_calibration object
#'
#' @export
ukb_ml_calibration <- function(object,
                               newdata = NULL,
                               n_bins = 10,
                               method = c("none", "loess", "isotonic"),
                               plot = TRUE,
                               ...) {
  
  method <- match.arg(method)
  
  if (object$task != "classification") {
    stop("Calibration curves are only for classification models")
  }
  
  # Get predictions
  if (is.null(newdata)) {
    y_true <- object$y_test
    pred <- ukb_ml_predict(object, type = "prob")
  } else {
    y_true <- as.factor(newdata[[object$outcome]])
    pred <- ukb_ml_predict(object, newdata = newdata, type = "prob")
  }
  
  prob <- if (is.matrix(pred)) pred[, 2] else pred
  y_binary <- as.numeric(y_true == levels(y_true)[2])
  
  # Create bins
  breaks <- seq(0, 1, length.out = n_bins + 1)
  bins <- cut(prob, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # Calculate observed vs predicted by bin
  cal_data <- data.frame(
    bin = seq_len(n_bins),
    predicted = numeric(n_bins),
    observed = numeric(n_bins),
    count = integer(n_bins),
    lower = numeric(n_bins),
    upper = numeric(n_bins)
  )
  
  for (i in seq_len(n_bins)) {
    idx <- which(bins == i)
    if (length(idx) > 0) {
      cal_data$predicted[i] <- mean(prob[idx])
      cal_data$observed[i] <- mean(y_binary[idx])
      cal_data$count[i] <- length(idx)
      
      # Wilson CI for proportion
      if (length(idx) >= 5) {
        prop <- mean(y_binary[idx])
        n <- length(idx)
        z <- 1.96
        denom <- 1 + z^2 / n
        center <- (prop + z^2 / (2 * n)) / denom
        margin <- z * sqrt((prop * (1 - prop) + z^2 / (4 * n)) / n) / denom
        cal_data$lower[i] <- max(0, center - margin)
        cal_data$upper[i] <- min(1, center + margin)
      } else {
        cal_data$lower[i] <- NA
        cal_data$upper[i] <- NA
      }
    } else {
      cal_data$predicted[i] <- (breaks[i] + breaks[i + 1]) / 2
      cal_data$observed[i] <- NA
      cal_data$count[i] <- 0
      cal_data$lower[i] <- NA
      cal_data$upper[i] <- NA
    }
  }
  
  # Remove empty bins
  cal_data <- cal_data[cal_data$count > 0, ]
  
  # Calculate calibration metrics
  # Brier score
  brier <- mean((prob - y_binary)^2)
  
  # Expected Calibration Error (ECE)
  ece <- sum(abs(cal_data$observed - cal_data$predicted) * cal_data$count) / sum(cal_data$count)
  
  result <- list(
    calibration = cal_data,
    brier_score = brier,
    ece = ece,
    method = method
  )
  
  class(result) <- "ukb_ml_calibration"
  
  if (plot) {
    result$plot <- plot_ml_calibration(result)
  }
  
  result
}

#' @export
print.ukb_ml_calibration <- function(x, ...) {
  cat("\nCalibration Analysis\n")
  cat("\n")
  cat(sprintf("Brier Score: %.4f\n", x$brier_score))
  cat(sprintf("Expected Calibration Error: %.4f\n", x$ece))
  cat("\nCalibration by bin:\n")
  print(x$calibration[, c("bin", "predicted", "observed", "count")], row.names = FALSE)
  invisible(x)
}

# Confusion Matrix

#' Confusion Matrix
#'
#' @description
#' Generate confusion matrix for classification model.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data
#' @param threshold Classification threshold (default 0.5)
#' @param plot Whether to create confusion matrix plot
#' @param ... Additional arguments
#'
#' @return ukb_ml_confusion object
#'
#' @export
ukb_ml_confusion <- function(object,
                             newdata = NULL,
                             threshold = 0.5,
                             plot = TRUE,
                             ...) {
  
  if (object$task != "classification") {
    stop("Confusion matrix is only for classification models")
  }
  
  # Get predictions
  if (is.null(newdata)) {
    y_true <- object$y_test
    pred <- ukb_ml_predict(object, type = "prob")
  } else {
    y_true <- as.factor(newdata[[object$outcome]])
    pred <- ukb_ml_predict(object, newdata = newdata, type = "prob")
  }
  
  prob <- if (is.matrix(pred)) pred[, 2] else pred
  pred_class <- ifelse(prob > threshold, levels(y_true)[2], levels(y_true)[1])
  pred_class <- factor(pred_class, levels = levels(y_true))
  
  # Create confusion matrix
  cm <- table(Predicted = pred_class, Actual = y_true)
  
  # Calculate metrics
  tp <- cm[2, 2]
  tn <- cm[1, 1]
  fp <- cm[2, 1]
  fn <- cm[1, 2]
  
  metrics <- c(
    accuracy = (tp + tn) / sum(cm),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv = tp / (tp + fp),
    npv = tn / (tn + fn),
    f1 = 2 * tp / (2 * tp + fp + fn)
  )
  
  result <- list(
    confusion_matrix = cm,
    metrics = metrics,
    threshold = threshold
  )
  
  class(result) <- "ukb_ml_confusion"
  
  if (plot) {
    result$plot <- plot_ml_confusion(result)
  }
  
  result
}

#' @export
print.ukb_ml_confusion <- function(x, ...) {
  cat("\nConfusion Matrix (threshold =", x$threshold, ")\n")
  cat("\n")
  print(x$confusion_matrix)
  cat("\nMetrics:\n")
  for (m in names(x$metrics)) {
    cat(sprintf("  %s: %.4f\n", m, x$metrics[m]))
  }
  invisible(x)
}
