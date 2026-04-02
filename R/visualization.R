#' @importFrom stats relevel
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh geom_vline labs
#'   theme_minimal theme element_text element_blank scale_x_log10
#'   geom_step geom_ribbon geom_hline annotate scale_color_manual
#'   scale_fill_manual coord_cartesian scale_x_continuous scale_y_continuous
#'   geom_histogram geom_density geom_segment geom_errorbar geom_smooth
#'   geom_abline coord_equal scale_size_continuous
NULL

#' Plot Forest Plot for Subgroup Analysis
#'
#' @description
#' Create a forest plot to visualize subgroup analysis results with effect estimates
#' and confidence intervals.
#'
#' @param results A data.frame from \code{run_subgroup_analysis()} or \code{run_multi_subgroup()}.
#' @param estimate_col Character string specifying the column name for effect estimates.
#'   Default "estimate".
#' @param lower_col Character string specifying the column for lower CI. Default "lower95".
#' @param upper_col Character string specifying the column for upper CI. Default "upper95".
#' @param label_col Character string specifying the column for subgroup labels. Default "subgroup".
#' @param pvalue_col Character string specifying the column for p-values. Default "pvalue".
#' @param p_interaction_col Character string for interaction p-value column. Default "p_interaction".
#' @param null_value Numeric value for the null effect line. Default 1 (for HR/OR).
#' @param log_scale Logical; whether to use log scale for x-axis. Default TRUE.
#' @param colors Character vector of colors. Default NULL uses ggplot2 defaults.
#' @param title Character string for plot title. Default "Subgroup Analysis".
#' @param xlab Character string for x-axis label. Default "Hazard Ratio (95\% CI)".
#' @param show_n Logical; whether to show sample size. Default TRUE.
#' @param show_events Logical; whether to show event count. Default TRUE.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' results <- run_multi_subgroup(data, exposure = "smoking",
#'                               subgroup_vars = c("sex", "age_group"),
#'                               model_type = "cox", endpoint = c("time", "status"))
#' plot_forest(results)
#' }
#'
#' @import ggplot2
#' @export
plot_forest <- function(results,
                         estimate_col = "estimate",
                         lower_col = "lower95",
                         upper_col = "upper95",
                         label_col = "subgroup",
                         pvalue_col = "pvalue",
                         p_interaction_col = "p_interaction",
                         null_value = 1,
                         log_scale = TRUE,
                         colors = NULL,
                         title = "Subgroup Analysis",
                         xlab = "Hazard Ratio (95% CI)",
                         show_n = TRUE,
                         show_events = TRUE) {

  # Validate inputs
  if (!is.data.frame(results)) {
    stop("'results' must be a data.frame.")
  }

  required_cols <- c(estimate_col, lower_col, upper_col, label_col)
  missing_cols <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0) {
    stop(sprintf("Required columns not found: %s", paste(missing_cols, collapse = ", ")))
  }

  # Prepare data for plotting
  plot_data <- data.frame(
    subgroup = results[[label_col]],
    estimate = results[[estimate_col]],
    lower95 = results[[lower_col]],
    upper95 = results[[upper_col]],
    stringsAsFactors = FALSE
  )

  # Add optional columns
  if ("n" %in% names(results)) {
    plot_data$n_total <- results$n
  } else {
    plot_data$n_total <- NA
  }

  if ("n_event" %in% names(results)) {
    plot_data$n_event <- results$n_event
  } else {
    plot_data$n_event <- NA
  }

  if (pvalue_col %in% names(results)) {
    plot_data$pvalue <- results[[pvalue_col]]
  }

  if (p_interaction_col %in% names(results)) {
    plot_data$p_interaction <- results[[p_interaction_col]]
  }

  # Create CI label
  plot_data$ci_label <- sprintf("%.2f (%.2f-%.2f)",
                                 plot_data$estimate,
                                 plot_data$lower95,
                                 plot_data$upper95)

  # Reverse order for plotting (top to bottom)
  plot_data$subgroup <- factor(plot_data$subgroup,
                                levels = rev(unique(plot_data$subgroup)))

  # Remove rows with NA estimates
  plot_data <- plot_data[!is.na(plot_data$estimate), ]

  # Create forest plot
  p <- ggplot(plot_data, aes(x = estimate, y = subgroup)) +
    # Point estimates
    geom_point(size = 3, shape = 18) +
    # Confidence intervals
    geom_errorbarh(aes(xmin = lower95, xmax = upper95),
                   height = 0.2, linewidth = 0.5) +
    # Null effect line
    geom_vline(xintercept = null_value, linetype = "dashed", color = "gray50") +
    # Labels
    labs(title = title, x = xlab, y = NULL) +
    # Theme
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(hjust = 0)
    )

  # Log scale
  if (log_scale && null_value > 0) {
    # Calculate reasonable breaks
    x_range <- range(c(plot_data$lower95, plot_data$upper95), na.rm = TRUE)
    breaks <- c(0.25, 0.5, 1, 2, 4)
    breaks <- breaks[breaks >= x_range[1] * 0.8 & breaks <= x_range[2] * 1.2]
    if (length(breaks) < 3) breaks <- c(0.5, 1, 2)

    p <- p + scale_x_log10(breaks = breaks)
  }

  # Add custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

  return(p)
}


#' Plot Kaplan-Meier Survival Curve
#'
#' @description
#' Create a Kaplan-Meier survival curve with optional risk table and log-rank p-value.
#'
#' @param data A data.frame or data.table containing survival data.
#' @param time_col Character string specifying the time column name.
#' @param status_col Character string specifying the event status column name.
#' @param group_col Character string specifying the grouping variable. Default NULL for overall curve.
#' @param conf_int Logical; whether to show confidence intervals. Default TRUE.
#' @param risk_table Logical; whether to show number at risk table. Default TRUE.
#' @param censor_marks Logical; whether to show censoring marks. Default TRUE.
#' @param palette Character string specifying color palette. Default "jco".
#'   Options: "jco", "nejm", "lancet", "npg", or custom color vector.
#' @param title Character string for plot title. Default NULL.
#' @param xlab Character string for x-axis label. Default "Time (years)".
#' @param ylab Character string for y-axis label. Default "Survival Probability".
#' @param legend_title Character string for legend title. Default "Group".
#' @param median_line Logical; whether to show median survival line. Default TRUE.
#' @param pvalue Logical; whether to show log-rank p-value. Default TRUE.
#' @param xlim Numeric vector of length 2 for x-axis limits. Default NULL.
#' @param break_time Numeric value for x-axis tick interval. Default NULL.
#'
#' @return A ggplot2 object (or a list with plot and risk table if risk_table = TRUE).
#'
#' @examples
#' \dontrun{
#' library(survival)
#' plot_km_curve(lung, time_col = "time", status_col = "status", group_col = "sex")
#' }
#'
#' @import ggplot2
#' @importFrom stats as.formula
#' @export
plot_km_curve <- function(data,
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
                           break_time = NULL) {

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for Kaplan-Meier curves.")
  }

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }

  if (!time_col %in% names(data) || !status_col %in% names(data)) {
    stop("time_col and status_col must exist in data.")
  }

  if (!is.null(group_col) && !group_col %in% names(data)) {
    stop(sprintf("group_col '%s' not found in data.", group_col))
  }

  # Build survival formula
  if (is.null(group_col)) {
    formula_obj <- stats::as.formula(
      paste0("survival::Surv(", time_col, ", ", status_col, ") ~ 1")
    )
  } else {
    formula_obj <- stats::as.formula(
      paste0("survival::Surv(", time_col, ", ", status_col, ") ~ ", group_col)
    )
  }

  # Fit survival model
  fit <- survival::survfit(formula_obj, data = data)

  # Extract survival data for plotting
  if (is.null(group_col)) {
    surv_data <- data.frame(
      time = fit$time,
      surv = fit$surv,
      upper = fit$upper,
      lower = fit$lower,
      n.risk = fit$n.risk,
      n.event = fit$n.event,
      strata = "Overall"
    )
  } else {
    # Multiple strata
    strata_names <- names(fit$strata)
    surv_list <- lapply(seq_along(fit$strata), function(i) {
      idx_start <- if (i == 1) 1 else sum(fit$strata[1:(i-1)]) + 1
      idx_end <- sum(fit$strata[1:i])
      data.frame(
        time = fit$time[idx_start:idx_end],
        surv = fit$surv[idx_start:idx_end],
        upper = fit$upper[idx_start:idx_end],
        lower = fit$lower[idx_start:idx_end],
        n.risk = fit$n.risk[idx_start:idx_end],
        n.event = fit$n.event[idx_start:idx_end],
        strata = strata_names[i]
      )
    })
    surv_data <- do.call(rbind, surv_list)

    # Clean strata names
    surv_data$strata <- gsub(paste0("^", group_col, "="), "", surv_data$strata)
  }

  # Get colors
  n_groups <- length(unique(surv_data$strata))
  if (is.character(palette) && length(palette) == 1) {
    if (palette == "jco") {
      colors <- c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC", "#003C67", "#8F7700", "#3B3B3B")
    } else if (palette == "nejm") {
      colors <- c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97")
    } else if (palette == "lancet") {
      colors <- c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6")
    } else {
      colors <- scales::hue_pal()(n_groups)
    }
  } else {
    colors <- palette
  }
  colors <- colors[1:n_groups]

  # Calculate log-rank p-value
  p_val <- NA
  if (pvalue && !is.null(group_col)) {
    logrank <- survival::survdiff(formula_obj, data = data)
    p_val <- 1 - stats::pchisq(logrank$chisq, df = length(logrank$n) - 1)
  }

  # Create main plot
  p <- ggplot(surv_data, aes(x = time, y = surv, color = strata)) +
    geom_step(linewidth = 1) +
    scale_color_manual(values = colors, name = legend_title) +
    labs(title = title, x = xlab, y = ylab) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  # Add confidence interval
  if (conf_int) {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, fill = strata),
                         alpha = 0.2, color = NA) +
      scale_fill_manual(values = colors, guide = "none")
  }

  # Add censoring marks
  if (censor_marks) {
    censor_data <- surv_data[surv_data$n.event == 0, ]
    if (nrow(censor_data) > 0) {
      p <- p + geom_point(data = censor_data, aes(x = time, y = surv),
                          shape = 3, size = 2)
    }
  }

  # Add median survival line
  if (median_line) {
    p <- p + geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50")
  }

  # Set x-axis limits and breaks
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  if (!is.null(break_time)) {
    max_time <- max(surv_data$time, na.rm = TRUE)
    p <- p + scale_x_continuous(breaks = seq(0, max_time, by = break_time))
  }

  # Add p-value annotation
  if (pvalue && !is.na(p_val)) {
    p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
    p <- p + annotate("text", x = max(surv_data$time) * 0.8, y = 0.9,
                      label = paste("Log-rank", p_label),
                      hjust = 1, size = 4)
  }

  # Create risk table if requested
  if (risk_table) {
    # Calculate number at risk at regular intervals
    max_time <- max(surv_data$time, na.rm = TRUE)
    time_points <- if (!is.null(break_time)) {
      seq(0, max_time, by = break_time)
    } else {
      seq(0, max_time, length.out = 6)
    }

    risk_data <- expand.grid(
      time = time_points,
      strata = unique(surv_data$strata),
      stringsAsFactors = FALSE
    )

    risk_data$n.risk <- sapply(1:nrow(risk_data), function(i) {
      t <- risk_data$time[i]
      s <- risk_data$strata[i]
      subset_data <- surv_data[surv_data$strata == s & surv_data$time >= t, ]
      if (nrow(subset_data) > 0) {
        max(subset_data$n.risk)
      } else {
        0
      }
    })

    risk_table_plot <- ggplot(risk_data, aes(x = time, y = strata, label = n.risk)) +
      geom_text(size = 3) +
      scale_x_continuous(limits = c(0, max_time)) +
      labs(x = NULL, y = "Number at Risk") +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()
      )

    # Return both plots (user can combine with cowplot/patchwork)
    return(list(km_plot = p, risk_table = risk_table_plot))
  }

  return(p)
}


#' Plot Propensity Score Distribution
#'
#' @description
#' Visualize the distribution of propensity scores by treatment group.
#'
#' @param data A data.frame or data.table containing propensity scores.
#' @param ps_col Character string specifying the PS column name. Default "ps".
#' @param treatment Character string specifying the treatment variable name.
#' @param type Character string specifying plot type: "histogram", "density", or "mirror".
#' @param matched Logical; whether to show matched vs unmatched. Default FALSE.
#' @param match_col Character string for the matching indicator column. Default NULL.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' ps_data <- estimate_propensity_score(mydata, "treated", c("age", "sex"))
#' plot_ps_distribution(ps_data, treatment = "treated", type = "mirror")
#' }
#'
#' @import ggplot2
#' @export
plot_ps_distribution <- function(data,
                                   ps_col = "ps",
                                   treatment,
                                   type = c("histogram", "density", "mirror"),
                                   matched = FALSE,
                                   match_col = NULL) {

  type <- match.arg(type)

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }

  if (!ps_col %in% names(data)) {
    stop(sprintf("PS column '%s' not found in data.", ps_col))
  }

  if (!treatment %in% names(data)) {
    stop(sprintf("Treatment variable '%s' not found in data.", treatment))
  }

  # Prepare data
  plot_data <- data.frame(
    ps = data[[ps_col]],
    group = factor(data[[treatment]], levels = c(0, 1),
                   labels = c("Control", "Treated"))
  )

  if (matched && !is.null(match_col)) {
    plot_data$method <- ifelse(!is.na(data[[match_col]]), "Matched", "Unmatched")
  }

  # Colors
  colors <- c("Control" = "#0073C2", "Treated" = "#CD534C")

  if (type == "histogram") {
    p <- ggplot(plot_data, aes(x = ps, fill = group)) +
      geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
      scale_fill_manual(values = colors) +
      labs(title = "Propensity Score Distribution",
           x = "Propensity Score", y = "Count", fill = "Group") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )

  } else if (type == "density") {
    p <- ggplot(plot_data, aes(x = ps, fill = group, color = group)) +
      geom_density(alpha = 0.4) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      labs(title = "Propensity Score Distribution",
           x = "Propensity Score", y = "Density", fill = "Group", color = "Group") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )

  } else {
    # Mirror histogram
    plot_data$count <- 1
    plot_data$count[plot_data$group == "Control"] <- -1

    p <- ggplot(plot_data, aes(x = ps, fill = group)) +
      geom_histogram(data = plot_data[plot_data$group == "Treated", ],
                     aes(y = after_stat(count)), bins = 30, alpha = 0.8) +
      geom_histogram(data = plot_data[plot_data$group == "Control", ],
                     aes(y = -after_stat(count)), bins = 30, alpha = 0.8) +
      scale_fill_manual(values = colors) +
      scale_y_continuous(labels = abs) +
      labs(title = "Propensity Score Distribution (Mirror Plot)",
           x = "Propensity Score", y = "Count", fill = "Group") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  }

  return(p)
}


#' Plot Covariate Balance (Love Plot)
#'
#' @description
#' Create a Love plot comparing standardized mean differences before and after matching/weighting.
#'
#' @param balance_before A data.frame from \code{assess_balance()} for unmatched data.
#' @param balance_after A data.frame from \code{assess_balance()} for matched/weighted data.
#' @param threshold Numeric threshold for balance (vertical lines). Default 0.1.
#' @param title Character string for plot title. Default "Covariate Balance".
#' @param xlab Character string for x-axis label. Default "Standardized Mean Difference".
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' balance_before <- assess_balance(data, "treated", c("age", "sex"), method = "unmatched")
#' balance_after <- assess_balance(matched_data, "treated", c("age", "sex"), method = "matched")
#' plot_balance(balance_before, balance_after)
#' }
#'
#' @import ggplot2
#' @export
plot_balance <- function(balance_before,
                          balance_after,
                          threshold = 0.1,
                          title = "Covariate Balance",
                          xlab = "Standardized Mean Difference") {

  # Validate inputs
  if (!is.data.frame(balance_before) || !is.data.frame(balance_after)) {
    stop("Both balance_before and balance_after must be data.frames.")
  }

  if (!"variable" %in% names(balance_before) || !"smd" %in% names(balance_before)) {
    stop("balance_before must contain 'variable' and 'smd' columns.")
  }

  if (!"variable" %in% names(balance_after) || !"smd" %in% names(balance_after)) {
    stop("balance_after must contain 'variable' and 'smd' columns.")
  }

  # Merge data
  plot_data <- merge(
    balance_before[, c("variable", "smd")],
    balance_after[, c("variable", "smd")],
    by = "variable",
    suffixes = c("_before", "_after")
  )

  names(plot_data) <- c("variable", "smd_before", "smd_after")

  # Create plot
  p <- ggplot(plot_data) +
    # Before matching (hollow circles)
    geom_point(aes(x = smd_before, y = variable), shape = 1, size = 3, color = "#CD534C") +
    # After matching (filled circles)
    geom_point(aes(x = smd_after, y = variable), shape = 16, size = 3, color = "#0073C2") +
    # Connect with lines
    geom_segment(aes(x = smd_before, xend = smd_after, y = variable, yend = variable),
                 color = "gray60", linewidth = 0.5) +
    # Threshold lines
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, color = "black") +
    # Labels
    labs(title = title, x = xlab, y = NULL) +
    # Theme
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    # Legend (manual)
    annotate("point", x = Inf, y = Inf, shape = 1, size = 3, color = "#CD534C") +
    annotate("point", x = Inf, y = Inf, shape = 16, size = 3, color = "#0073C2")

  return(p)
}


#' Plot Calibration Curve
#'
#' @description
#' Create a calibration plot comparing predicted probabilities to observed outcomes.
#'
#' @param data A data.frame or data.table.
#' @param predicted Character string specifying the column with predicted probabilities.
#' @param observed Character string specifying the column with observed binary outcomes.
#' @param n_bins Integer number of bins for calibration. Default 10.
#' @param smooth Logical; whether to add a smooth calibration line. Default TRUE.
#' @param conf_int Logical; whether to show confidence intervals. Default TRUE.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' plot_calibration(data, predicted = "pred_prob", observed = "outcome")
#' }
#'
#' @import ggplot2
#' @importFrom stats quantile binom.test
#' @export
plot_calibration <- function(data,
                              predicted,
                              observed,
                              n_bins = 10,
                              smooth = TRUE,
                              conf_int = TRUE) {

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }

  if (!predicted %in% names(data)) {
    stop(sprintf("predicted column '%s' not found in data.", predicted))
  }

  if (!observed %in% names(data)) {
    stop(sprintf("observed column '%s' not found in data.", observed))
  }

  # Create bins
  pred_vals <- data[[predicted]]
  obs_vals <- data[[observed]]

  breaks <- stats::quantile(pred_vals, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  breaks <- unique(breaks)
  bins <- cut(pred_vals, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  # Calculate observed proportion and mean predicted in each bin
  calib_data <- data.frame(
    bin = 1:length(unique(bins[!is.na(bins)])),
    predicted = tapply(pred_vals, bins, mean, na.rm = TRUE),
    observed = tapply(obs_vals, bins, mean, na.rm = TRUE),
    count = tapply(obs_vals, bins, length)
  )

  # Calculate confidence intervals
  if (conf_int) {
    n_events <- tapply(obs_vals, bins, sum, na.rm = TRUE)
    n_total <- calib_data$count

    calib_data$lower <- mapply(function(x, n) {
      if (n > 0) stats::binom.test(x, n)$conf.int[1] else NA
    }, n_events, n_total)

    calib_data$upper <- mapply(function(x, n) {
      if (n > 0) stats::binom.test(x, n)$conf.int[2] else NA
    }, n_events, n_total)
  }

  # Create plot
  p <- ggplot(calib_data, aes(x = predicted, y = observed)) +
    # Perfect calibration line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    # Points
    geom_point(aes(size = count), color = "#0073C2") +
    scale_size_continuous(name = "N", range = c(2, 6)) +
    # Labels
    labs(
      title = "Calibration Plot",
      x = "Mean Predicted Probability",
      y = "Observed Proportion"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  # Add confidence intervals
  if (conf_int) {
    p <- p + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, color = "#0073C2")
  }

  # Add smooth line
  if (smooth && nrow(calib_data) >= 3) {
    p <- p + geom_smooth(method = "loess", se = FALSE, color = "#CD534C", linewidth = 1)
  }

  return(p)
}


#' Plot Mediation Analysis Results
#'
#' @description
#' Create visualizations for mediation analysis results, including path diagrams,
#' effect bar charts, and decomposition plots.
#'
#' @param mediation_result An object of class "mediation_result" from \code{run_mediation()}.
#' @param type Character string specifying plot type:
#'   \describe{
#'     \item{"path"}{Path diagram showing exposure -> mediator -> outcome}
#'     \item{"effects"}{Bar chart of effect estimates with confidence intervals}
#'     \item{"decomposition"}{Pie/bar chart showing effect decomposition (NDE vs NIE)}
#'   }
#' @param show_ci Logical; whether to show confidence intervals. Default TRUE.
#' @param show_pvalue Logical; whether to show p-values. Default TRUE.
#' @param exponentiate Logical; whether to exponentiate estimates (for HR/OR). Default FALSE.
#' @param title Character string for plot title. Default NULL (auto-generated).
#' @param colors Character vector of colors. Default NULL uses package defaults.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' med_result <- run_mediation(data, "exposure", "mediator", "outcome")
#' plot_mediation(med_result, type = "effects")
#' plot_mediation(med_result, type = "decomposition")
#' }
#'
#' @import ggplot2
#' @export
plot_mediation <- function(mediation_result,
                            type = c("effects", "path", "decomposition"),
                            show_ci = TRUE,
                            show_pvalue = TRUE,
                            exponentiate = FALSE,
                            title = NULL,
                            colors = NULL) {

  type <- match.arg(type)

  if (!inherits(mediation_result, "mediation_result")) {
    stop("'mediation_result' must be an object from run_mediation().")
  }

  effects <- mediation_result$effects
  params <- mediation_result$params

  # Default colors
  if (is.null(colors)) {
    colors <- c(
      "direct" = "#0073C2",
      "indirect" = "#CD534C",
      "total" = "#868686",
      "controlled" = "#EFC000"
    )
  }

  if (type == "effects") {
    # Bar chart of effect estimates
    plot_data <- effects[effects$effect %in% c("cde", "pnde", "tnie", "te"), ]

    # Labels
    plot_data$effect_label <- factor(
      plot_data$effect,
      levels = c("cde", "pnde", "tnie", "te"),
      labels = c("CDE\n(Controlled Direct)",
                 "PNDE\n(Natural Direct)",
                 "TNIE\n(Natural Indirect)",
                 "TE\n(Total Effect)")
    )

    # Assign colors
    plot_data$color_group <- c("controlled", "direct", "indirect", "total")

    if (exponentiate && params$outcome_type != "linear") {
      plot_data$est <- exp(plot_data$est)
      plot_data$lower <- exp(plot_data$lower)
      plot_data$upper <- exp(plot_data$upper)
      ylab <- ifelse(params$outcome_type == "cox", "Hazard Ratio", "Odds Ratio")
      null_val <- 1
    } else {
      ylab <- "Effect Estimate"
      null_val <- 0
    }

    p <- ggplot(plot_data, aes(x = effect_label, y = est, fill = color_group)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = null_val, linetype = "dashed", color = "gray50") +
      scale_fill_manual(values = colors, guide = "none") +
      labs(
        title = if (is.null(title)) "Mediation Effect Estimates" else title,
        x = NULL,
        y = ylab
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(size = 10)
      )

    if (show_ci) {
      p <- p + geom_errorbar(
        aes(ymin = lower, ymax = upper),
        width = 0.2,
        color = "black"
      )
    }

    if (show_pvalue) {
      plot_data$p_label <- ifelse(
        plot_data$p < 0.001, "p<0.001",
        sprintf("p=%.3f", plot_data$p)
      )
      p <- p + geom_text(
        aes(label = p_label),
        vjust = -0.5,
        size = 3
      )
    }

    return(p)

  } else if (type == "decomposition") {
    # Effect decomposition plot
    te_row <- effects[effects$effect == "te", ]
    pnde_row <- effects[effects$effect == "pnde", ]
    tnie_row <- effects[effects$effect == "tnie", ]

    te <- te_row$est
    nde <- pnde_row$est
    nie <- tnie_row$est

    # Calculate proportions (use absolute values for visualization)
    total_abs <- abs(nde) + abs(nie)
    if (total_abs == 0) total_abs <- 1  # Avoid division by zero

    decomp_data <- data.frame(
      effect_type = c("Direct Effect (NDE)", "Indirect Effect (NIE)"),
      value = c(nde, nie),
      abs_value = c(abs(nde), abs(nie)),
      proportion = c(abs(nde) / total_abs, abs(nie) / total_abs),
      color = c("direct", "indirect"),
      stringsAsFactors = FALSE
    )

    # Stacked bar chart (easier to interpret than pie)
    p <- ggplot(decomp_data, aes(x = 1, y = abs_value, fill = effect_type)) +
      geom_col(width = 0.5, position = "stack") +
      geom_text(
        aes(label = sprintf("%.1f%%", proportion * 100)),
        position = position_stack(vjust = 0.5),
        color = "white",
        fontface = "bold",
        size = 5
      ) +
      scale_fill_manual(
        values = c("Direct Effect (NDE)" = colors["direct"],
                   "Indirect Effect (NIE)" = colors["indirect"]),
        name = "Effect Type"
      ) +
      coord_flip() +
      labs(
        title = if (is.null(title)) "Effect Decomposition" else title,
        subtitle = sprintf("Total Effect = %.3f (NDE: %.3f + NIE: %.3f)", te, nde, nie),
        x = NULL,
        y = "Absolute Effect Size"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom"
      )

    return(p)

  } else {
    # Path diagram
    # Create a simplified path diagram using ggplot2
    exposure_name <- params$exposure
    mediator_name <- params$mediator
    outcome_name <- params$outcome

    # Get path coefficients (from mediator and outcome models)
    tnie_row <- effects[effects$effect == "tnie", ]
    pnde_row <- effects[effects$effect == "pnde", ]

    # Create path diagram data
    nodes <- data.frame(
      name = c(exposure_name, mediator_name, outcome_name),
      x = c(0, 1, 2),
      y = c(0, 1, 0),
      stringsAsFactors = FALSE
    )

    # Edges
    edges <- data.frame(
      from = c(exposure_name, mediator_name, exposure_name),
      to = c(mediator_name, outcome_name, outcome_name),
      path_type = c("a path", "b path", "c' path (direct)"),
      stringsAsFactors = FALSE
    )

    p <- ggplot() +
      # Draw edges (arrows)
      geom_segment(
        aes(x = 0.15, y = 0.1, xend = 0.85, yend = 0.9),
        arrow = grid::arrow(length = unit(0.3, "cm")),
        linewidth = 1.2, color = colors["indirect"]
      ) +
      geom_segment(
        aes(x = 1.15, y = 0.9, xend = 1.85, yend = 0.1),
        arrow = grid::arrow(length = unit(0.3, "cm")),
        linewidth = 1.2, color = colors["indirect"]
      ) +
      geom_segment(
        aes(x = 0.15, y = -0.1, xend = 1.85, yend = -0.1),
        arrow = grid::arrow(length = unit(0.3, "cm")),
        linewidth = 1.2, color = colors["direct"]
      ) +
      # Draw nodes
      geom_point(data = nodes, aes(x = x, y = y),
                 size = 20, shape = 21, fill = "white", color = "black", stroke = 2) +
      geom_text(data = nodes, aes(x = x, y = y, label = name),
                size = 4, fontface = "bold") +
      # Labels for paths
      annotate("text", x = 0.5, y = 0.7, label = "a", size = 5, fontface = "italic") +
      annotate("text", x = 1.5, y = 0.7, label = "b", size = 5, fontface = "italic") +
      annotate("text", x = 1, y = -0.3, label = "c' (direct)", size = 4, fontface = "italic") +
      # Effect annotations
      annotate("text", x = 1, y = 1.3,
               label = sprintf("Indirect Effect (NIE): %.3f", tnie_row$est),
               size = 4, color = colors["indirect"]) +
      annotate("text", x = 1, y = -0.5,
               label = sprintf("Direct Effect (NDE): %.3f", pnde_row$est),
               size = 4, color = colors["direct"]) +
      # Theme
      coord_cartesian(xlim = c(-0.5, 2.5), ylim = c(-0.8, 1.5)) +
      labs(
        title = if (is.null(title)) "Mediation Path Diagram" else title
      ) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      )

    return(p)
  }
}


#' Plot Forest Plot for Multiple Mediator Analysis
#'
#' @description
#' Create a forest plot to visualize results from multiple mediator analysis.
#'
#' @param multi_mediation_result A data.frame from \code{run_multi_mediator()}.
#' @param effect_type Character string specifying which effect to display:
#'   "tnie" (indirect effect), "pnde" (direct effect), "te" (total effect),
#'   or "pm" (proportion mediated). Default "tnie".
#' @param exponentiate Logical; whether to exponentiate estimates. Default FALSE.
#' @param null_value Numeric; null effect value for reference line. Default 0.
#' @param title Character string for plot title.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' multi_results <- run_multi_mediator(data, "exposure",
#'                                      c("bmi", "bp", "ldl"), "outcome")
#' plot_mediation_forest(multi_results, effect_type = "tnie")
#' }
#'
#' @import ggplot2
#' @export
plot_mediation_forest <- function(multi_mediation_result,
                                   effect_type = c("tnie", "pnde", "te", "pm"),
                                   exponentiate = FALSE,
                                   null_value = 0,
                                   title = "Mediation Analysis: Multiple Mediators") {

  effect_type <- match.arg(effect_type)

  if (!is.data.frame(multi_mediation_result)) {
    stop("'multi_mediation_result' must be a data.frame from run_multi_mediator().")
  }

  # Select columns based on effect type
  est_col <- effect_type
  se_col <- paste0(effect_type, "_se")
  p_col <- paste0(effect_type, "_p")

  if (!est_col %in% names(multi_mediation_result)) {
    stop(sprintf("Column '%s' not found in results.", est_col))
  }

  # Prepare plot data
  plot_data <- data.frame(
    mediator = multi_mediation_result$mediator,
    estimate = multi_mediation_result[[est_col]],
    se = multi_mediation_result[[se_col]],
    pvalue = multi_mediation_result[[p_col]],
    stringsAsFactors = FALSE
  )

  # Calculate CI
  plot_data$lower <- plot_data$estimate - 1.96 * plot_data$se
  plot_data$upper <- plot_data$estimate + 1.96 * plot_data$se

  # Exponentiate if requested
  if (exponentiate) {
    plot_data$estimate <- exp(plot_data$estimate)
    plot_data$lower <- exp(plot_data$lower)
    plot_data$upper <- exp(plot_data$upper)
    null_value <- 1
  }

  # Remove NA rows
  plot_data <- plot_data[!is.na(plot_data$estimate), ]

  # Order by effect size
  plot_data$mediator <- factor(
    plot_data$mediator,
    levels = plot_data$mediator[order(plot_data$estimate)]
  )

  # Effect type labels
  effect_labels <- c(
    "tnie" = "Indirect Effect (NIE)",
    "pnde" = "Direct Effect (NDE)",
    "te" = "Total Effect",
    "pm" = "Proportion Mediated"
  )

  xlab <- effect_labels[effect_type]
  if (exponentiate && effect_type != "pm") {
    xlab <- paste(xlab, "(exp)")
  }

  # Create forest plot
  p <- ggplot(plot_data, aes(x = estimate, y = mediator)) +
    geom_point(size = 3, shape = 18, color = "#CD534C") +
    geom_errorbarh(
      aes(xmin = lower, xmax = upper),
      height = 0.2,
      linewidth = 0.5,
      color = "#CD534C"
    ) +
    geom_vline(xintercept = null_value, linetype = "dashed", color = "gray50") +
    labs(
      title = title,
      x = xlab,
      y = "Mediator"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )

  # Add significance markers
  plot_data$sig <- ifelse(plot_data$pvalue < 0.05, "*", "")
  plot_data$sig[plot_data$pvalue < 0.01] <- "**"
  plot_data$sig[plot_data$pvalue < 0.001] <- "***"

  p <- p + geom_text(
    data = plot_data,
    aes(x = upper + (max(upper, na.rm = TRUE) - min(lower, na.rm = TRUE)) * 0.05,
        label = sig),
    hjust = 0, size = 5
  )

  return(p)
}


# Multiple Imputation Visualization Functions

#' Plot Multiple Imputation Pooled Results
#'
#' @description
#' Creates a forest plot for pooled estimates from multiple imputation analysis.
#'
#' @param mi_result An object of class \code{mi_pooled_result} from \code{pool_mi_models()}.
#' @param terms Character vector of terms to include. If NULL, all terms except
#'   intercept are shown.
#' @param exponentiate Logical; whether to exponentiate estimates. If NULL,
#'   uses the setting from the mi_result object.
#' @param null_value Numeric; reference line value. If NULL, automatically set
#'   based on exponentiation (0 for linear scale, 1 for exp scale).
#' @param title Character string for plot title.
#' @param colors Named character vector for colors. Default uses package palette.
#' @param show_fmi Logical; whether to display FMI (Fraction of Missing Information)
#'   as point size or annotation. Default TRUE.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' pooled <- pool_mi_models(models, model_type = "cox")
#' plot_mi_pooled(pooled)
#' plot_mi_pooled(pooled, exponentiate = TRUE, show_fmi = TRUE)
#' }
#'
#' @import ggplot2
#' @export
plot_mi_pooled <- function(mi_result,
                            terms = NULL,
                            exponentiate = NULL,
                            null_value = NULL,
                            title = "Pooled Estimates (Multiple Imputation)",
                            colors = NULL,
                            show_fmi = TRUE) {

  if (!inherits(mi_result, "mi_pooled_result")) {
    stop("'mi_result' must be an object from pool_mi_models().")
  }

  # Get pooled data
  plot_data <- mi_result$pooled

  # Determine exponentiation
  if (is.null(exponentiate)) {
    exponentiate <- mi_result$exponentiated
  }

  # Handle exponentiation transformation
  if (exponentiate && !mi_result$exponentiated && mi_result$model_type != "lm") {
    plot_data$estimate <- exp(plot_data$estimate)
    plot_data$conf.low <- exp(plot_data$conf.low)
    plot_data$conf.high <- exp(plot_data$conf.high)
  } else if (!exponentiate && mi_result$exponentiated) {
    plot_data$estimate <- log(plot_data$estimate)
    plot_data$conf.low <- log(plot_data$conf.low)
    plot_data$conf.high <- log(plot_data$conf.high)
  }

  # Filter terms
  if (is.null(terms)) {
    # Exclude intercept by default
    plot_data <- plot_data[!grepl("^\\(Intercept\\)$", plot_data$term), ]
  } else {
    plot_data <- plot_data[plot_data$term %in% terms, ]
  }

  if (nrow(plot_data) == 0) {
    stop("No terms to plot after filtering.")
  }

  # Set null value
  if (is.null(null_value)) {
    null_value <- if (exponentiate) 1 else 0
  }

  # Order by estimate
  plot_data$term <- factor(
    plot_data$term,
    levels = plot_data$term[order(plot_data$estimate)]
  )

  # Default colors
  if (is.null(colors)) {
    colors <- c("point" = "#0073C2", "error" = "#0073C2", "null" = "gray50")
  }

  # Axis label
  xlab <- if (exponentiate && mi_result$model_type != "lm") {
    switch(mi_result$model_type,
           "logistic" = "Odds Ratio",
           "cox" = "Hazard Ratio",
           "poisson" = "Rate Ratio",
           "negbin" = "Rate Ratio",
           "exp(Estimate)")
  } else {
    "Estimate"
  }

  # Create base plot
  if (show_fmi && any(!is.na(plot_data$fmi))) {
    # Use FMI for point size
    plot_data$fmi_pct <- ifelse(is.na(plot_data$fmi), 0, plot_data$fmi * 100)

    p <- ggplot(plot_data, aes(x = estimate, y = term)) +
      geom_point(aes(size = fmi_pct), color = colors["point"], shape = 18) +
      geom_errorbarh(
        aes(xmin = conf.low, xmax = conf.high),
        height = 0.2,
        linewidth = 0.6,
        color = colors["error"]
      ) +
      scale_size_continuous(
        name = "FMI (%)",
        range = c(2, 6),
        guide = guide_legend(title.position = "top")
      )
  } else {
    p <- ggplot(plot_data, aes(x = estimate, y = term)) +
      geom_point(size = 3, color = colors["point"], shape = 18) +
      geom_errorbarh(
        aes(xmin = conf.low, xmax = conf.high),
        height = 0.2,
        linewidth = 0.6,
        color = colors["error"]
      )
  }

  p <- p +
    geom_vline(xintercept = null_value, linetype = "dashed", color = colors["null"]) +
    labs(
      title = title,
      subtitle = sprintf("Model: %s | Imputations: %d",
                         mi_result$model_type, mi_result$n_imputations),
      x = xlab,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )

  # Use log scale for exponentiated estimates
  if (exponentiate && mi_result$model_type != "lm") {
    p <- p + scale_x_log10()
  }

  return(p)
}


#' Plot Multiple Imputation Diagnostics
#'
#' @description
#' Creates diagnostic plots for multiple imputation results, including
#' fraction of missing information (FMI), variance ratios, and degrees of freedom.
#'
#' @param mi_result An object of class \code{mi_pooled_result} from \code{pool_mi_models()}.
#' @param type Character string specifying the diagnostic plot type:
#'   \describe{
#'     \item{"fmi"}{Bar plot of FMI for each coefficient}
#'     \item{"variance_ratio"}{Ratio of between- to within-imputation variance}
#'     \item{"df"}{Degrees of freedom for each coefficient}
#'   }
#' @param title Character string for plot title. If NULL, auto-generated.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' pooled <- pool_mi_models(models, model_type = "logistic")
#' plot_mi_diagnostics(pooled, type = "fmi")
#' }
#'
#' @import ggplot2
#' @export
plot_mi_diagnostics <- function(mi_result,
                                 type = c("fmi", "variance_ratio", "df"),
                                 title = NULL) {

  if (!inherits(mi_result, "mi_pooled_result")) {
    stop("'mi_result' must be an object from pool_mi_models().")
  }

  type <- match.arg(type)
  pooled <- mi_result$pooled

  # Exclude intercept
  pooled <- pooled[!grepl("^\\(Intercept\\)$", pooled$term), ]

  if (type == "fmi") {
    # FMI bar plot
    if (all(is.na(pooled$fmi))) {
      stop("FMI values not available in the result.")
    }

    plot_data <- data.frame(
      term = pooled$term,
      value = pooled$fmi * 100,
      stringsAsFactors = FALSE
    )
    plot_data <- plot_data[!is.na(plot_data$value), ]

    # Order by FMI
    plot_data$term <- factor(
      plot_data$term,
      levels = plot_data$term[order(plot_data$value)]
    )

    # Color based on FMI level
    plot_data$level <- cut(
      plot_data$value,
      breaks = c(-Inf, 10, 30, 50, Inf),
      labels = c("Low (<10%)", "Moderate (10-30%)", "High (30-50%)", "Very High (>50%)")
    )

    if (is.null(title)) title <- "Fraction of Missing Information (FMI)"

    p <- ggplot(plot_data, aes(x = term, y = value, fill = level)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = c(10, 30, 50), linetype = "dashed",
                 color = c("green4", "orange", "red"), alpha = 0.7) +
      scale_fill_manual(
        values = c("Low (<10%)" = "#4DAF4A",
                   "Moderate (10-30%)" = "#FF7F00",
                   "High (30-50%)" = "#E41A1C",
                   "Very High (>50%)" = "#984EA3"),
        name = "FMI Level"
      ) +
      coord_flip() +
      labs(
        title = title,
        x = NULL,
        y = "FMI (%)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )

  } else if (type == "variance_ratio") {
    # Variance ratio: B / U_bar (between / within)
    mi_res <- mi_result$mi_result

    # Calculate variance components
    total_var <- diag(mi_res$variance)
    m <- mi_result$n_imputations

    # FMI ~ (1 + 1/m) * B / T, so B = FMI * T / (1 + 1/m)
    fmi <- pooled$fmi
    between_var <- fmi * total_var / (1 + 1/m)
    within_var <- total_var - (1 + 1/m) * between_var

    # Ratio
    ratio <- between_var / within_var
    ratio[!is.finite(ratio)] <- NA

    plot_data <- data.frame(
      term = pooled$term,
      value = ratio,
      stringsAsFactors = FALSE
    )
    plot_data <- plot_data[!is.na(plot_data$value), ]

    plot_data$term <- factor(
      plot_data$term,
      levels = plot_data$term[order(plot_data$value)]
    )

    if (is.null(title)) title <- "Between/Within Variance Ratio"

    p <- ggplot(plot_data, aes(x = term, y = value)) +
      geom_col(width = 0.7, fill = "#377EB8") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(
        title = title,
        x = NULL,
        y = "Variance Ratio (B/W)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

  } else {
    # Degrees of freedom plot
    plot_data <- data.frame(
      term = pooled$term,
      value = pooled$df,
      stringsAsFactors = FALSE
    )

    # Cap infinite df for display
    plot_data$value[!is.finite(plot_data$value)] <- NA
    plot_data <- plot_data[!is.na(plot_data$value), ]

    if (nrow(plot_data) == 0) {
      stop("No finite degrees of freedom to plot.")
    }

    plot_data$term <- factor(
      plot_data$term,
      levels = plot_data$term[order(plot_data$value)]
    )

    if (is.null(title)) title <- "Degrees of Freedom (Rubin's Adjustment)"

    p <- ggplot(plot_data, aes(x = term, y = value)) +
      geom_col(width = 0.7, fill = "#4DAF4A") +
      coord_flip() +
      labs(
        title = title,
        x = NULL,
        y = "Degrees of Freedom"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  }

  return(p)
}

# Machine Learning Visualization Functions

#' Plot Variable Importance
#'
#' @description
#' Create a bar plot of variable importance from a trained ML model.
#'
#' @param object A ukb_ml object from ukb_ml_model()
#' @param n_features Number of top features to display (default 20)
#' @param type Plot type: "bar" or "dot"
#' @param color Bar color (default "#3182BD")
#' @param title Plot title
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' ml <- ukb_ml_model(outcome ~ ., data, model = "rf")
#' plot_ml_importance(ml, n_features = 15)
#' }
#'
#' @import ggplot2
#' @export
plot_ml_importance <- function(object,
                               n_features = 20,
                               type = c("bar", "dot"),
                               color = "#3182BD",
                               title = "Variable Importance",
                               ...) {
  
  type <- match.arg(type)
  
  # Get importance
  imp <- ukb_ml_importance(object)
  
  if (is.null(imp)) {
    stop("Variable importance not available for this model type")
  }
  
  # Limit features
  if (nrow(imp) > n_features) {
    imp <- imp[seq_len(n_features), ]
  }
  
  # Order factor levels
  imp$variable <- factor(imp$variable, levels = rev(imp$variable))
  
  if (type == "bar") {
    p <- ggplot(imp, aes(x = variable, y = importance)) +
      geom_col(fill = color, width = 0.7) +
      coord_flip() +
      labs(
        title = title,
        x = NULL,
        y = "Importance"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  } else {
    p <- ggplot(imp, aes(x = variable, y = importance)) +
      geom_segment(aes(xend = variable, yend = 0), color = "gray70") +
      geom_point(color = color, size = 3) +
      coord_flip() +
      labs(
        title = title,
        x = NULL,
        y = "Importance"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  }
  
  p
}

#' Plot ROC Curves
#'
#' @description
#' Create ROC curve plot for one or more ML models.
#'
#' @param object A ukb_ml_roc object from ukb_ml_roc()
#' @param ci_alpha Alpha for confidence interval ribbon (default 0.2)
#' @param title Plot title
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' ml <- ukb_ml_model(outcome ~ ., data, model = "rf")
#' roc_result <- ukb_ml_roc(ml)
#' plot(roc_result)
#' }
#'
#' @import ggplot2
#' @export
plot_ml_roc <- function(object,
                        ci_alpha = 0.2,
                        title = "ROC Curve",
                        ...) {
  
  if (!inherits(object, "ukb_ml_roc")) {
    stop("object must be a ukb_ml_roc object")
  }
  
  # Extract ROC data
  roc_data <- data.frame()
  
  for (name in names(object$roc)) {
    roc_obj <- object$roc[[name]]
    auc_val <- object$auc$auc[object$auc$model == name]
    
    df <- data.frame(
      sensitivity = roc_obj$sensitivities,
      specificity = roc_obj$specificities,
      model = sprintf("%s (AUC = %.3f)", name, auc_val),
      stringsAsFactors = FALSE
    )
    
    roc_data <- rbind(roc_data, df)
  }
  
  p <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = model)) +
    geom_line(linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    coord_equal() +
    labs(
      title = title,
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  p
}

#' Plot Calibration Curve
#'
#' @description
#' Create calibration curve plot showing predicted vs observed probabilities.
#'
#' @param object A ukb_ml_calibration object from ukb_ml_calibration()
#' @param title Plot title
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @export
plot_ml_calibration <- function(object,
                                title = "Calibration Curve",
                                ...) {
  
  if (!inherits(object, "ukb_ml_calibration")) {
    stop("object must be a ukb_ml_calibration object")
  }
  
  cal_data <- object$calibration
  
  p <- ggplot(cal_data, aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, color = "gray70") +
    geom_point(aes(size = count), color = "#E34A33") +
    scale_size_continuous(range = c(2, 8), name = "N") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    coord_equal() +
    labs(
      title = title,
      subtitle = sprintf("Brier Score: %.4f, ECE: %.4f", object$brier_score, object$ece),
      x = "Predicted Probability",
      y = "Observed Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  p
}

#' Plot Confusion Matrix
#'
#' @description
#' Create heatmap visualization of confusion matrix.
#'
#' @param object A ukb_ml_confusion object from ukb_ml_confusion()
#' @param normalize Whether to show percentages (default TRUE)
#' @param colors Color gradient (default c("white", "#E34A33"))
#' @param title Plot title
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @export
plot_ml_confusion <- function(object,
                              normalize = TRUE,
                              colors = c("white", "#E34A33"),
                              title = "Confusion Matrix",
                              ...) {
  
  if (!inherits(object, "ukb_ml_confusion")) {
    stop("object must be a ukb_ml_confusion object")
  }
  
  cm <- object$confusion_matrix
  
  # Convert to data frame
  cm_df <- as.data.frame(as.table(cm))
  names(cm_df) <- c("Predicted", "Actual", "Freq")
  
  if (normalize) {
    total <- sum(cm_df$Freq)
    cm_df$Proportion <- cm_df$Freq / total
    cm_df$Label <- sprintf("%d\n(%.1f%%)", cm_df$Freq, cm_df$Proportion * 100)
    fill_var <- "Proportion"
  } else {
    cm_df$Label <- as.character(cm_df$Freq)
    fill_var <- "Freq"
  }
  
  p <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = .data[[fill_var]])) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Label), size = 5) +
    scale_fill_gradient(low = colors[1], high = colors[2]) +
    labs(
      title = title,
      subtitle = sprintf("Accuracy: %.1f%%", object$metrics["accuracy"] * 100),
      x = "Actual",
      y = "Predicted"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  p
}

#' Plot Model Comparison
#'
#' @description
#' Create comparison plot for multiple ML models.
#'
#' @param object A ukb_ml_compare object from ukb_ml_compare()
#' @param metric Metric to highlight (default first available)
#' @param type Plot type: "bar", "dot", or "radar"
#' @param title Plot title
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @export
plot_ml_compare <- function(object,
                            metric = NULL,
                            type = c("bar", "dot"),
                            title = "Model Comparison",
                            ...) {
  
  if (!inherits(object, "ukb_ml_compare")) {
    stop("object must be a ukb_ml_compare object")
  }
  
  type <- match.arg(type)
  comparison <- object$comparison
  
  # Reshape for plotting
  metrics_cols <- setdiff(names(comparison), "model")
  
  plot_data <- data.frame()
  for (m in metrics_cols) {
    df <- data.frame(
      model = comparison$model,
      metric = m,
      value = comparison[[m]],
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, df)
  }
  
  if (type == "bar") {
    p <- ggplot(plot_data, aes(x = model, y = value, fill = metric)) +
      geom_col(position = "dodge", width = 0.7) +
      labs(
        title = title,
        x = NULL,
        y = "Value",
        fill = "Metric"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
  } else {
    p <- ggplot(plot_data, aes(x = model, y = value, color = metric, group = metric)) +
      geom_point(size = 3) +
      geom_line() +
      labs(
        title = title,
        x = NULL,
        y = "Value",
        color = "Metric"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
  }
  
  p
}

#' Plot SHAP Summary
#'
#' @description
#' Create SHAP summary plot (beeswarm or bar) for feature importance.
#'
#' @param object A ukb_shap object from ukb_shap()
#' @param max_features Maximum features to display (default 20)
#' @param type Plot type: "beeswarm" or "bar"
#' @param color_palette Color palette for beeswarm (default "viridis")
#' @param title Plot title
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' ml <- ukb_ml_model(outcome ~ ., data, model = "rf")
#' shap <- ukb_shap(ml)
#' plot_shap_summary(shap)
#' }
#'
#' @import ggplot2
#' @export
plot_shap_summary <- function(object,
                              max_features = 20,
                              type = c("beeswarm", "bar"),
                              color_palette = "viridis",
                              title = "SHAP Summary",
                              ...) {
  
  if (!inherits(object, "ukb_shap")) {
    stop("object must be a ukb_shap object")
  }
  
  type <- match.arg(type)
  
  # Get importance summary
  summary_df <- ukb_shap_summary(object, n = max_features)
  
  if (type == "bar") {
    summary_df$feature <- factor(
      summary_df$feature,
      levels = rev(summary_df$feature)
    )
    
    p <- ggplot(summary_df, aes(x = feature, y = mean_abs_shap)) +
      geom_col(fill = "#3182BD", width = 0.7) +
      coord_flip() +
      labs(
        title = title,
        x = NULL,
        y = "Mean |SHAP Value|"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else {
    # Beeswarm plot
    shap_vals <- object$shap_values
    feature_vals <- object$feature_values
    
    # Get top features
    top_features <- summary_df$feature[seq_len(min(max_features, length(summary_df$feature)))]
    
    # Prepare data for plotting
    plot_data <- data.frame()
    
    for (feat in top_features) {
      idx <- which(object$feature_names == feat)
      
      # Normalize feature values for color
      fv <- as.numeric(feature_vals[[feat]])
      fv_norm <- (fv - min(fv, na.rm = TRUE)) / (max(fv, na.rm = TRUE) - min(fv, na.rm = TRUE) + 1e-10)
      
      df <- data.frame(
        feature = feat,
        shap = shap_vals[, idx],
        value_norm = fv_norm,
        stringsAsFactors = FALSE
      )
      
      plot_data <- rbind(plot_data, df)
    }
    
    # Order features
    plot_data$feature <- factor(
      plot_data$feature,
      levels = rev(top_features)
    )
    
    p <- ggplot(plot_data, aes(x = feature, y = shap, color = value_norm)) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
      scale_color_viridis_c(option = "D", name = "Feature Value\n(normalized)") +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(
        title = title,
        x = NULL,
        y = "SHAP Value"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
  }
  
  p
}

#' Plot SHAP Dependence
#'
#' @description
#' Create SHAP dependence plot showing the relationship between
#' a feature's value and its SHAP value.
#'
#' @param object A ukb_shap object
#' @param feature Feature name to analyze
#' @param color_feature Optional feature for coloring points (interaction)
#' @param alpha Point transparency (default 0.5)
#' @param smooth Add smooth line (default TRUE)
#' @param title Plot title
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @export
plot_shap_dependence <- function(object,
                                 feature,
                                 color_feature = NULL,
                                 alpha = 0.5,
                                 smooth = TRUE,
                                 title = NULL,
                                 ...) {
  
  if (!inherits(object, "ukb_shap")) {
    stop("object must be a ukb_shap object")
  }
  
  # Get dependence data
  dep_data <- ukb_shap_dependence(object, feature, color_feature)
  
  if (is.null(title)) {
    title <- sprintf("SHAP Dependence: %s", feature)
  }
  
  if (is.null(color_feature)) {
    p <- ggplot(dep_data, aes(x = feature_value, y = shap_value)) +
      geom_point(alpha = alpha, color = "#3182BD") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
    
    if (smooth) {
      p <- p + geom_smooth(method = "loess", color = "#E34A33", se = FALSE)
    }
  } else {
    p <- ggplot(dep_data, aes(x = feature_value, y = shap_value, color = color_value)) +
      geom_point(alpha = alpha) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_viridis_c(name = color_feature)
    
    if (smooth) {
      p <- p + geom_smooth(method = "loess", color = "black", se = FALSE)
    }
  }
  
  p <- p +
    labs(
      title = title,
      x = feature,
      y = "SHAP Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  p
}

#' Plot SHAP Force (Waterfall)
#'
#' @description
#' Create a waterfall plot showing feature contributions for a single prediction.
#'
#' @param object A ukb_shap object
#' @param row_id Row index to explain (default 1)
#' @param max_features Maximum features to show (default 10)
#' @param title Plot title
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @export
plot_shap_force <- function(object,
                            row_id = 1,
                            max_features = 10,
                            title = NULL,
                            ...) {
  
  if (!inherits(object, "ukb_shap")) {
    stop("object must be a ukb_shap object")
  }
  
  # Get force data
  force_data <- ukb_shap_force(object, row_id, max_features)
  baseline <- attr(force_data, "baseline")
  prediction <- attr(force_data, "prediction")
  
  if (is.null(title)) {
    title <- sprintf("SHAP Explanation (Observation %d)", row_id)
  }
  
  # Prepare waterfall data
  force_data <- force_data[order(abs(force_data$shap), decreasing = FALSE), ]
  force_data$direction <- ifelse(force_data$shap >= 0, "Positive", "Negative")
  force_data$feature_label <- paste0(force_data$feature, " = ", force_data$value)
  force_data$feature_label <- factor(force_data$feature_label, levels = force_data$feature_label)
  
  p <- ggplot(force_data, aes(x = feature_label, y = shap, fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Positive" = "#E34A33", "Negative" = "#3182BD")) +
    coord_flip() +
    geom_hline(yintercept = 0, color = "black") +
    labs(
      title = title,
      subtitle = sprintf("Baseline: %.3f -> Prediction: %.3f", baseline, prediction),
      x = NULL,
      y = "SHAP Value",
      fill = "Effect"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  p
}

