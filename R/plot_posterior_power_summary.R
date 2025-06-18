#' Plot Posterior Power Summary
#'
#' Visualize posterior probability that the top dose exceeds placebo by a clinically meaningful threshold,
#' stratified by different `MaxEff` levels (on x-axis) and faceted by model.
#'
#' @param results A named list from `simulate_td_with_power()` with `MaxEff_` keys and associated outputs.
#' @param power_threshold Probability threshold (e.g., 0.8) to display reference line.
#' @param effect_threshold Effect threshold used in power calculation (included for plot annotation).
#' @param metric Character. One of: "prob_top_gt_placebo", "prob_top_gt_placebo_threshold", "prob_top_diff_placebo_threshold".
#' @param digits Integer. Number of digits for display in y-axis.
#' @param show_smooth Logical. If TRUE, adds smoothing line.
#' @param point_size Numeric. Size of plot points.
#' @param line_size Numeric. Size of line geometry.
#' @param legend_position Character. Position of the legend (e.g., "right", "bottom").
#' @param color_palette Optional named vector of colors to apply to models.
#' @param return_data Logical. If TRUE, returns both the plot and the data as a list.
#' @return A ggplot object (or list with plot and data if return_data = TRUE).
#' @export
#'
#' @examples
#' plot_posterior_power_summary(results, metric = "prob_top_gt_placebo")
plot_posterior_power_summary <- function(
    results,
    power_threshold = 0.8,
    effect_threshold = NULL,
    metric = "prob_top_gt_placebo",
    digits = 1,
    show_smooth = FALSE,
    point_size = 2,
    line_size = 0.8,
    legend_position = "right",
    color_palette = NULL,
    return_data = FALSE
) {
  library(dplyr)
  library(ggplot2)

  allowed_metrics <- c("prob_top_gt_placebo", "prob_top_gt_placebo_threshold", "prob_top_diff_placebo_threshold")
  if (!metric %in% allowed_metrics) {
    stop("Invalid metric. Choose one of: ", paste(allowed_metrics, collapse = ", "))
  }

  metric_labels <- c(
    prob_top_gt_placebo = "P(Top dose > Placebo)",
    prob_top_gt_placebo_threshold = paste0("P(Top dose > Placebo + ", effect_threshold, ")"),
    prob_top_diff_placebo_threshold = paste0("P(|Top - Placebo| > ", effect_threshold, ")")
  )

  posterior_power <- purrr::map_dfr(names(results), function(name) {
    data <- results[[name]]$posterior_power
    if (is.null(data)) return(NULL)
    maxEff_val <- suppressWarnings(as.numeric(gsub("MaxEff_", "", name)))
    data$MaxEff <- maxEff_val
    data
  }, .id = "scenario")

  if (nrow(posterior_power) == 0) {
    stop("No posterior power data found in input.")
  }
  if (!(metric %in% names(posterior_power))) {
    stop("Selected metric not found in posterior_power data.")
  }

  p <- ggplot(posterior_power, aes(x = MaxEff, y = .data[[metric]])) +
    geom_line(aes(group = model, color = model), size = line_size) +
    geom_point(aes(color = model), size = point_size) +
    facet_wrap(~ model) +
    geom_hline(yintercept = power_threshold, linetype = "dashed", color = "gray40") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = digits)) +
    labs(
      x = "Max Effect (MaxEff)",
      y = metric_labels[[metric]],
      title = "Posterior Power by Max Effect Level",
      caption = paste0("Metric: ", metric_labels[[metric]])
    ) +
    theme_minimal() +
    theme(legend.position = legend_position)

  if (show_smooth) {
    p <- p + geom_smooth(aes(color = model), method = "loess", se = FALSE, linetype = "dotted")
  }

  if (!is.null(color_palette)) {
    p <- p + scale_color_manual(values = color_palette)
  }

  if (return_data) {
    return(list(plot = p, data = posterior_power))
  } else {
    return(p)
  }
}
