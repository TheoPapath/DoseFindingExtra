#' Plot Posterior Power Summary
#'
#' Visualize posterior probability that the top dose exceeds placebo by a clinically meaningful threshold,
#' stratified by different `MaxEff` levels (on x-axis) and faceted by model.
#'
#' @param results A named list from `simulate_td_with_power()` with `MaxEff_` keys and associated outputs.
#' @param power_threshold Probability threshold (e.g., 0.8) to display reference line.
#' @param effect_threshold Effect threshold used in power calculation (included for plot annotation). If NULL, will extract from input_arguments in `results`.
#' @param metric Character. One of: "prob_top_gt_placebo", "prob_top_gt_placebo_threshold", "prob_top_diff_placebo_threshold", or "all".
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
    go_threshold = NULL,
    separation_threshold = NULL,
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

  allowed_metrics <- c("prob_top_gt_placebo", "prob_top_gt_placebo_threshold", "prob_top_diff_placebo_threshold", "all")
  if (!metric %in% allowed_metrics) {
    stop("Invalid metric. Choose one of: ", paste(allowed_metrics, collapse = ", "))
  }

  if (is.null(go_threshold)) {
    first_args <- results[[1]]$input_arguments
    go_threshold <- if (!is.null(first_args$go_threshold)) first_args$go_threshold else NA
  }

  if (is.null(separation_threshold)) {
    first_args <- results[[1]]$input_arguments
    separation_threshold <- if (!is.null(first_args$separation_threshold)) first_args$separation_threshold else NA
  }


  metric_labels <- c(
    prob_top_gt_placebo = "P(Top dose > Placebo)",
    prob_top_gt_placebo_threshold = paste0("P(Top dose > Placebo + ", go_threshold, ")"),
    prob_top_diff_placebo_threshold = paste0("P(|Top - Placebo| > ", separation_threshold, ")")
  )

  posterior_power <- purrr::map_dfr(names(results), function(name) {
    data <- results[[name]]$posterior_probabilities$data
    if (is.null(data)) return(NULL)
    maxEff_val <- suppressWarnings(as.numeric(gsub(".*?([0-9.]+)$", "\\1", name)))
    data$MaxEff <- maxEff_val
    data
  }, .id = "scenario")

  if (nrow(posterior_power) == 0) {
    stop("No posterior power data found in input.")
  }

  if (metric == "all") {
    posterior_power_long <- posterior_power |>
      tidyr::pivot_longer(
        cols = c("prob_top_gt_placebo", "prob_top_gt_placebo_threshold", "prob_top_diff_placebo_threshold"),
        names_to = "metric", values_to = "value"
      ) |>
      mutate(metric_label = metric_labels[metric])

    p <- ggplot(posterior_power_long, aes(x = MaxEff, y = value)) +
      geom_line(aes(group = model, color = model), size = line_size) +
      geom_point(aes(color = model), size = point_size) +
      # facet_grid(metric_label ~ model, scales = "free_y") +
      facet_wrap( ~ metric_label, scales = "free_y", ncol = 2) +
      geom_hline(yintercept = power_threshold, linetype = "dashed", color = "gray40") +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_y_continuous(labels = scales::percent_format(accuracy = digits)) +
      labs(
        x = "Max Effect (MaxEff)",
        y = "Posterior Probability",
        title = "Posterior Power by Max Effect Level",
        color = "Model"
      ) +
      theme_minimal() +
      theme(legend.position = legend_position)
  } else {
    if (!(metric %in% names(posterior_power))) {
      stop("Selected metric not found in posterior power data.")
    }

    p <- ggplot(posterior_power, aes(x = MaxEff, y = .data[[metric]])) +
      geom_line(aes(group = model, color = model), size = line_size) +
      geom_point(aes(color = model), size = point_size) +
      # facet_wrap(~ model) +
      geom_hline(yintercept = power_threshold, linetype = "dashed", color = "gray40") +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_y_continuous(labels = scales::percent_format(accuracy = digits)) +
      labs(
        x = "Max Effect (MaxEff)",
        y = metric_labels[[metric]],
        title = "Posterior Power by Max Effect Level",
        caption = paste0("Metric: ", metric_labels[[metric]]),
        color = "Model"
      ) +
      theme_minimal() +
      theme(legend.position = legend_position)
  }

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
