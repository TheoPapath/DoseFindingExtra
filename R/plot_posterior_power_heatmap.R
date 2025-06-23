#' Plot Posterior Power Comparison Heatmap
#'
#' This function visualizes pairwise posterior comparison probabilities from the
#' `simulate_td_with_power()` output.
#'
#' @param results A list object returned by `simulate_td_with_power()`.
#' @param maxeff The value of the MaxEff scenario to visualize (e.g., maxeff=3 for "MaxEff_3").
#' @param metric Character. Which metric to plot. Must be one of:
#'   - "prob_top_gt_placebo_threshold": Probability that top dose > placebo + threshold
#'   - "prob_top_diff_placebo_threshold": Probability that |top - placebo| > threshold
#'   - "all": Plot both metrics in facets.
#' @param midpoint Midpoint value for the heatmap gradient. Default is 0.5.
#' @param title Optional title for the plot.
#' @param digits Number of decimal places shown in labels. Default is 2.
#' @param model Optional character vector. Subset to specific model(s).
#' @param text_size Size of the numbers in each tile. Default is 3.
#' @param tile_border_color Color for the tile border. Default is "white".
#' @param legend_position Legend position. Default is "right".
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' plot_posterior_power_heatmap(results, maxeff = 3, metric = "prob_top_gt_placebo_threshold")
plot_posterior_power_heatmap <- function(
    results,
    maxeff,
    metric = "prob_top_gt_placebo_threshold",
    midpoint = 0.5,
    title = NULL,
    digits = 2,
    model = NULL,
    text_size = 3,
    tile_border_color = "white",
    legend_position = "right"
) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  maxeff_name <- paste0("MaxEff_", maxeff)

  if (!maxeff_name %in% names(results)) {
    stop("Invalid maxeff. Available options are: ", paste(names(results), collapse = ", "))
  }

  entry <- results[[maxeff_name]]
  input_args <- entry$input_arguments
  posterior_power_grand <- entry$pairwise_comparisons$data
  attr(posterior_power_grand, "input_arguments") <- input_args

  allowed_metrics <- c("prob_top_gt_placebo_threshold", "prob_top_diff_placebo_threshold", "all")
  if (!metric %in% allowed_metrics) {
    stop("Invalid 'metric' specified. Must be one of: ", paste(allowed_metrics, collapse = ", "))
  }

  if (!is.null(model)) {
    posterior_power_grand <- posterior_power_grand %>% filter(model %in% model)
  }

  go_threshold <- input_args$go_threshold %||% ""
  separation_threshold <- input_args$separation_threshold %||% ""

  label_format <- paste0("%.", digits, "f")

  metric_labels <- c(
    prob_top_gt_placebo_threshold = paste0("P(Top dose > Placebo + ", go_threshold, ")"),
    prob_top_diff_placebo_threshold = paste0("P(|Top - Placebo| > ", separation_threshold, ")")
  )

  if (metric == "all") {
    df_long <- posterior_power_grand %>%
      pivot_longer(cols = c("prob_top_gt_placebo_threshold", "prob_top_diff_placebo_threshold"),
                   names_to = "metric", values_to = "value") %>%
      filter(from_dose-to_dose >= 0) %>%
      mutate(metric_label = metric_labels[metric])

    caption_text <- paste0("Posterior probabilities for dose comparisons (upper triangle only).\n",
                           "Thresholds: GO = ", go_threshold, ", Separation = ", separation_threshold)

    p <- ggplot(df_long, aes(x = factor(to_dose), y = factor(from_dose), fill = value)) +
      geom_tile(color = tile_border_color) +
      geom_text(aes(label = sprintf(label_format, value)), size = text_size) +
      facet_grid(metric_label ~ model, scales = "free") +
      scale_fill_gradientn(colors = c("darkred", "white", "darkgreen"), limits = c(0, 1)) +
      labs(
        x = "Compared To",
        y = "Reference Dose",
        fill = "Posterior Probability",
        title = title %||% "Pairwise Posterior Comparison (All Metrics)",
        caption = caption_text
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = legend_position,
        strip.text = element_text(size = 10)
      ) +
      coord_flip()
  } else {
    if (!(metric %in% names(posterior_power_grand))) {
      stop("Selected metric not found in pairwise comparisons data.")
    }

    posterior_power_grand2 <- posterior_power_grand |> dplyr::filter(from_dose-to_dose >= 0)
    caption_text <- switch(
      metric,
      "prob_top_gt_placebo_threshold" = paste0("Posterior probability that top dose exceeds placebo + threshold (", go_threshold, ")."),
      "prob_top_diff_placebo_threshold" = paste0("Posterior probability that |top - placebo| > threshold (", separation_threshold, ").")
    )

    p <- ggplot(posterior_power_grand2, aes(x = factor(to_dose), y = factor(from_dose), fill = .data[[metric]])) +
      geom_tile(color = tile_border_color) +
      geom_text(aes(label = sprintf(label_format, .data[[metric]])), size = text_size) +
      facet_wrap(~ model) +
      scale_fill_gradientn(colors = c("darkred","white","darkgreen"), limits = c(0, 1)) +
      labs(
        x = "Compared To",
        y = "Reference Dose",
        fill = metric_labels[[metric]],
        title = title %||% paste("Posterior Power Comparison:", metric_labels[[metric]], "-", maxeff_name),
        caption = caption_text
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = legend_position,
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
      coord_flip()
  }
  return(p)
}
