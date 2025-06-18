#' Plot Posterior Power Comparison Heatmap
#'
#' This function visualizes pairwise posterior comparison probabilities from the
#' `posterior_power_grand` object returned by `simulate_td_with_power()`.
#'
#' @param posterior_power_grand A data frame with posterior power comparisons including columns:
#'   `model`, `from_dose`, `to_dose`, and the comparison metrics.
#' @param metric Character. Which metric to plot. Must be one of:
#'   - "prob_A_gt_B": Probability that A > B
#'   - "prob_A_gt_B_threshold": Probability that A > B by at least the separation threshold
#'   - "prob_diff_gt_threshold": Probability that |A - B| > separation threshold
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
#' plot_posterior_power_heatmap(results$MaxEff_3$posterior_power_grand, metric = "prob_A_gt_B")
plot_posterior_power_heatmap <- function(
    posterior_power_grand,
    metric = "prob_A_gt_B",
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

  allowed_metrics <- c("prob_A_gt_B", "prob_A_gt_B_threshold", "prob_diff_gt_threshold")
  if (!metric %in% allowed_metrics) {
    warning("Invalid 'metric' specified. Must be one of: ", paste(allowed_metrics, collapse = ", "))
    stop("Metric not supported.")
  }

  required_cols <- c("model", "from_dose", "to_dose", metric)
  missing_cols <- setdiff(required_cols, names(posterior_power_grand))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  if (!is.null(model)) {
    posterior_power_grand <- posterior_power_grand %>% filter(model %in% model)
  }

  caption_text <- switch(
    metric,
    "prob_A_gt_B" = "Posterior probability that the reference dose has greater effect than the comparator dose.",
    "prob_A_gt_B_threshold" = "Posterior probability that the reference dose exceeds the comparator dose by at least the pre-specified threshold.",
    "prob_diff_gt_threshold" = "Posterior probability that the difference between doses exceeds the pre-specified threshold (regardless of direction)."
  )

  label_format <- paste0("%.", digits, "f")

  p <- ggplot(posterior_power_grand, aes(x = factor(to_dose), y = factor(from_dose), fill = .data[[metric]])) +
    geom_tile(color = tile_border_color) +
    geom_text(aes(label = sprintf(label_format, .data[[metric]])), size = text_size) +
    facet_wrap(~ model) +
    scale_fill_gradient2(low = "darkred", high = "darkgreen", mid = "white", midpoint = midpoint) +
    labs(
      x = "Compared To",
      y = "Reference Dose",
      fill = metric,
      title = title %||% paste("Posterior Power Comparison:", metric),
      caption = caption_text
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = legend_position
    )

  return(p)
}
