#' Compute Recommended Group Separation for Signal-Seeking Studies
#'
#' This function computes the recommended effect size for a second group
#' given a reference effect, standard deviation, and a desired effect size category
#' (small, moderate, large, or very large), based on Cohen's d thresholds.
#'
#' Cohen's d is a standardized measure of effect size:
#'   d = (Delta_2 - Delta_1) / SD
#'
#' Typical Cohen's d thresholds are:
#' - Small effect:    d = 0.2
#' - Moderate effect: d = 0.5
#' - Large effect:    d = 0.8
#' - Very large:      d = 1.0+
#'
#' In signal-seeking settings (e.g., early-phase studies, biomarker exploration),
#' moderate to large effect sizes (d = 0.5 to 0.8) are often recommended to ensure
#' a detectable signal, particularly under noisy conditions.
#'
#' @param sd Numeric. The standard deviation of the outcome.
#' @param delta_ref Numeric. The expected mean effect level for the reference group (e.g., Group 1).
#' @param effect_size Character. One of: "small", "moderate", "large", "very_large".
#'
#' @return A list with:
#' \describe{
#'   \item{delta_ref}{The reference group effect (input).}
#'   \item{sd}{The standard deviation (input).}
#'   \item{cohen_d}{The Cohen's d threshold used.}
#'   \item{delta_target}{The recommended target group effect to achieve the specified effect size.}
#'   \item{delta_difference}{The required difference between groups.}
#' }
#'
#' @examples
#' calculate_group_separation_thresholds(sd = 3.3, delta_ref = 1, effect_size = "moderate")
#' calculate_group_separation_thresholds(sd = 3.3, delta_ref = 1, effect_size = "large")
#'
#' @export
calculate_group_separation_thresholds <- function(sd, delta_ref, effect_size = c("small", "moderate", "large", "very_large")) {
  effect_size <- match.arg(effect_size)

  # Cohen's d thresholds
  d_thresholds <- c(
    small = 0.2,
    moderate = 0.5,
    large = 0.8,
    very_large = 1.0
  )

  d_value <- d_thresholds[[effect_size]]

  # Calculate minimum required separation
  delta_difference <- d_value * sd
  delta_target <- delta_ref + delta_difference

  result <- list(
    delta_ref = delta_ref,
    sd = sd,
    cohen_d = d_value,
    delta_difference = delta_difference,
    delta_target = delta_target
  )

  class(result) <- "group_separation_recommendation"
  return(result)
}

# Custom print method for clean display
print.group_separation_recommendation <- function(x, ...) {
  cat("Group Separation Recommendation:\n")
  cat("  Reference effect (Delta_1):", x$delta_ref, "\n")
  cat("  Standard deviation (SD):   ", x$sd, "\n")
  cat("  Cohen's d threshold:       ", x$cohen_d, "\n")
  cat("  Required separation (Δ₂ - Δ₁):", round(x$delta_difference, 2), "\n")
  cat("  Target effect for Group 2 (Delta_2):", round(x$delta_target, 2), "\n")
}
