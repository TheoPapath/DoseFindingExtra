#' Plot Simulated Dose Estimates with Confidence Intervals
#'
#' This function creates a horizontal interval plot showing simulated Effective Doses (ED) or Toxic Doses (TD),
#' with confidence intervals at 95%, 80%, and 50% levels, and overlayed true values.
#'
#' @param x A simulation result object, typically from `planDoseFinding()`, with attributes `altModels` and `doses`.
#' @param type Character, either "ED" or "TD" to specify whether to plot effective or toxic dose estimates.
#' @param p Numeric, the target quantile (used if `type = "ED"`).
#' @param Delta Numeric, the minimal relevant effect size (used if `type = "TD"`).
#' @param xlab Label for the x-axis.
#' @param ci_levels Numeric vector of confidence levels to display. Must be subset of c(0.95, 0.80, 0.50).
#'
#' @return A ggplot object showing interval estimates and reference markers.
#'
#' @examples
#' # Assuming `x` is a simulation object from planDoseFinding()
#' # ggplotDoseSims(x, type = "ED", p = 0.8)
#'
#' @export

ggplotDoseSims <- function(x, type = c("ED", "TD"), p, Delta, xlab = "Dose", ci_levels = c(0.95, 0.80, 0.50)) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DoseFinding)

  altMods <- attr(x, "altModels")
  direction <- attr(altMods, "direction")

  if (type == "ED") {
    out <- DoseFinding:::getSimEst(x, "ED", p = p)
    trueDoses <- DoseFinding:::ED(altMods, p = p, EDtype = "continuous")
  } else {
    out <- DoseFinding:::getSimEst(x, "TD", Delta = Delta, direction = direction)
    trueDoses <- DoseFinding:::TD(altMods, Delta = Delta, TDtype = "continuous", direction = direction)
  }

  nams <- names(out)
  group <- factor(rep(nams, each = length(out[[1]])))
  pdat <- data.frame(est = unlist(out), group = group)

  probs <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
  summary_data <- pdat %>%
    group_by(group) %>%
    summarise(
      q025 = quantile(est, probs[1], na.rm = TRUE),
      q10 = quantile(est, probs[2], na.rm = TRUE),
      q25 = quantile(est, probs[3], na.rm = TRUE),
      q50 = quantile(est, probs[4], na.rm = TRUE),
      q75 = quantile(est, probs[5], na.rm = TRUE),
      q90 = quantile(est, probs[6], na.rm = TRUE),
      q975 = quantile(est, probs[7], na.rm = TRUE),
      na_pct = mean(is.na(est)) * 100
    ) %>%
    ungroup()

  trueDoses_df <- data.frame(group = nams, trueDose = trueDoses)
  summary_data <- left_join(summary_data, trueDoses_df, by = "group")

  maxdose <- max(attr(x, "doses"))
  ylim_range <- c(max(-0.05 * maxdose, min(summary_data$q025) - 0.04),
                  min(2 * maxdose, max(summary_data$q975) + 0.04))

  parVal <- ifelse(type == "ED", paste("p=", p, sep=""), paste("Delta=", Delta, sep=""))
  maintxt <- paste0(
    paste0(ci_levels * 100, "%", collapse = ", "),
    " intervals and median of simulated ", type,
    " estimates (", parVal, ")")

  ggplot(summary_data, aes(y = group)) +
    {if (0.95 %in% ci_levels) geom_linerange(aes(xmin = q025, xmax = q975), size = 5, alpha = 0.2, color = "#1f77b4")} +
    {if (0.80 %in% ci_levels) geom_linerange(aes(xmin = q10, xmax = q90), size = 5, alpha = 0.4, color = "#1f77b4")} +
    {if (0.50 %in% ci_levels) geom_linerange(aes(xmin = q25, xmax = q75), size = 5, alpha = 0.6, color = "#1f77b4")} +
    geom_point(aes(x = trueDose), shape = 18, color = "darkred", size = 5) +
    geom_point(aes(x = q50), shape = 21, fill = "black", color = "black", size = 2) +
    geom_vline(xintercept = c(0, maxdose), linetype = "dotted", color = "grey") +
    labs(title = maintxt,
         y = NULL,
         x = xlab,
         caption = paste0(
           "Horizontal blue bars represent the ",
           paste0(ci_levels * 100, "%", collapse = ", "),
           " confidence intervals (light to dark). ",
           "Black points indicate medians. Red diamonds show true dose values.")) +
    coord_cartesian(xlim = ylim_range) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    if (type == "TD") geom_text(aes(label = sprintf("%% No TD: %.1f%%", na_pct), x = ylim_range[2]),
                                hjust = 1.1, vjust = -0.5, size = 3.5)
}
