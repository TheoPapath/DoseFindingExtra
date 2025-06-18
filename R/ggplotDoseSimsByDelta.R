#' Plot Simulated Dose Estimates Across a Range of Delta Values
#'
#' This function visualizes simulated ED or TD estimates across multiple Delta values
#' and model groups. Each facet represents a model group. Within each facet,
#' the x-axis is dose, the y-axis is Delta, and confidence intervals and medians are shown.
#'
#' @param x A simulation result object with attributes `altModels` and `doses`.
#' @param delta_grid A numeric vector of Delta values to evaluate (e.g., seq(0.1, 1.0, by = 0.1)).
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param ci_levels Numeric vector of confidence levels to display.
#'
#' @return A ggplot object with facets per model group and Delta.
#'
#' @examples
#' # Assuming `x` is a simulation object
#' # ggplotDoseSimsByDelta(x, delta_grid = seq(0.1, 1.5, by = 0.1))
#'
#' @export

# Work in progress

ggplotDoseSimsByDelta <- function(x, type = c("ED", "TD"), p = NULL, delta_grid = NULL, xlab = "Dose", ylab = "Delta", ci_levels = c(0.95, 0.80, 0.50)) {
  type <- match.arg(type)
  altMods <- attr(x, "altModels")
  direction <- attr(altMods, "direction")
  doses <- attr(x, "doses")

  # Prepare a dataframe based on type
  if (type == "ED") {
    if (is.null(p)) stop("Please provide `p` for type = 'ED'.")
    sim_data <- purrr::map_dfr(p, function(p_val) {
      out <- DoseFinding:::getSimEst(x, "ED", p = p_val)
      trueDoses <- DoseFinding:::ED(altMods, p = p_val, EDtype = "continuous")
      nams <- names(out)
      group <- factor(rep(nams, each = length(out[[1]])))
      data.frame(
        est = unlist(out),
        group = group,
        Delta = p_val,
        trueDose = rep(trueDoses, each = length(out[[1]]))
      )
    })
  } else {
    if (is.null(delta_grid)) stop("Please provide `delta_grid` for type = 'TD'.")

  }
  sim_data <- purrr::map_dfr(delta_grid, function(Delta) {
    out <- DoseFinding:::getSimEst(x, "TD", Delta = Delta, direction = direction)
    trueDoses <- DoseFinding:::TD(altMods, Delta = Delta, TDtype = "continuous", direction = direction)
    nams <- names(out)
    group <- factor(rep(nams, each = length(out[[1]])))
    data.frame(
      est = unlist(out),
      group = group,
      Delta = Delta,
      trueDose = rep(trueDoses, each = length(out[[1]]))
    )
  })

  # Summarise by group and Delta
  summary_data <- sim_data %>%
    group_by(group, Delta) %>%
    summarise(
      q025 = quantile(est, 0.025, na.rm = TRUE),
      q10 = quantile(est, 0.10, na.rm = TRUE),
      q25 = quantile(est, 0.25, na.rm = TRUE),
      q50 = quantile(est, 0.50, na.rm = TRUE),
      q75 = quantile(est, 0.75, na.rm = TRUE),
      q90 = quantile(est, 0.90, na.rm = TRUE),
      q975 = quantile(est, 0.975, na.rm = TRUE),
      trueDose = unique(trueDose),
      na_pct = mean(is.na(est)) * 100,
      .groups = "drop"
    )

  maxdose <- max(attr(x, "doses"))
  ylim_range <- c(max(-0.05 * maxdose, min(summary_data$q025) - 0.04),
                  min(2 * maxdose, max(summary_data$q975) + 0.04))

  ggplot(summary_data, aes(x = trueDose, y = Delta, group = group)) +
    {if (0.95 %in% ci_levels) geom_ribbon(aes(xmin = q025, xmax = q975), height = 0.15, alpha = 0.2, fill = "#1f77b4")} +
    {if (0.80 %in% ci_levels) geom_ribbon(aes(xmin = q10, xmax = q90), height = 0.15, alpha = 0.4, fill = "#1f77b4")} +
    {if (0.50 %in% ci_levels) geom_ribbon(aes(xmin = q25, xmax = q75), height = 0.15, alpha = 0.6, fill = "#1f77b4")} +

    geom_line(linewidth = 0.8, color = "black") +
    geom_point(aes(x = q50), shape = 21, fill = "black", size = 2) +
    geom_line(aes(x = trueDose, y = Delta, group = group), color = "red", linewidth = 0.8) +
    geom_point(aes(x = trueDose, y = Delta), shape = 18, color = "red", size = 2) +
    facet_wrap(~ group) +
    labs(x = xlab, y = ylab,
         caption = paste0("Black lines connect median dose estimates across Delta. Shaded bands represent ",
                          paste0(ci_levels * 100, "%", collapse = ", "),
                          " confidence intervals. Red diamonds = true values.")) +
    coord_flip(xlim = ylim_range) +
    # Depends on how I want to see the data. I think now make more sense
    # coord_cartesian(ylim = ylim_range) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) # + coord_flip()

}
