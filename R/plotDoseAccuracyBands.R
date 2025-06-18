#' Plot Probability of Estimated Dose Accuracy
#'
#' This function visualizes the proportion of simulated estimated doses falling within
#' 30%, 50%, and 100% of the true dose for each model group and Delta/p value.
#'
#' @param x A simulation result object with attributes `altModels` and `doses`.
#' @param type Character, either "ED" or "TD".
#' @param p Optional vector of p values (for ED).
#' @param delta_grid Optional vector of Delta values (for TD).
#'
#' @return A ggplot object showing proportion of accurate estimates.
#'
#' @export

plotDoseAccuracyBands <- function(x, type = c("ED", "TD"), p = NULL, delta_grid = NULL) {
  type <- match.arg(type)
  altMods <- attr(x, "altModels")
  direction <- attr(altMods, "direction")

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
        effect_level = p_val,
        trueDose = rep(trueDoses, each = length(out[[1]]))
      )
    })
  } else {
    if (is.null(delta_grid)) stop("Please provide `delta_grid` for type = 'TD'.")
    sim_data <- purrr::map_dfr(delta_grid, function(Delta) {
      out <- DoseFinding:::getSimEst(x, "TD", Delta = Delta, direction = direction)
      trueDoses <- DoseFinding:::TD(altMods, Delta = Delta, TDtype = "continuous", direction = direction)
      nams <- names(out)
      group <- factor(rep(nams, each = length(out[[1]])))
      data.frame(
        est = unlist(out),
        group = group,
        effect_level = Delta,
        trueDose = rep(trueDoses, each = length(out[[1]]))
      )
    })
  }

  # # This is where I stopped
  # out <- DoseFinding:::getSimEst(x, "dose-response", doseSeq = c(80,160))
  # DoseFinding:::TD(altMods, Delta = Delta, TDtype = "continuous", direction = direction)
  #

  sim_data <- sim_data %>%
    mutate(
      band = case_when(
        abs(est - trueDose)/trueDose <= 0.3 ~ "Within 30%",
        abs(est - trueDose)/trueDose <= 0.5 ~ "Within 50%",
        abs(est - trueDose)/trueDose <= 0.75 ~ "Within 75%",
        TRUE ~ "Above 75%"
      )
    ) %>%
    filter(!is.na(band))

  accuracy_summary <- sim_data %>%
    group_by(group, effect_level, band) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(group, effect_level) %>%
    mutate(probability = n / sum(n)) %>%
    ungroup()

  accuracy_summary$band <- factor(accuracy_summary$band,
                                  levels = c("Above 75%", "Within 75%", "Within 50%", "Within 30%"),
                                  labels = c("Above 75%", "Within 75%", "Within 50%", "Within 30%"))

  ggplot(accuracy_summary, aes(x = effect_level, y = probability, fill = band)) +
    geom_area(position = "stack", color = "white", linewidth = 0.2) +
    facet_wrap(~ group, scales = "free_y") +
    scale_fill_manual(values = c("Within 30%" = "#2ca02c", "Within 50%" = "#ff7f0e", "Within 75%" = "#d62728","Above 75%" = "darkred")) +
    labs(
      x = ifelse(type == "ED", "p (Effect Proportion)", "Delta (Effect Size)"),
      y = "Proportion of Accurate Estimates",
      fill = "Accuracy Band",
      title = "Accuracy of Estimated Doses",
      caption = "Stacked area shows the proportion of simulated dose estimates within 30%, 50%, and 75% of true dose."
    ) +
    theme_minimal(base_size = 14)
}
