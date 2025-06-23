#' Plot GO/NO GO Outcomes
#'
#' Visualizes GO/NO GO/Consider outcomes from the simulation results
#' using the difference between each selected dose and placebo. This produces
#' a stacked bar plot representing the frequency of each outcome.
#'
#' The function compares predictive quantiles (controlled by `go_prob` and `nogo_prob`)
#' between the active doses and placebo, and categorizes simulation outcomes into:
#'   - GO: both lower and upper quantiles exceed thresholds (go_crit and nogo_crit)
#'   - STOP: both are below their respective thresholds
#'   - Consider: mixed decision boundaries crossed
#'   - Intermediate: reverse mixed boundary condition
#'   - ERROR: fallback if none of the above apply
#'
#' @param results A named list from `simulate_td_with_power()` with MaxEff keys and associated outputs.
#' @param go_crit Numeric. Clinically meaningful threshold for GO decision (e.g., 1.0).
#' @param nogo_crit Numeric. Futility threshold for NO GO decision (e.g., 1.5).
#' @param go_prob Numeric between 0 and 1. Lower quantile probability used for GO decision.
#' @param nogo_prob Numeric between 0 and 1. Upper quantile probability used for STOP decision.
#' @param doses_to_compare Numeric vector of doses to compare against placebo. Default is the maximum dose.
#' @param bar_width Numeric. Width of the bars in the plot.
#' @param fill_colors Named vector of colors for each outcome category.
#' @return A ggplot2 stacked bar plot.
#' @export
#'
#' @examples
#' doses <- c(0, 5, 40, 80, 160)
#' MaxEffs <- c(3)
#' pMods <- setNames(lapply(MaxEffs, function(maxEff) {
#'   DoseFinding::Mods(
#'     sigEmax = c(50, 3.5),
#'     placEff = 4,
#'     maxEff = maxEff,
#'     doses = doses
#'   )
#' }), paste0("MaxEff_", MaxEffs))
#'
#' results <- simulate_td_with_power(
#'   pMods = pMods,
#'   sigma = 0.3,
#'   doses = doses,
#'   delta_grid = c(0.5, 1, 1.5),
#'   go_threshold = 1.5,
#'   separation_threshold = 0.3,
#'   nSim = 10
#' )
#' plot_go_nogo_outcomes(results, go_crit = 1.0, nogo_crit = 1.5)
plot_go_nogo_outcomes <- function(
    results,
    go_crit = 3,
    nogo_crit = 2,
    go_prob = 0.2,
    nogo_prob = 0.9,
    doses_to_compare = NULL,
    bar_width = 0.6,
    fill_colors = c(
      GO = "darkgreen",
      STOP = "darkred",
      Consider = "gold",
      Intermediate = "steelblue",
      ERROR = "gray"
    )
) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)

  go_q <- paste0("q", round(go_prob * 100))
  nogo_q <- paste0("q", round(nogo_prob * 100))

  # Collapse all prediction quantiles into long format
  predictions_df <- purrr::imap_dfr(results, function(res, maxeff_name) {
    pred_list <- res$predictions_df_var$data
    if (is.null(pred_list)) return(NULL)

    purrr::imap_dfr(pred_list, function(df, model_name) {
      if (nrow(df) == 0) return(NULL)
      df %>%
        mutate(quantile = stringr::str_replace(quantile, "(\\d+)%", "q\\1")) |>
        filter(quantile %in% c(go_q, nogo_q)) %>%
        pivot_wider(names_from = quantile, values_from = value) %>%
        mutate(
          model = model_name,
          MaxEff = gsub("MaxEff_", "", maxeff_name)
        )
    })
  })

  if (nrow(predictions_df) == 0) {
    stop("No simulation quantile data found in input.")
  }

  placebo_val <- min(predictions_df$dose, na.rm = TRUE)
  all_doses <- unique(predictions_df$dose)

  if (is.null(doses_to_compare)) {
    doses_to_compare <- max(all_doses, na.rm = TRUE)
  }

  invalid_doses <- setdiff(doses_to_compare, all_doses)
  if (length(invalid_doses) > 0) {
    stop("Invalid doses provided in doses_to_compare: ", paste(invalid_doses, collapse = ", "))
  }

  # Compute GO/STOP outcomes for selected doses vs placebo
  placebo_df <- predictions_df %>%
    filter(dose == placebo_val) %>%
    rename(pred20_placebo = !!go_q, pred90_placebo = !!nogo_q) %>%
    select(simulation, model, MaxEff, pred20_placebo, pred90_placebo)

  outcome_df <- predictions_df %>%
    filter(dose %in% doses_to_compare) %>%
    left_join(placebo_df, by = c("simulation", "model", "MaxEff")) %>%
    mutate(
      diff20 = !!sym(go_q) - pred20_placebo,
      diff90 = !!sym(nogo_q) - pred90_placebo,
      outcome = case_when(
        diff20 > go_crit & diff90 > nogo_crit ~ "GO",
        diff20 < go_crit & diff90 < nogo_crit ~ "STOP",
        diff20 < go_crit & diff90 > nogo_crit ~ "Consider",
        diff20 > go_crit & diff90 < nogo_crit ~ "Intermediate",
        TRUE ~ "ERROR"
      )
    )

  # Count occurrences per outcome for stacked bar plot
  summary_df <- outcome_df %>%
    mutate(outcome = factor(outcome, levels = c("STOP", "Consider", "Intermediate", "GO", "ERROR"))) %>%
    count(model, MaxEff, dose, outcome, name = "count") %>%
    group_by(model, MaxEff, dose) %>%
    mutate(pct = count / sum(count)) %>%
    ungroup()

  ggplot(summary_df, aes(x = factor(MaxEff), y = pct, fill = outcome)) +
    geom_bar(stat = "identity", position = "stack", width = bar_width) +
    geom_text(aes(label = scales::percent(pct, accuracy = 1)),
              position = position_stack(vjust = 0.5), size = 3, color = "black") +
    facet_wrap(~ paste0("Dose: ", dose), scales = "free_y") +
    scale_fill_manual(values = fill_colors, name = "Decision") +
    labs(
      x = "Maximum Possible Effect",
      y = "Probability(%)",
      title = "GO/NO GO Outcome Probabilities vs Placebo",
      caption = "Each bar represents the proportion of simulations
      classified under each decision category based on predictive quantiles at the evaluated dose."
    ) +
    geom_hline(yintercept = 0.8) +
    scale_y_continuous(labels = scales::percent_format(), breaks = seq(0,1,.2)) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 16)  # Set title font size here
    )
}
