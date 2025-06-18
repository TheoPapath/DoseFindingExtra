#' Simulate Target Dose Estimation and Power Across Efficacy Scenarios
#'
#' This function runs simulation-based target dose estimation and posterior probability analysis
#' under multiple efficacy scenarios defined by varying MaxEff values in a `Mods` object.
#'
#' @param pMods Named list of `Mods` objects with different max effect levels.
#' @param sigma Standard deviation of the outcome.
#' @param doses Numeric vector of dose levels used in the study.
#' @param delta_grid Vector of Delta values for which target dose estimation will be evaluated.
#' @param models_to_estimate Character vector of model types to fit (default includes common parametric models).
#' @param n Sample size. Can be a vector. When a single number is specified for ‘⁠n⁠’ it is assumed this is the sample size per group and balanced allocations are used.
#' @param go_threshold Minimum meaningful effect (Delta) required to support a "go" decision.
#' @param separation_threshold Minimum difference used for pairwise posterior probability evaluations.
#' @param nSim Number of simulation replicates.
#'
#' @return A list where each entry corresponds to a `MaxEff` scenario and contains:
#'   - `summary_data`: quantiles of estimated target doses
#'   - `posterior_power`: summary of power for detecting meaningful differences vs placebo
#'   - `posterior_power_grand`: matrix of pairwise posterior comparisons across doses
#'   - `td_summary_formatted`: formatted table of true vs estimated target doses with CIs
#'
#' @importFrom DoseFinding planMod defBnds Mods
#' @importFrom purrr map_dfr map2_dbl
#' @importFrom dplyr summarise transmute group_by mutate select filter
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' doses <- c(0, 5, 40, 80, 160)
#' MaxEff_vals <- c(3, 4)
#' pMods <- setNames(lapply(MaxEff_vals, function(maxEff) {
#'   Mods(
#'     emax = 80,
#'     sigEmax = c(50, 3.5),
#'     quadratic = -0.004,
#'     placEff = 4,
#'     maxEff = maxEff,
#'     doses = doses
#'   )
#' }), paste0("MaxEff_", MaxEff_vals))
#'
#' results <- simulate_td_with_power(
#'   pMods = pMods,
#'   sigma = 5,
#'   doses = doses,
#'   delta_grid = c(0.5, 1, 1.5),
#'   go_threshold = 1.5,
#'   separation_threshold = 0.3,
#'   nSim = 100
#' )
#' results$MaxEff_3$td_summary_formatted
simulate_td_with_power <- function(
    pMods,  # <-- NEW: pass Mods() object as argument
    sigma,
    doses,
    delta_grid,
    models_to_estimate = c("betaMod", "emax", "sigEmax", "linlog"),
    n = 100,
    go_threshold = 1.0,
    separation_threshold = 0.3,
    nSim = 100
) {
  results <- list()
  for (MaxEff_name in names(pMods)) {
    mods <- pMods[[MaxEff_name]]
    # total_n <- n_per_arm * length(ratio)
    # n <- total_n * ratio / sum(ratio)
    pObj <- planMod(
      models_to_estimate, mods, n, sigma, doses = doses,
      simulation = TRUE, nSim = nSim, asyApprox = FALSE,
      bnds = defBnds(max(doses))
    )
    altMods <- attr(pObj, "altModels")
    direction <- attr(altMods, "direction")
    doses <- attr(pObj, "doses")
    sim_data <- purrr::map_dfr(delta_grid, function(Delta) {
      out <- DoseFinding:::getSimEst(pObj, "TD", Delta = Delta, direction = direction)
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
    dose_response_profiles <- DoseFinding:::getSimEst(pObj, type = "dose-response", doseSeq = doses)
    posterior_power <- map_dfr(dose_response_profiles, function(mat) {
      dose_names <- as.numeric(colnames(mat))
      placebo_col <- which.min(dose_names)
      top_col <- which.max(dose_names)
      diff <- mat[, top_col] - mat[, placebo_col]
      prob_gt_0 <- mean(diff > 0)
      prob_gt_threshold <- mean(diff > go_threshold)
      prob_diff <- mean(abs(diff) > separation_threshold)
      tibble(
        prob_top_gt_placebo = prob_gt_0,
        prob_top_gt_placebo_threshold = prob_gt_threshold,
        prob_top_diff_placebo_threshold = prob_diff
      )
    }, .id = "model")
    posterior_power_grand <- map_dfr(dose_response_profiles, function(mat) {
      dose_names <- as.numeric(colnames(mat))
      comparisons <- expand.grid(
        comp_from = seq_along(dose_names),
        comp_to = seq_along(dose_names)
      ) %>%
        filter(comp_from != comp_to) %>%
        mutate(
          from_dose = dose_names[comp_from],
          to_dose = dose_names[comp_to],
          contrast = paste0(from_dose, " vs ", to_dose),
          prob_A_gt_B = map2_dbl(comp_from, comp_to, ~ mean(mat[, .x] > mat[, .y])),
          prob_A_gt_B_threshold = map2_dbl(comp_from, comp_to, ~ mean((mat[, .x] - mat[, .y]) > separation_threshold)),
          prob_diff_gt_threshold = map2_dbl(comp_from, comp_to, ~ mean(abs(mat[, .x] - mat[, .y]) > separation_threshold))
        )
      comparisons %>%
        select(contrast, from_dose, to_dose, prob_A_gt_B, prob_A_gt_B_threshold, prob_diff_gt_threshold)
    }, .id = "model")
    est_td_data <- DoseFinding:::getSimEst(pObj, "TD", Delta = go_threshold, direction = direction)
    true_td <- DoseFinding:::TD(altMods, Delta = go_threshold, TDtype = "continuous", direction = direction)
    td_summary <- map_dfr(est_td_data, function(x) {
      tibble(
        q025 = quantile(x, 0.025, na.rm = TRUE),
        q05  = quantile(x, 0.05, na.rm = TRUE),
        q10  = quantile(x, 0.10, na.rm = TRUE),
        q50  = quantile(x, 0.50, na.rm = TRUE),
        q90  = quantile(x, 0.90, na.rm = TRUE),
        q95  = quantile(x, 0.95, na.rm = TRUE),
        q975 = quantile(x, 0.975, na.rm = TRUE),
        mean = mean(x, na.rm = TRUE)
      )
    }, .id = "model") %>%
      mutate(true_target_dose = true_td[model])
    td_summary_formatted <- td_summary %>%
      transmute(
        Model = model,
        MaxEff = MaxEff_name,
        `True Target Dose` = true_target_dose,
        `Mean Estimated Dose (95% CI)` = sprintf("%.1f (%.1f–%.1f)", mean, q025, q975),
        `Mean Estimated Dose (90% CI)` = sprintf("%.1f (%.1f–%.1f)", mean, q05, q95),
        `Mean Estimated Dose (80% CI)` = sprintf("%.1f (%.1f–%.1f)", mean, q10, q90)
      )
    results[[MaxEff_name]] <- list(
      summary_data = summary_data,
      posterior_power = posterior_power,
      posterior_power_grand = posterior_power_grand,
      td_summary_formatted = td_summary_formatted
    )
  }
  return(results)
}
