#' Simulate Target Dose Estimation and Power Across Efficacy Scenarios
#'
#' This function runs simulation-based target dose estimation and posterior probability analysis
#' under multiple efficacy scenarios defined by varying MaxEff values in a `Mods` object.
#'
#' @param pMods A single `Mods` object or a named list of `Mods` objects with different max effect levels.
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
#'   - `target_dose_estimation`: list with `data` and a `description`
#'   - `posterior_probabilities`: list with `data` and a `description`
#'   - `pairwise_comparisons`: list with `data` and a `description`
#'   - `target_dose_summary`: list with `table` and a `description`
#'   - `planMod_output`: raw output from planMod()
#'   - `predictions_df`: matrix of simulated dose-response predictions
#'   - `dose_response_quantiles_long`: predictive quantile summaries per simulation
#'
#' @importFrom DoseFinding planMod defBnds Mods
#' @importFrom purrr map_dfr map2_dbl
#' @importFrom dplyr summarise transmute group_by mutate select filter
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' doses <- c(0, 5, 40, 80, 160)
#' MaxEff_vals <- c(3, 4)
#' pMods <- setNames(lapply(MaxEff_vals, function(maxEff) {
#'   DoseFinding::Mods(
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
#'   sigma = 0.3,
#'   doses = doses,
#'   delta_grid = c(0.5, 1, 1.5),
#'   go_threshold = 1.5,
#'   separation_threshold = 0.3,
#'   nSim = 100
#' )
#' results$MaxEff_3$target_dose_summary$table
simulate_td_with_power <- function(
    pMods,
    sigma,
    doses,
    delta_grid,
    models_to_estimate = c("emax", "sigEmax"),
    n = 100,
    go_threshold = 1.0,
    separation_threshold = 0.3,
    nSim = 100
) {
  results <- list()

  if (inherits(pMods, "Mods")) {
    pMods <- list(SingleScenario = pMods)
  }

  for (MaxEff_name in names(pMods)) {
    # MaxEff_name = "MaxEff_3"
    mods <- pMods[[MaxEff_name]]

    pObj <- DoseFinding::planMod(
      models_to_estimate, mods, n, sigma, doses = doses,
      simulation = TRUE, nSim = nSim, asyApprox = FALSE,
      bnds = DoseFinding::defBnds(max(doses))
    )
    altMods <- attr(pObj, "altModels")
    direction <- attr(altMods, "direction")
    doses <- attr(pObj, "doses")

    # this gives a dataframe of the ground truth models
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

    summary_data <- sim_data |>
      dplyr::group_by(group, Delta) |>
      dplyr::summarise(
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

    # -------------------------------------------------------------------------------------------------
    # Predictions module (point wise with and variance)
    # -------------------------------------------------------------------------------------------------
    dose_response_profiles <- DoseFinding:::getSimEst(pObj, type = "dose-response", doseSeq = doses)

    # Add selected model info to dose_response_profiles
    modelSel <- attr(pObj$sim, "modelSel") |> as.data.frame()
    for (mod in names(dose_response_profiles)) {
      modelSelc <- modelSel |> data.frame() |> dplyr::select(mod)
      if (nrow(modelSelc) == nrow(dose_response_profiles[[mod]])) {
        dose_response_profiles[[mod]] <- dplyr::bind_cols(dose_response_profiles[[mod]], model_selected = modelSelc[[1]])
      }
    }

    # Initialize prediction storage
    models <- attr(pObj, "model")
    S <- attr(pObj, "S")
    off <- attr(pObj, "off")
    scal <- attr(pObj, "scal")
    coefs <- attr(pObj$sim, "coefs")

    # Setting code up for computing uncertainty predictions.
    nDraw <- 2000
    pred_var_profiles <- list()
    dose_response_simulated <- list()
    dose_response_quantiles_long <- list()
    quantile_probs <- seq(0.05, 0.95, by = 0.05)

    for (true_mod in names(modelSel)) {
      pred_var_profiles[[true_mod]] <- matrix(NA, nrow = nSim, ncol = length(doses), dimnames = list(NULL, doses))
      dose_response_simulated[[true_mod]] <- array(NA, dim = c(nSim, nDraw, length(doses)), dimnames = list(NULL, NULL, doses))
      dose_response_quantiles_long[[true_mod]] <- vector("list", length = nSim)
    }

    for (true_mod in names(modelSel)) {
      # true_mod = "emax"
      for (i in seq_len(nSim)) {
        # i = 1
        mod <- modelSel[[true_mod]][i]

        cf <- try(coefs[[mod]][[i]], silent = TRUE)
        if (inherits(cf, "try-error") || is.null(cf)) next

        F <- DoseFinding:::gradCalc(mod, cf, doses, off, scal)
        if (!is.matrix(F) || nrow(F) != length(doses)) next
        if (!is.matrix(S) || ncol(S) != nrow(F)) next

        V <- try(solve(t(F) %*% solve(S) %*% F), silent = TRUE)
        if (inherits(V, "try-error")) next

        pred_mean <- dose_response_profiles[[true_mod]][i, ]
        pred_var  <- DoseFinding:::getPredVar(mod, cf, V, pDose = doses, off = off, scal = scal)
        pred_var_profiles[[true_mod]][i, ] <- pred_var

        sim_mat <- sapply(seq_along(pred_mean)[-length(pred_mean)], function(j) {
          rnorm(nDraw, mean = pred_mean[[j]], sd = sqrt(pred_var[j]))
        })
        dose_response_simulated[[true_mod]][i, , ] <- sim_mat

        qmat <- apply(sim_mat, 2, quantile, probs = quantile_probs)
        qdf <- as.data.frame(t(qmat))
        qdf$dose <- doses#[-length(doses)]
        qdf$simulation <- i
        qdf_long <- tidyr::pivot_longer(qdf, -c(dose, simulation), names_to = "quantile", values_to = "value")
        dose_response_quantiles_long[[true_mod]][[i]] <- qdf_long
      }
    }

    # Combine simulation-level quantile results into a single dataframe for each true model
    dose_response_quantiles_long <- purrr::map(dose_response_quantiles_long, dplyr::bind_rows)

    # -------------------------------------------------------------------------------------------------
    # Posterior power Module
    # -------------------------------------------------------------------------------------------------
    dose_response_profiles <- DoseFinding:::getSimEst(pObj, type = "dose-response", doseSeq = doses)

    posterior_power <- purrr::map_dfr(dose_response_profiles, function(mat) {
      dose_names <- as.numeric(colnames(mat))
      placebo_col <- which.min(dose_names)
      top_col <- which.max(dose_names)
      diff <- mat[, top_col] - mat[, placebo_col]
      prob_gt_0 <- mean(diff > 0)
      prob_gt_threshold <- mean(diff > go_threshold)
      prob_diff <- mean(abs(diff) > separation_threshold)
      tibble::tibble(
        prob_top_gt_placebo = prob_gt_0,
        prob_top_gt_placebo_threshold = prob_gt_threshold,
        prob_top_diff_placebo_threshold = prob_diff
      )
    }, .id = "model")

    posterior_power_grand <- purrr::map_dfr(dose_response_profiles, function(mat) {
      dose_names <- as.numeric(colnames(mat))
      comparisons <- expand.grid(
        comp_from = seq_along(dose_names),
        comp_to = seq_along(dose_names)
      ) |>
        dplyr::filter(comp_from != comp_to) |>
        dplyr::mutate(
          from_dose = dose_names[comp_from],
          to_dose = dose_names[comp_to],
          contrast = paste0(from_dose, " vs ", to_dose),
          prob_A_gt_B = purrr::map2_dbl(comp_from, comp_to, ~ mean(mat[, .x] > mat[, .y])),
          prob_A_gt_B_threshold = purrr::map2_dbl(comp_from, comp_to, ~ mean((mat[, .x] - mat[, .y]) > separation_threshold)),
          prob_diff_gt_threshold = purrr::map2_dbl(comp_from, comp_to, ~ mean(abs(mat[, .x] - mat[, .y]) > separation_threshold))
        )
      comparisons |>
        dplyr::select(contrast, from_dose, to_dose, prob_A_gt_B, prob_A_gt_B_threshold, prob_diff_gt_threshold)
    }, .id = "model")

    est_td_data <- DoseFinding:::getSimEst(pObj, "TD", Delta = go_threshold, direction = direction)
    true_td <- DoseFinding:::TD(altMods, Delta = go_threshold, TDtype = "continuous", direction = direction)
    td_summary <- purrr::map_dfr(est_td_data, function(x) {
      tibble::tibble(
        q025 = quantile(x, 0.025, na.rm = TRUE),
        q05  = quantile(x, 0.05, na.rm = TRUE),
        q10  = quantile(x, 0.10, na.rm = TRUE),
        q50  = quantile(x, 0.50, na.rm = TRUE),
        q90  = quantile(x, 0.90, na.rm = TRUE),
        q95  = quantile(x, 0.95, na.rm = TRUE),
        q975 = quantile(x, 0.975, na.rm = TRUE),
        mean = mean(x, na.rm = TRUE)
      )
    }, .id = "model") |>
      dplyr::mutate(true_target_dose = true_td[model])

    td_summary_formatted <- td_summary |>
      dplyr::transmute(
        Model = model,
        MaxEff = MaxEff_name,
        `True Target Dose` = true_target_dose,
        `Mean Estimated Dose (95% CI)` = sprintf("%.1f (%.1f–%.1f)", mean, q025, q975),
        `Mean Estimated Dose (90% CI)` = sprintf("%.1f (%.1f–%.1f)", mean, q05, q95),
        `Mean Estimated Dose (80% CI)` = sprintf("%.1f (%.1f–%.1f)", mean, q10, q90)
      )

    # -------------------------------------------------------------------------------------------------
    # Output collection module
    # -------------------------------------------------------------------------------------------------

    results[[MaxEff_name]] <- list(
      target_dose_estimation = list(
        data = summary_data,
        description = "Quantile summary of estimated target doses by model and Delta"
      ),
      posterior_probabilities = list(
        data = posterior_power,
        description = "Posterior probabilities that top dose exceeds placebo by various thresholds"
      ),
      pairwise_comparisons = list(
        data = posterior_power_grand,
        description = "Pairwise comparisons between all dose levels across models"
      ),
      target_dose_summary = list(
        table = td_summary_formatted,
        description = "Formatted table of true and estimated target dose with multiple CIs"
      ),
      planMod_output = pObj,
      predictions_df = list(
        data = dose_response_profiles,
        description = "Predictions of simulated dose-response profiles across models."
      ),
      predictions_df_var = list(
        data = dose_response_quantiles_long,
        description = "Predictive quantile summaries of simulated dose-response profiles across models. Each entry contains per-simulation quantiles (e.g., 5%, 10%, ..., 95%) for every dose and fitted model. This allows uncertainty characterization around predicted mean responses."
      )
    )
  }
  class(results) <- "td_power_simulation"
  return(results)
}
