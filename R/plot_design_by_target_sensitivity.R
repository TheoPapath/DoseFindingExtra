#' Plot Optimal Dose Allocations by Target Effect Level
#'
#' Visualizes how optimal dose allocation proportions shift across a range of target effect levels (Delta),
#' using three design criteria: D-optimal, TD-optimal, and the hybrid Dopt&TD.
#'
#' @param fmods A model object of class `Mods` from the DoseFinding package.
#' @param probs A numeric vector of model probabilities. If NULL, defaults to 1 for a single model or equal weights for multiple models.
#' @param delta_seq A numeric sequence specifying the range of target effect levels (Delta) to evaluate. Default is from 0.5 to 1.5.
#'
#' @return An invisible list with three components:
#' \describe{
#'   \item{data}{A data frame of optimal allocation proportions for each dose, Delta, and design criterion.}
#'   \item{plot}{A `ggplot` object showing the stacked area plot, faceted by design criterion.}
#'   \item{animation}{A `gganimate` object showing allocation shifts across Delta (if supported).}
#' }
#'
#' @examples
#' library(DoseFinding)
#' doses <- c(0, 5, 40, 80, 160)
#' mods <- Mods(emax = 50, doses = doses, placEff = 0, maxEff = 1)
#' plot_design_by_target_sensitivity(mods)
#'
#' # With multiple models
#' mods <- Mods(
#'   emax = c(20, 50, 80, 200),
#'   sigEmax = c(30, 3.5),
#'   quadratic = -0.004,
#'   placEff = 0,
#'   maxEff = 1,
#'   doses = doses
#' )
#' probs <- c(0.3, rep(0.1, 5) + 0.2/5)
#' plot_design_by_target_sensitivity(mods, probs)
#'
#' @export
plot_design_by_target_sensitivity <- function(fmods, probs = NULL, delta_seq = seq(0.5, 1.5, by = 0.1)) {
  # Ensure required packages
  if (!requireNamespace("DoseFinding", quietly = TRUE)) stop("Package 'DoseFinding' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  # if (!requireNamespace("gganimate", quietly = TRUE)) stop("Package 'gganimate' is required for animation.")

  # Infer default probs if not provided
  if (is.null(probs)) {
    if (inherits(fmods, "Mods") && length(attr(fmods, "model")) > 1) {
      warning("Multiple models detected but 'probs' not provided. Assuming equal probability across models.")
      probs <- rep(1 / length(attr(fmods, "model")), length(attr(fmods, "model")))
    } else {
      probs <- 1
    }
  }

  # Initialize empty list to store results
  all_designs <- list()
  counter <- 1

  for (d in delta_seq) {
    for (crit in c("Dopt", "TD", "Dopt&TD")) {
      message(sprintf("Running optDesign with Delta = %.2f and criterion = %s", d, crit))
      result <- tryCatch({
        des <- DoseFinding::optDesign(fmods, probs = probs, designCrit = crit, Delta = d)

        if (length(des$design) == 0) {
          message(sprintf("Empty design returned for Delta = %.2f and criterion = %s", d, crit))
          NULL
        } else {
          data.frame(
            Delta = rep(d, length(des$design)),
            Dose = as.numeric(des$doses),
            Proportion = rndDesign(des, n = 100),
            Criterion = crit
          )
        }
      }, error = function(e) {
        message(sprintf("Error at Delta = %.2f with criterion = %s: %s", d, crit, e$message))
        NULL
      })

      if (!is.null(result)) {
        all_designs[[counter]] <- result
        counter <- counter + 1
      }
    }
  }

  # Combine all valid results
  if (length(all_designs) == 0) {
    stop("No valid designs returned across the Delta range.")
  }

  design_df <- do.call(rbind, all_designs)

  # Ensure Dose = 0 appears at the bottom
  design_df$Dose <- factor(design_df$Dose, levels = sort(unique(design_df$Dose)))

  nejm_colors <- c(
    "#2166AC", "#B2182B", "#4DAF4A", "#984EA3", "#FF7F00",  # NEJM base colors
    "#999999", "#E41A1C", "#377EB8", "#A65628",             # additional NEJM-like
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"   # extended scientific palette
  )

  # Static stacked area plot with facets for each criterion
  static_plot <- ggplot2::ggplot(design_df, aes(x = Delta, y = Proportion, fill = Dose)) +
    geom_area(position = "stack", alpha = 0.8) +
    facet_wrap(~ Criterion) +
    scale_fill_manual(values = nejm_colors) +
    labs(
      title = "Shifts in Optimal Dose Allocation vs. Target Effect (Delta)",
      x = expression(Target ~ Effect ~ (Delta)),
      y = "Allocation Proportion",
      fill = "Dose"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )

  # # Animated version
  # animated_plot <- static_plot +
  #   gganimate::transition_states(Delta, transition_length = 2, state_length = 1) +
  #   labs(title = "Delta: {closest_state}")

  # Return both plot and data invisibly
  print(static_plot)
  invisible(list(data = design_df, plot = static_plot
                 # , animation = animated_plot
  ))
}

# Example usage:
# library(DoseFinding)
# doses <- c(0, 5, 40, 80, 160)
# mods <- Mods(emax = 50, doses = doses, placEff = 0, maxEff = 1)
# plot_design_by_target_sensitivity(mods)

# With multiple models:
# mods <- Mods(
#   emax = c(20, 50, 80, 200),
#   sigEmax = c(30, 3.5),
#   quadratic = -0.004,
#   placEff = 0,
#   maxEff = 1,
#   doses = doses
# )
# probs <- c(0.3, rep(0.1, 5) + 0.2/5)
# plot_design_by_target_sensitivity(mods, probs)
