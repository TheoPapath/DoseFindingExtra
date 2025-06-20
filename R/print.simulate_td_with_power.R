#' Print method for `simulate_td_with_power` output
#'
#' Provides a clear summary of the results from `simulate_td_with_power()`,
#' including the evaluated scenarios and available result components.
#'
#' @param x A `simulate_td_with_power` object as returned by `simulate_td_with_power()`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#' @export
print.simulate_td_with_power <- function(x, ...) {
  cat("\nTarget Dose Simulation Results Summary\n")
  cat(strrep("-", 42), "\n")

  if (!is.list(x)) {
    cat("[Error] Input object is not a valid simulation result.\n")
    return(invisible(x))
  }

  scenario_names <- names(x)
  cat("\nScenarios Evaluated (e.g., MaxEff values):\n")
  for (sc in scenario_names) cat("  • ", sc, "\n")

  for (sc in scenario_names) {
    cat("\n--- Results for Scenario:", sc, "---\n")
    result <- x[[sc]]
    for (comp_name in names(result)) {
      entry <- result[[comp_name]]
      if (is.list(entry) && !is.null(entry$description)) {
        cat(sprintf("  • %-25s : %s\n", comp_name, entry$description))
      } else {
        cat(sprintf("  • %-25s : [no description]\n", comp_name))
      }
    }
  }
  cat("\nTip: Use `results$<scenario>$<component>$data` to extract detailed outputs.\n\n")
  invisible(x)
}
