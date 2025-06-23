#' Convert Between Randomization Ratios and Percentages
#'
#' Converts allocation ratios (e.g., 5:1:1:5) to percentage allocations and vice versa.
#'
#' @param x A numeric vector. If `type = "ratio"`, this should be a vector of allocation ratios (e.g., `c(5, 1, 1, 5)`). If `type = "percent"`, this should be a vector of percentages (e.g., `c(41.67, 8.33, 8.33, 41.67)`).
#' @param type A character string specifying the type of input. Either `"ratio"` (default) or `"percent"`.
#' @param digits Number of digits to round percentages to when `type = "ratio"`. Default is 2.
#'
#' @return A numeric vector:
#' - If `type = "ratio"`, returns percentage allocations summing to ~100%.
#' - If `type = "percent"`, returns the simplified whole-number allocation ratio.
#'
#' @examples
#' # Convert a ratio to percentages
#' allocation_converter(c(5, 1, 1, 5), type = "ratio")
#'
#' # Convert percentages to simplified ratio
#' allocation_converter(c(41.67, 8.33, 8.33, 41.67), type = "percent")
#'
#' # More irregular percentages
#' allocation_converter(c(30, 10, 20, 40), type = "percent")
#'
#' @export
allocation_converter <- function(x, type = c("ratio", "percent"), digits = 2) {
  type <- match.arg(type)

  if (type == "ratio") {
    total <- sum(x)
    percent <- round(100 * x / total, digits)
    return(percent)
  }

  if (type == "percent") {
    # Normalize to smallest component = 1
    scaled <- round(x / min(x))

    # Greatest common divisor function
    gcd_vec <- function(v) Reduce(function(a, b) if (b == 0) a else Recall(b, a %% b), v)

    # Simplify
    divisor <- gcd_vec(scaled)
    ratio <- scaled / divisor
    return(as.integer(ratio))
  }
}
