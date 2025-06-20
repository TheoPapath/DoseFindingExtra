% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplotDoseSims.R
\name{ggplotDoseSims}
\alias{ggplotDoseSims}
\title{Plot Simulated Dose Estimates with Confidence Intervals}
\usage{
ggplotDoseSims(
  x,
  type = c("ED", "TD"),
  p,
  Delta,
  xlab = "Dose",
  ci_levels = c(0.95, 0.8, 0.5)
)
}
\arguments{
\item{x}{A simulation result object, typically from \code{planDoseFinding()}, with attributes \code{altModels} and \code{doses}.}

\item{type}{Character, either "ED" or "TD" to specify whether to plot effective or toxic dose estimates.}

\item{p}{Numeric, the target quantile (used if \code{type = "ED"}).}

\item{Delta}{Numeric, the minimal relevant effect size (used if \code{type = "TD"}).}

\item{xlab}{Label for the x-axis.}

\item{ci_levels}{Numeric vector of confidence levels to display. Must be subset of c(0.95, 0.80, 0.50).}
}
\value{
A ggplot object showing interval estimates and reference markers.
}
\description{
This function creates a horizontal interval plot showing simulated Effective Doses (ED) or Toxic Doses (TD),
with confidence intervals at 95\%, 80\%, and 50\% levels, and overlayed true values.
}
\examples{
# Assuming `x` is a simulation object from planDoseFinding()
# ggplotDoseSims(x, type = "ED", p = 0.8)

}
