% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resistances.R
\name{resistances}
\alias{resistances}
\title{Soil-plant resistances}
\usage{
resistances(
  x,
  cohort,
  relative = FALSE,
  draw = FALSE,
  cumulative = FALSE,
  xlab = NULL,
  ylab = NULL
)
}
\arguments{
\item{x}{An object of class \code{\link{spwb}}, \code{\link{pwb}}, \code{\link{growth}} or \code{\link{fordyn}}.
The function only works with the result of simulations with \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}.}

\item{cohort}{An string indicating the cohort for which resistances are desired.}

\item{relative}{A boolean flag to indicate that relative percentages are desired as output.}

\item{draw}{A boolean flag to indicate that a plot should be drawn (only pathway resistances, without discriminating between soil layers).}

\item{cumulative}{A flag to indicate that drawn series should be cumulative.}

\item{xlab}{x-axis label.}

\item{ylab}{y-axis label.}
}
\value{
If \code{draw = FALSE}, the function returns list with three items:
\itemize{
\item{\code{pathway}: A data frame with dates in rows and resistance segments in columns (Rhizosphere, Root, Stem and Leaf). }
\item{\code{root}: A data frame with dates in rows and root resistances for soil layers in columns.}
\item{\code{rhizosphere}: A data frame with dates in rows and rhizosphere resistances for soil layers in columns. }
}
Values depend on whether \code{relative = TRUE} (percentages) or \code{relative = FALSE} (absolute resistance values).

If \code{draw = TRUE}, a plot object is returned showing the time series of pathway resistances.
}
\description{
Calculates and draws rhizosphere, root, stem and leaf resistances for simulation time steps
}
\details{
The function makes internal calls to \code{\link{hydraulics_soilPlantResistancesWeibull}} or
\code{\link{hydraulics_soilPlantResistancesSigmoid}} depending on the value of \code{transpirationMode} in \code{x}.
}
\seealso{
\code{\link{waterUseEfficiency}}, \code{\link{droughtStress}}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF

\enc{Léa}{Lea} Veuillen, INRAE-URFM
}
