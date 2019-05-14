\name{spwb_resistances}
\alias{spwb_resistances}
\title{
Soil-plant resistances
}
\description{
Calculates and draws rhizosphere, root, stem and leaf resistances for simulation time steps
}
\usage{
spwb_resistances(x, cohort = 1, relative = FALSE, draw = FALSE, 
                 cumulative = FALSE, xlab = NULL, ylab = NULL)
}
\arguments{
  \item{x}{An object of class \code{spwb}. The function only works with the result of simulations with \code{transpirationMode = "Sperry"}.}
  \item{cohort}{An integer index indicating the cohort for which resistances are desired (by default the first cohort).}
  \item{relative}{A boolean flag to indicate that relative percentages are desired as output}
  \item{draw}{A boolean flag to indicate that a plot should be drawn.}
  \item{cumulative}{A flag to indicate that drawn series should be cumulative.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
}
\details{
The function makes internal calls to \code{\link{hydraulics_soilPlantResistances}}.
}
\value{
A data frame with dates in rows and resistance segments in columns (Rhizosphere, Root, Stem and Leaf). Values depend on whether \code{relative = TRUE} (percentages) or \code{relative = FALSE} (absolute resistance values). If \code{draw = TRUE} then a plot object is returned.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{spwb}}, \code{\link{plot.spwb}}
}