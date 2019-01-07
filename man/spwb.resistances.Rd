\name{spwb.resistances}
\alias{spwb.resistances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Soil-plant resistances
}
\description{
Calculates and draws rhizosphere, root, stem and leaf resistances for simulation time steps
}
\usage{
spwb.resistances(x, cohort = 1, relative = FALSE, draw = FALSE, cumulative = FALSE, yearAxis = FALSE, xlab = NULL, ylab = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{spwb}. The function only works with the result of simulations with \code{transpirationMode = "Complex".}
}
  \item{cohort}{
An integer index indicating the cohort for which resistances are desired (by default the first cohort).
}
  \item{relative}{
A boolean flag to indicate that relative percentages are desired as output
}
  \item{draw}{
A boolean flag to indicate that a plot should be drawn.
}
  \item{cumulative}{
A flag to indicate that drawn series should be cumulative.
}
  \item{yearAxis}{
A boolean to indicate whether the units of the x-axis are years (by default they are dates).}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
}
\details{
The function makes internal calls to \code{\link{hydraulics.soilPlantResistances}}.
}
\value{
A data frame with dates in rows and resistance segments in columns (Rhizosphere, Root, Stem and Leaf). Values depend on whether \code{relative = TRUE} (percentages) or \code{relative = FALSE} (absolute resistance values). If \code{draw = TRUE} the calculated resistances are returned as an invisible object.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{spwb}}, \code{\link{plot.spwb}}
}