\name{spwb_stress}
\alias{spwb_stress}
\title{
Drought stress indices
}
\description{
Allows calculating annual-based or monthly-based drought stress indices from \code{\link{spwb}} objects.
}
\usage{
spwb_stress(x, index = "NDD", freq = "years", bySpecies = FALSE, draw = TRUE)
}
\arguments{
  \item{x}{An object of class \code{spwb}.}
  \item{index}{A string with the index to be calculated, either \code{"DI"}, \code{"NDD"}, \code{"ADS"}, \code{"MDS"} or \code{"WSI"} (see details).}
  \item{freq}{Frequency of stress statistics (see \code{\link{cut.Date}}). Normally, either \code{"years"} or \code{"months"} for yearly-based or monthly-based indices.}
  \item{bySpecies}{Allows aggregating output by species.}
  \item{draw}{A boolean flag to indicate that a plot should be returned.}
}
\details{
The currently available drought stress indices are:
  \itemize{
    \item{\code{"ADS"}:}{ Average of daily drought stress values for the period considered.}
    \item{\code{"MDS"}:}{ Maximum daily drought stress during the period considered.}
    \item{\code{"DI"}:}{ Drought intensity, as defined in De \enc{Cáceres}{Caceres} et al. (2015).}
    \item{\code{"NDD"}:}{ Number of drought days, as defined in De \enc{Cáceres}{Caceres} et al. (2015).}
    \item{\code{"WSI"}:}{ Water stress integral, as defined in Myers (1988).}
  }
}
\value{
A data frame with periods (e.g., years or months) in rows and plant cohorts (or species) in columns. Values are the calculated stress index. If \code{draw=TRUE} a ggplot is returned instead.
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).

Myers BJ (1988) Water stress integral - a link between short-term stress and long-term growth. Tree Physiology 4: 315–323 (doi: 10.1093/treephys/4.4.315)
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}

\seealso{
\code{\link{spwb}}, \code{\link{summary.spwb}}
}

