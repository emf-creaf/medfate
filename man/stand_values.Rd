\encoding{UTF-8}
\name{Stand values}
\alias{stand_basalArea}
\alias{stand_foliarBiomass}
\alias{stand_fuel}
\alias{stand_LAI}
\alias{stand_phytovolume}


\title{Stand values}
\description{
Functions to calculate stand attributes of a \code{\link{forest}} object.
}
\usage{
stand_basalArea(x)
stand_foliarBiomass(x, SpParams, gdd = NA, mode = "MED")
stand_fuel(x, SpParams, gdd = NA, includeDead = TRUE, mode = "MED")
stand_LAI(x, SpParams, gdd = NA, mode = "MED")
stand_phytovolume(x, SpParams)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{gdd}{Growth degree days (to account for leaf phenology effects).}
  \item{includeDead}{A flag to indicate that standing dead fuels (dead branches) are included.}
  \item{mode}{Calculation mode, either "MED" or "US".}
}
\value{
\itemize{
\item{\code{stand_basalArea}: Stand basal area (m2/ha).}
\item{\code{stand_foliarBiomass}: Standing biomass of leaves (in kg/m2).}
\item{\code{stand_fuel}: Stand fine fuel load (in kg/m2).}
\item{\code{stand_LAI}: Stand leaf area index (m2/m2).}
\item{\code{stand_phytovolume}: Stand shrub phytovolume (m3/m2).}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{forest}}, \code{\link{plant_basalArea}}, \code{\link{summary.forest}}
}
\examples{
#Default species parameterization
data(SpParamsMED)

#Load example plot
data(exampleforestMED)

#A short way to obtain total basal area
stand_basalArea(exampleforestMED)
}