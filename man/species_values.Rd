\encoding{UTF-8}
\name{Species values}
\alias{species_basalArea}
\alias{species_cover}
\alias{species_density}
\alias{species_foliarBiomass}
\alias{species_fuel}
\alias{species_LAI}
\alias{species_phytovolume}
\alias{species_parameter}
\alias{species_characterParameter}


\title{Species description functions}
\description{
Functions to calculate attributes of a \code{\link{forest}} object by species or to extract species parameters from a species parameter table (\code{\link{SpParamsMED}}).
}
\usage{
species_basalArea(x, SpParams)
species_cover(x, SpParams)
species_density(x, SpParams, mode = "MED")
species_foliarBiomass(x, SpParams, gdd = NA, mode = "MED")
species_fuel(x, SpParams, gdd = NA, includeDead = TRUE, mode = "MED")
species_LAI(x, SpParams, gdd = NA, mode = "MED")
species_phytovolume(x, SpParams)
species_parameter(SP, SpParams, parName)
species_characterParameter(SP, SpParams, parName)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{gdd}{Growth degree days (to account for leaf phenology effects).}
  \item{includeDead}{A flag to indicate that standing dead fuels (dead branches) are included.}
  \item{mode}{Calculation mode, either "MED" or "US".}
  \item{SP}{An integer vector of species codes.}
  \item{parName}{A string with a parameter name.}
}
\value{
A vector with values for each species in \code{SpParams}:
\itemize{
\item{\code{species_basalArea}: Species basal area (m2/ha).}
\item{\code{species_cover}: Shrub cover (in percent).}
\item{\code{species_density}: Plant density (ind/ha). Tree density is directly taken from the forest object, while the shrub density is estimated from cover and height by calculating the area of a single individual.}
\item{\code{species_foliarBiomass}: Standing biomass of leaves (in kg/m2).}
\item{\code{species_fuel}: Fine fuel load (in kg/m2).}
\item{\code{species_LAI}: Leaf area index (m2/m2).}
\item{\code{species_phytovolume}: Shrub phytovolume (m3/m2).}
\item{\code{species_parameter}: A numeric vector with the parameter values of each input species.}
\item{\code{species_characterParameter}: A character vector with the parameter values of each input species.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{spwb}}, \code{\link{forest}}, \code{\link{plant_basalArea}}, \code{\link{summary.forest}}
}
\examples{
# Default species parameterization
data(SpParamsMED)

# Load example plot
data(exampleforestMED)

# Species basal area in the forest plot
species_basalArea(exampleforestMED, SpParamsMED)

# Value of parameter `Psi_Extract` for species 54 (Pinus halepensis)
# and 68 (Quercus ilex)
species_parameter(c(54, 68), SpParamsMED, "Psi_Extract")
}
