\encoding{UTF-8}
\name{Plant values}
\alias{plant_basalArea}
\alias{plant_characterParameter}
\alias{plant_cover}
\alias{plant_crownBaseHeight}
\alias{plant_crownLength}
\alias{plant_crownRatio}
\alias{plant_density}
\alias{plant_equilibriumLeafLitter}
\alias{plant_equilibriumSmallBranchLitter}
\alias{plant_foliarBiomass}
\alias{plant_fuel}
\alias{plant_height}
\alias{plant_ID}
\alias{plant_parameter}
\alias{plant_LAI}
\alias{plant_largerTreeBasalArea}
\alias{plant_phytovolume}
\alias{plant_shrubIndividualArea}
\alias{plant_species}
\alias{plant_speciesName}

\title{Plant description functions}
\description{
Functions to calculate attributes of plants in a \code{\link{forest}} object.
}
\usage{
plant_basalArea(x)
plant_largerTreeBasalArea(x)
plant_characterParameter(x, SpParams, parName)
plant_cover(x)
plant_crownBaseHeight(x, SpParams, mode = "MED")
plant_crownLength(x, SpParams, mode = "MED")
plant_crownRatio(x, SpParams, mode = "MED")
plant_density(x, SpParams, mode = "MED")
plant_equilibriumLeafLitter(x, SpParams, AET = 800, mode = "MED")
plant_equilibriumSmallBranchLitter(x, SpParams, 
                                   smallBranchDecompositionRate = 0.81, mode = "MED")
plant_foliarBiomass(x, SpParams, gdd = NA, mode = "MED")
plant_fuel(x, SpParams, gdd = NA, includeDead = TRUE, mode = "MED")
plant_height(x)
plant_ID(x, treeOffset = 0, shrubOffset = 0)
plant_LAI(x, SpParams, gdd = NA, mode = "MED")
plant_parameter(x, SpParams, parName, fillMissing = TRUE)
plant_phytovolume(x, SpParams)
plant_species(x)
plant_speciesName(x, SpParams)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{parName}{A string with a parameter name.}
  \item{mode}{Calculation mode, either "MED" or "US".}
  \item{gdd}{Growth degree days (to account for leaf phenology effects).}
  \item{AET}{Actual annual evapotranspiration (in mm).}
  \item{smallBranchDecompositionRate}{Decomposition rate of small branches.}
  \item{includeDead}{A flag to indicate that standing dead fuels (dead branches) are included.}
  \item{treeOffset, shrubOffset}{Integers to offset cohort IDs.}
  \item{fillMissing}{A boolean flag to try imputation on missing values.}
}
\value{
A vector with values for each plant of the input \code{\link{forest}} object:
\itemize{
\item{\code{plant_basalArea}: Tree basal area (m2/ha).}
\item{\code{plant_largerTreeBasalArea}: Basal area (m2/ha) of trees larger (in diameter) than the tree. Half of the trees of the same record are included.}
\item{\code{plant_characterParameter}: The parameter values of each plant, as strings.}
\item{\code{plant_cover}: Shrub cover (in percent).}
\item{\code{plant_crownBaseHeight}: The height corresponding to the start of the crown (in cm).}
\item{\code{plant_crownLength}: The difference between crown base height and total height (in cm).}
\item{\code{plant_crownRatio}: The ratio between crown length and total height (between 0 and 1).}
\item{\code{plant_density}: Plant density (ind/ha). Tree density is directly taken from the forest object, while the shrub density is estimated from cover and height by calculating the area of a single individual.}
\item{\code{plant_equilibriumLeafLitter}: Litter biomass of leaves at equilibrium (in kg/m2).}
\item{\code{plant_equilibriumSmallBranchLitter}: Litter biomass of small branches (< 6.35 mm diameter) at equilibrium (in kg/m2).}
\item{\code{plant_foliarBiomass}: Standing biomass of leaves (in kg/m2).}
\item{\code{plant_fuel}: Fine fuel load (in kg/m2).}
\item{\code{plant_height}: Total height (in cm).}
\item{\code{plant_ID}: Cohort coding for simulation functions (concatenation of 'T' (Trees) or 'S' (Shrub), cohort index and species index).}
\item{\code{plant_LAI}: Leaf area index (m2/m2).}
\item{\code{plant_parameter}: The parameter values of each plant, as numeric.}
\item{\code{plant_phytovolume}: Shrub phytovolume (m3/m2).}
\item{\code{plant_species}: Species identity integer (indices start with 0).}
\item{\code{plant_speciesName}: String with species taxonomic name (or a functional group).}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{spwb}}, \code{\link{forest}}, \code{\link{summary.forest}}
}
\examples{
#Default species parameterization
data(SpParamsMED)

#Load example plot
data(exampleforestMED)

#A short way to obtain total basal area
sum(plant_basalArea(exampleforestMED), na.rm=TRUE)

#The same forest level function for LAI
sum(plant_LAI(exampleforestMED, SpParamsMED))

#The same forest level function for fuel loading
sum(plant_fuel(exampleforestMED, SpParamsMED))

#Summary function for 'forest' objects can be also used
summary(exampleforestMED, SpParamsMED)

plant_speciesName(exampleforestMED, SpParamsMED)

plant_ID(exampleforestMED)
}
