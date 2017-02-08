\encoding{UTF-8}
\name{Plant values}
\alias{plant.BasalArea}
\alias{plant.CharacterParameter}
\alias{plant.Cover}
\alias{plant.CrownBaseHeight}
\alias{plant.CrownLength}
\alias{plant.CrownRatio}
\alias{plant.Density}
\alias{plant.EquilibriumLeafLitter}
\alias{plant.EquilibriumSmallBranchLitter}
\alias{plant.FoliarBiomass}
\alias{plant.Fuel}
\alias{plant.Height}
\alias{plant.Parameter}
\alias{plant.LAI}
\alias{plant.LargerTreeBasalArea}
\alias{plant.Phytovolume}
\alias{plant.Species}
\alias{plant.SpeciesName}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plant description functions}
\description{
Functions to calculate attributes of plants in a \code{\link{forest}} object.
}
\usage{
plant.BasalArea(x)
plant.LargerTreeBasalArea(x)
plant.CharacterParameter(x, SpParams, parName)
plant.Cover(x)
plant.CrownBaseHeight(x, SpParams)
plant.CrownLength(x, SpParams)
plant.CrownRatio(x, SpParams)
plant.Density(x, SpParams)
plant.EquilibriumLeafLitter(x, SpParams, AET = 800)
plant.EquilibriumSmallBranchLitter(x, SpParams, 
                                   smallBranchDecompositionRate = 0.81)
plant.FoliarBiomass(x, SpParams, gdd = NA)
plant.Fuel(x, SpParams, gdd = NA, includeDead = TRUE)
plant.Height(x)
plant.LAI(x, SpParams, gdd = NA)
plant.Parameter(x, SpParams, parName)
plant.Phytovolume(x, SpParams)
plant.Species(x)
plant.SpeciesName(x, SpParams)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{parName}{A string with a parameter name.}
  \item{gdd}{Growth degree days (to account for leaf phenology effects).}
  \item{AET}{Actual annual evapotranspiration (in mm).}
  \item{smallBranchDecompositionRate}{Decomposition rate of small branches.}
  \item{includeDead}{A flag to indicate that standing dead fuels (dead branches) are included.}
}
\value{
A vector with values for each plant of the input \code{\link{forest}} object:
\itemize{
\item{\code{plant.BasalArea}: Tree basal area (m2/ha).}
\item{\code{plant.LargerTreeBasalArea}: Basal area (m2/ha) of trees larger (in diameter) than the tree. Half of the trees of the same record are included.}
\item{\code{plant.CharacterParameter}: The parameter values of each plant, as strings.}
\item{\code{plant.Cover}: Shrub cover (in percent).}
\item{\code{plant.CrownBaseHeight}: The height corresponding to the start of the crown (in cm).}
\item{\code{plant.CrownLength}: The difference between crown base height and total height (in cm).}
\item{\code{plant.CrownRatio}: The ratio between crown length and total height (between 0 and 1).}
\item{\code{plant.Density}: Plant density (ind/ha). Tree density is directly taken from the forest object, while the shrub density is estimated from cover and height by calculating the area of a single individual.}
\item{\code{plant.EquilibriumLeafLitter}: Litter biomass of leaves at equilibrium (in kg/m2).}
\item{\code{plant.EquilibriumSmallBranchLitter}: Litter biomass of small branches (< 6.35 mm diameter) at equilibrium (in kg/m2).}
\item{\code{plant.FoliarBiomass}: Standing biomass of leaves (in kg/m2).}
\item{\code{plant.Fuel}: Fine fuel load (in kg/m2).}
\item{\code{plant.Height}: Total height (in cm).}
\item{\code{plant.LAI}: Leaf area index (m2/m2).}
\item{\code{plant.Parameter}: The parameter values of each plant, as numeric.}
\item{\code{plant.Phytovolume}: Shrub phytovolume (m3/m2).}
\item{\code{plant.Species}: Species identity integer (indices start with 0).}
\item{\code{plant.SpeciesName}: String with species taxonomic name (or a functional group).}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swb}}, \code{\link{forest}}, \code{\link{summary.forest}}
}
\examples{
#Default species parameterization
data(SpParamsMED)

#Load example plot
data(exampleforest)

#A short way to obtain total basal area
sum(plant.BasalArea(exampleforest), na.rm=TRUE)

#The same forest level function for LAI
sum(plant.LAI(exampleforest, SpParamsMED))

#The same forest level function for fuel loading
sum(plant.Fuel(exampleforest, SpParamsMED))

#Summary function for 'forest' objects can be also used
summary(exampleforest, SpParamsMED)

plant.SpeciesName(exampleforest, SpParamsMED)
}