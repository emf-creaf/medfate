\encoding{UTF-8}
\name{growth}
\alias{growth}
\title{Forest growth}
\description{
Function \code{growth} is a process-based model that performs energy, water and carbon balances; and determines changes in water/carbon pools, functional variables (leaf area, sapwood area, root area) and structural ones (tree diameter, tree height, shrub cover) for woody plant cohorts in a given forest stand during a period specified in the input climatic data. 
}
\usage{
growth(x, meteo,  latitude, elevation = NA, slope = NA, aspect = NA) 
}
\arguments{
  \item{x}{An object of class \code{\link{growthInput}}.}
  \item{meteo}{A data frame with daily meteorological data series. Row names of the data frame should correspond to date strings with format "yyyy-mm-dd" (see \code{\link{Date}}). When  \code{x$TranspirationMode = "Granier"} the following columns are required:
   \itemize{
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{Radiation}: Solar radiation (in MJ/m2/day), required only if \code{snowpack = TRUE}.}
     \item{\code{PET}: Potential evapotranspiration (in mm).}
     \item{\code{WindSpeed}: Wind speed (in m/s). If not available, this column can be left with \code{NA} values.}
    }
  When \code{x$TranspirationMode = "Sperry"} the following columns are required:
   \itemize{
     \item{\code{MeanTemperature}: Mean temperature(in degrees Celsius).}
     \item{\code{MinTemperature}: Minimum temperature (in degrees Celsius).}
     \item{\code{MaxTemperature}: Maximum temperature (in degrees Celsius).}
     \item{\code{MinRelativeHumidity}: Minimum relative humidity (in percent).}
     \item{\code{MaxRelativeHumidity}: Maximum relative humidity (in percent).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{Radiation}: Solar radiation (in MJ/m2/day).}
     \item{\code{WindSpeed}: Wind speed (in m/s). If not available, this column can be left with \code{NA} values.}
     \item{\code{CO2}: Atmospheric (abovecanopy) CO2 concentration (in ppm). This column may not exist, or can be left with \code{NA} values. In both cases simulations will assume a constant value specified in \code{\link{defaultControl}}.}
    }
  }
  \item{latitude}{Latitude (in degrees). Required when \code{x$TranspirationMode = "Sperry"}.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when \code{x$TranspirationMode = "Sperry"}. Elevation is also required for 'Granier' if snowpack dynamics are simulated.}
}
\details{
  Detailed model description is available in the vignettes section. Simulations using the 'Sperry' transpiration mode are computationally much more expensive than those using the simple transpiration mode. 
}
\value{
  A list of class 'growth' with the following elements:
  \itemize{
  \item{\code{"latitude"}: Latitude (in degrees) given as input.} 
  \item{\code{"topography"}: Vector with elevation, slope and aspect given as input.} 
  \item{\code{"growthInput"}: A copy of the object \code{x} of class \code{\link{growthInput}} given as input.}
  \item{\code{"WaterBalance"}: A data frame where different water balance variables (see \code{\link{spwb}}).}
  \item{\code{"EnergyBalance"}: A data frame with the daily values of energy balance components for the soil and the canopy (only for \code{transpirationMode = "Sperry"}; see \code{\link{spwb}}).}
  \item{\code{"Temperature"}: A data frame with the daily values of minimum/mean/maximum temperatures for the atmosphere (input), canopy and soil (only for \code{transpirationMode = "Sperry"}; see \code{\link{spwb}}).}
  \item{\code{"Soil"}: A data frame where different soil variables  (see \code{\link{spwb}}).}
  \item{\code{"Stand"}: A data frame where different stand-level variables (see \code{\link{spwb}}).}
  \item{\code{"Plants"}: A list of daily results for plant cohorts (see \code{\link{spwb}}).}
  \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: A list with daily results for sunlit and shade leaves (only for \code{transpirationMode = "Sperry"}; see \code{\link{spwb}}).}
  \item{\code{"PlantCarbonBalance"}: A list of daily carbon balance results for plant cohorts, with elements:}
    \itemize{
        \item{\code{"GrossPhotosynthesis"}: Daily gross photosynthesis per dry weight of living biomass (g gluc · g dry-1).}
        \item{\code{"MaintentanceRespiration"}: Daily maintenance respiration per dry weight of living biomass (g gluc · g dry-1).}
        \item{\code{"GrowthCosts"}: Daily growth costs per dry weight of living biomass (g gluc · g dry-1).}
        \item{\code{"RootExudation"}: Root exudation per dry weight of living biomass (g gluc · g dry-1).}    
        \item{\code{"CarbonBalance"}: Daily plant carbon balance (photosynthesis - maintenance respiration - growth costs - root exudation) per dry weight of living biomass (g gluc · g dry-1).}
        \item{\code{"SugarLeaf"}: Sugar concentration (mol·l-1) in leaves.}
        \item{\code{"StarchLeaf"}: Starch concentration (mol·l-1) in leaves.}
        \item{\code{"SugarSapwood"}: Sugar concentration (mol·l-1) in sapwood.}
        \item{\code{"StarchSapwood"}: Starch concentration (mol·l-1) in sapwood.}
        \item{\code{"SugarTransport"}:  Average instantaneous rate of carbon transferred between leaves and stem compartments via floem (mol gluc·s-1).}
        \item{\code{"LeafPI0"}: Osmotic potential at full turgor of symplastic leaf tissue (MPa).}
        \item{\code{"StemPI0"}: Osmotic potential at full turgor of symplastic stem tissue (MPa).}
      }
  \item{\code{"PlantStructure"}: A list of daily area and biomass values for compartments of plant cohorts, with elements:}
    \itemize{
        \item{\code{"LeafArea"}: Daily amount of leaf area (in m2) for an average individual of each plant cohort.}
        \item{\code{"SapwoodArea"}: Daily amount of sapwood area (in cm2) for an average individual of each plant cohort.}
        \item{\code{"FineRootArea"}: Daily amount of fine root area (in m2) for an average individual of each plant cohort (only for \code{transpirationMode = "Sperry"}).}
        \item{\code{"SapwoodBiomass"}: Daily amount of sapwood biomass (in g dry) for an average individual of each plant cohort.}
        \item{\code{"LeafBiomass"}: Daily amount of leaf biomass (in g dry) for an average individual of each plant cohort.}
        \item{\code{"FineRootBiomass"}: Daily amount of fine root biomass (in g dry) for an average individual of each plant cohort.}
        \item{\code{"LabileBiomass"}: Daily amount of labile C biomass, i.e. starch and sugars (in g dry) for an average individual of each plant cohort.}
        \item{\code{"TotalLivingBiomass"}: Daily amount of total living biomass, i.e. excluding heartwood, (in g dry) for an average individual of each plant cohort.}
    }
  \item{\code{"PlantGrowth"}: A list of daily growth results for plant cohorts, with elements:}
    \itemize{
        \item{\code{"LAgrowth"}: Leaf area growth relative to sapwood area (in m2·cm-2·day-1) for an average individual of each plant cohort.}
        \item{\code{"SAgrowth"}: Sapwood area relative growth rate (in cm2·cm-2·day-1) for an average individual of each plant cohort.}
        \item{\code{"FRAgrowth"}: Fine root area growth relative to sapwood area (in m2·cm-2·day-1) for an average individual of each plant cohort (only for \code{transpirationMode = "Sperry"}).}
    }
    \item{\code{"subdaily"}: A list of objects of class \code{\link{growth_day}}, one per day simulated (only if required in \code{control} parameters, see \code{\link{defaultControl}}).}
  }
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\note{
Objects \code{x} and \code{soil} are modified during the simulation.   
}
\seealso{
\code{\link{growthInput}}, \code{\link{growth_day}}, \code{\link{plot.growth}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control = defaultControl("Granier")

#Initialize soil with default soil params (4 layers)
examplesoil1 = soil(defaultSoilParams(4))

#Initialize vegetation input
x1 = forest2growthInput(exampleforestMED, examplesoil1, SpParamsMED, control)

#Call simulation function
G1<-growth(x1, examplemeteo, latitude = 41.82592, elevation = 100)

\dontrun{
#Switch to 'Sperry' transpiration mode
control = defaultControl("Sperry")

#Initialize soil with default soil params (4 layers)
examplesoil2 = soil(defaultSoilParams(4))

#Initialize vegetation input
x2 = forest2growthInput(exampleforestMED,examplesoil2, SpParamsMED, control)

#Call simulation function (11 days)
d = 100:110
G2<-growth(x2, examplemeteo[d,], latitude = 41.82592, elevation = 100)
}
}
