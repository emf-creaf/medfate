\encoding{UTF-8}
\name{growth}
\alias{growth}
\title{Forest growth}
\description{
Function \code{growth} is a forest growth model that calculates changes in leaf area, sapwood area and structural variables for all plant cohorts in a given forest stand during a period specified in the input climatic data. 
}
\usage{
growth(x, soil, meteo,  latitude, elevation = NA, slope = NA, aspect = NA) 
}
\arguments{
  \item{x}{An object of class \code{\link{growthInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
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
  \item{\code{"spwbInput"}: A copy of the object \code{x} of class \code{\link{spwbInput}} given as input.}
  \item{\code{"soilInput"}: A copy of the object \code{soil} of class \code{\link{soil}} given as input.}
  \item{\code{"WaterBalance"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
  \itemize{
  \item{\code{"PET"}: Potential evapotranspiration (in mm).}
  \item{\code{"Rainfall"}: Input precipitation (in mm).}
  \item{\code{"NetPrec"}: Net precipitation, after accounting for interception (in mm).}
  \item{\code{"Infiltration"}: The amount of water infiltrating into the soil (in mm).}
  \item{\code{"Runoff"}: The amount of water exported via surface runoff (in mm).}
  \item{\code{"DeepDrainage"}: The amount of water exported via deep drainage (in mm).}
  \item{\code{"Evapotranspiration"}: Evapotranspiration (in mm).}
  \item{\code{"SoilEvaporation"}: Bare soil evaporation (in mm).}
  \item{\code{"Transpiration"}: Plant transpiration (considering all soil layers) (in mm).}
  }
  \item{\code{"Soil"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
  \itemize{
  \item{\code{"W.1"}, \code{...}, \code{"W.k"}: Relative soil moisture content (relative to field capacity) in each soil layer.}
  \item{\code{"ML.1"}, \code{...}, \code{"ML.k"}: Soil water volume in each soil layer (in L/m2).}
  \item{\code{"MLTot"}: Total soil water volume (in L/m2).}
  \item{\code{"PlantExt.1"}, \code{...}, \code{"PlantExt.k"}: Plant extraction from each soil layer (in mm).}
  \item{\code{"psi.1"}, \code{...}, \code{"psi.k"}: Soil water potential in each soil layer (in MPa).}
  }
  \item{\code{"Stand"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
       \itemize{
          \item{\code{"GDD"}: Cumulative growth degree days.}
          \item{\code{"LAIcell"}: LAI of the stand (accounting for leaf phenology) (in m2/m2).}
          \item{\code{"LAIcelldead"}: LAI of the stand corresponding to dead leaves (in m2/m2).}
          \item{\code{"Cm"}: Water retention capacity of the canopy (in mm) (accounting for leaf phenology).}
          \item{\code{"LgroundPAR"}: The proportion of PAR that reaches the ground (accounting for leaf phenology).}
          \item{\code{"LgroundSWR"}: The proportion of SWR that reaches the ground (accounting for leaf phenology).}
      }
  \item{\code{"PlantTranspiration"}: A data frame with the amount of daily transpiration (in mm) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantGrossPhotosynthesis"}: A data frame with the amount of daily gross photosynthesis (in g C · m-2) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantRespiration"}: A data frame with the amount of daily maintenance respiration (in g C · m-2) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantCStorageFast"}: A data frame with the daily amount of fast-dynamics carbon reserves (as proportion of total storage capacity) for an average individual of each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantCStorageSlow"}: A data frame with the daily amount of slow-dynamics carbon reserves (as proportion of total storage capacity) for an average individual of each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantSAgrowth"}: A data frame with the daily amount of newly created sapwood area (in cm2) for an average individual of each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantSA"}: A data frame with the daily amount of sapwood area (in cm2) for an average individual of each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantStress"}: A data frame with the amount of daily stress suffered by each plant cohort (relative whole-plant conductance). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"PlantPsi"}: A data frame with the average daily water potential of each plant (in MPa). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"PlantLAIdead"}: A data frame with the dead leaf area index (in m2 · m-2) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"PlantLAIlive"}: A data frame with the live leaf area index (in m2 · m-2) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
 \item{\code{"StandStructures"}: A list with the aboveground structure of the stand at the beginning of the simulation and after each simulated year.}
 \item{\code{"StandSummary"}: A data frame with stand-level summaries (leaf area index, tree basal area, tree density, shrub cover, etc.) at the beginning of the simulation and after each simulated year.}
  }
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
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

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl()

#Initialize input
x = forest2growthInput(exampleforestMED,examplesoil, SpParamsMED, control)

#Call simulation function
G1<-growth(x, examplesoil, examplemeteo, latitude = 41.82592, elevation = 100)
 
}