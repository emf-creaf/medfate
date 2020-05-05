\encoding{UTF-8}
\name{transpiration}
\alias{transp_profitMaximization}
\alias{transp_stomatalRegulationPlot}
\alias{transp_transpirationSperry}
\alias{transp_transpirationGranier}
\title{
Transpiration submodel functions
}
\description{
Set of high-level functions used in the calculation of stomatal conductance and transpiration. Function \code{transp_profitMaximization} calculates gain and cost functions, as well as profit maximization from supply and photosynthesis input functions. Function \code{transp_stomatalRegulationPlot} produces a plot with the cohort supply functions against water potential and a plot with the cohort photosynthesis functions against water potential, both with the maximum profit values indicated.
}
\usage{
transp_profitMaximization(supplyFunction, photosynthesisFunction, Gwmin, Gwmax, 
                          gainModifier = 1.0, costModifier = 1.0, costWater = "dEdP")
transp_transpirationGranier(x, soil, meteo, day, 
                            modifyInputX = TRUE, modifyInputSoil = TRUE)
transp_transpirationSperry(x, soil, meteo, day,
                           latitude, elevation, slope, aspect, 
                           canopyEvaporation = 0, snowMelt = 0, soilEvaporation = 0, 
                           stepFunctions = NA, 
                           modifyInputX = TRUE, modifyInputSoil = TRUE)
transp_stomatalRegulationPlot(x, soil, meteo, day, timestep, 
                              latitude, elevation, slope = NA, aspect = NA,
                              type="E")
}

\arguments{
  \item{supplyFunction}{Water supply function (see \code{\link{hydraulics_supplyFunctionNetwork}}).}
  \item{photosynthesisFunction}{Function returned by \code{photo_photosynthesisFunction()}.}
  \item{Gwmin, Gwmax}{Minimum and maximum stomatal conductance (mol·m-2·s-1).}
  \item{gainModifier, costModifier}{Modifiers (exponents) of the gain and cost functions defined in Sperry et al. (2016).}
  \item{costWater}{Cost water function, either "dEdP" (default) or "kterm" (experimental).}
  \item{x}{An object of class \code{\link{spwbInput}} built using the 'Sperry' transpiration mode.}
  \item{soil}{An object of class \code{\link{soil}}.}
  \item{meteo}{A data frame with daily meteorological data series:
   \itemize{
     \item{\code{DOY}: Day of the year (Julian day).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{MinTemperature}: Minimum temperature (in degrees Celsius).}
     \item{\code{MaxTemperature}: Maximum temperature (in degrees Celsius).}
     \item{\code{MinRelativeHumidity}: Minimum relative humidity (in percent).}
     \item{\code{MaxRelativeHumidity}: Maximum relative humidity (in percent).}
     \item{\code{Radiation}: Solar radiation (in MJ/m2/day).}
     \item{\code{WindSpeed}: Wind speed (in m/s). If not available, this column can be left with \code{NA} values.}
    }
  }  
  \item{day}{An integer to identify a day within \code{meteo}.}
  \item{timestep}{An integer between 1 and \code{ndailysteps} specified in \code{x} (see \code{\link{defaultControl}}).}
  \item{latitude}{Latitude (in degrees).}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North).} 
  \item{canopyEvaporation}{Canopy evaporation (from interception) for \code{day} (mm).}
  \item{soilEvaporation}{Bare soil evaporation for \code{day} (mm).}
  \item{snowMelt}{Snow melt values  for \code{day} (mm).}
  \item{stepFunctions}{An integer to indicate a simulation step for which photosynthesis and profit maximization functions are desired.}
  \item{type}{A string with plot type, either \code{"E"} (transpiration flow), \code{"Ag"} (gross photosynthesis), \code{"An"} (net photosynthesis), \code{"Gw"} (stomatal conductance), \code{"T"} (temperature) or \code{"VPD"} (vapour pressure deficit).}
  \item{modifyInputX, modifyInputSoil}{Boolean flags to indicate that the input \code{x} and \code{soil} objects are allowed to be modified during the simulation.}
}

\value{
Function \code{transp_transpirationGranier} and \code{transp_transpirationSperry} return a list with the following elements:
\itemize{
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
  \item{\code{"Stand"}: A vector of stand-level variables.}
  \item{\code{"Plants"}: A data frame of results for each plant cohort. When using \code{transp_transpirationGranier}, element \code{"Plants"} includes:
    \itemize{
      \item{\code{"LAI"}: Leaf area index of the plant cohort.}
      \item{\code{"AbsorbedSWRFraction"}: Fraction of SWR absorbed by each cohort.}
      \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
      \item{\code{"Photosynthesis"}: Net photosynthesis (in gC/m2) corresponding to each cohort.}
      \item{\code{"psi"}: Water potential (in MPa) of the plant cohort (average over soil layers).}
      \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
    }
    When using \code{transp_transpirationSperry}, element \code{"Plants"} includes:
    \itemize{
      \item{\code{"LAI"}: Leaf area index of the plant cohort.}
      \item{\code{"Extraction"}: Water extracted from the soil (in mm) for each cohort.}
      \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
      \item{\code{"GrossPhotosynthesis"}: Gross photosynthesis (in gC/m2) corresponding to each cohort.}
      \item{\code{"NetPhotosynthesis"}: Net photosynthesis (in gC/m2) corresponding to each cohort.}
      \item{\code{"RootPsi"}: Minimum water potential (in MPa) at the root collar.}
      \item{\code{"StemPsi"}: Minimum water potential (in MPa) at the stem.}
      \item{\code{"StemPLC"}: Proportion of conductance loss in stem.}
      \item{\code{"LeafPsiMin"}: Minimum (predawn) water potential (in MPa) at the leaf (representing an average leaf).}
      \item{\code{"LeafPsiMax"}: Maximum (midday) water potential (in MPa) at the leaf (representing an average leaf).}
      \item{\code{"LeafPsiMin_SL"}: Minimum (predawn) water potential (in MPa) at sunlit leaves.}
      \item{\code{"LeafPsiMax_SL"}: Maximum (midday) water potential (in MPa) at sunlit leaves.}
      \item{\code{"LeafPsiMin_SH"}: Minimum (predawn) water potential (in MPa) at shade leaves.}
      \item{\code{"LeafPsiMax_SH"}: Maximum (midday) water potential (in MPa) at shade leaves.}
      \item{\code{"dEdP"}: Overall soil-plant conductance (derivative of the supply function).}
      \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
      \item{\code{"StemRWC"}: Relative water content of stem tissue (including symplasm and apoplasm).}
      \item{\code{"LeafRWC"}: Relative water content of leaf tissue (including symplasm and apoplasm).}
      \item{\code{"StemSympRWC"}: Relative water content of symplastic stem tissue.}
      \item{\code{"LeafSympRWC"}: Relative water content of symplastic leaf tissue.}
      \item{\code{"WaterBalance"}: Plant water balance (extraction - transpiration).}
    }
  }
  \item{\code{"Extraction"}: A data frame with mm of water extracted from each soil layer (in columns) by each cohort (in rows).}
  
  The remaining items are only given by \code{transp_transpirationSperry}:
  \item{\code{"EnergyBalance"}: When using the 'Sperry' transpiration mode, the model performs energy balance of the stand and 'EnergyBalance' is a list with the following:
    \itemize{
      \item{\code{"Temperature"}: A data frame with the temperature of the atmosphere ('Tatm'), canopy ('Tcan') and soil ('Tsoil.1', 'Tsoil.2', ...) for each time step.}
      \item{\code{"CanopyEnergyBalance"}: A data frame with the components of the canopy energy balance (in W/m2) for each time step.}
      \item{\code{"SoilEnergyBalance"}: A data frame with the components of the soil energy balance (in W/m2) for each time step.}
    }  
   }
   \item{\code{"RhizoPsi"}: Minimum water potential (in MPa) inside roots, after crossing rhizosphere, per cohort and soil layer.}
   \item{\code{"Sunlitleaves"} and \code{"ShadeLeaves"}: Data frames for sunlit leaves and shade leaves and the following columns per cohort:
      \itemize{
      \item{\code{"LAI"}: Cumulative leaf area index of sunlit/shade leaves.}
      \item{\code{"Vmax298"}: Average maximum carboxilation rate for sunlit/shade leaves.}
      \item{\code{"Jmax298"}: Average maximum electron transport rate for sunlit/shade leaves.}
    }  
   }
   \item{\code{"ExtractionInst"}: Water extracted by each plant cohort during each time step.}
   \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort:
    \itemize{
      \item{\code{"E"}: A data frame with the cumulative transpiration (mm) for each plant cohort during each time step. }
      \item{\code{"Ag"}: A data frame with the cumulative gross photosynthesis (gC/m2) for each plant cohort during each time step. }
      \item{\code{"An"}: A data frame with the cumulative net photosynthesis (gC/m2) for each plant cohort during each time step. }
      \item{\code{"Sunlitleaves"} and \code{"ShadeLeaves"}: Lists with instantaneous (for each time step) results for sunlit leaves and shade leaves and the following items:
        \itemize{
            \item{\code{"Abs_SWR"}: A data frame with instantaneous absorbed short-wave radiation (SWR).} 
            \item{\code{"Abs_LWR"}: A data frame with instantaneous absorbed long-wave radiation (LWR).} 
            \item{\code{"An"}: A data frame with instantaneous net photosynthesis (in micromol/m2/s). }
            \item{\code{"Ci"}: A data frame with instantaneous intercellular CO2 concentration (in ppm). }
            \item{\code{"GW"}: A data frame with instantaneous stomatal conductance (in mmol/m2/s). }
            \item{\code{"VPD"}: A data frame with instantaneous vapour pressure deficit (in kPa). }
            \item{\code{"Temp"}: A data frame with leaf temperature (in degrees Celsius). }
            \item{\code{"Psi"}: A data frame with leaf water potential (in MPa). }
        }
      }
      \item{\code{"dEdP"}: A data frame with the slope of the plant supply function (an estimation of whole-plant conductance).}
      \item{\code{"RootPsi"}: A data frame with root crown water potential (in MPa) for each plant cohort during each time step.}
      \item{\code{"StemPsi"}: A data frame with stem water potential (in MPa) for each plant cohort during each time step.}
      \item{\code{"LeafPsi"}: A data frame with leaf (average) water potential (in MPa) for each plant cohort during each time step. }
      \item{\code{"StemPLC"}: A data frame with the proportion loss of conductance [0-1] for each plant cohort during each time step. }
      \item{\code{"StemRWC"}: A data frame with the (average) relative water content of stem tissue [0-1] for each plant cohort during each time step. }
      \item{\code{"LeafRWC"}: A data frame with the relative water content of leaf tissue [0-1] for each plant cohort during each time step. }
      \item{\code{"StemSympRWC"}: A data frame with the (average) relative water content of symplastic stem tissue [0-1] for each plant cohort during each time step. }
      \item{\code{"LeafSympRWC"}: A data frame with the relative water content of symplastic leaf tissue [0-1] for each plant cohort during each time step. }
      \item{\code{"PWB"}: A data frame with plant water balance (extraction - transpiration).}
    }
  }
  \item{\code{"LightExtinction"}: A list of information regarding radiation balance through the canopy, as returned by function \code{\link{light_instantaneousLightExtinctionAbsortion}}.}
  \item{\code{"WindExtinction"}: A numeric vector with the wind speed (in m/s) for each plant cohort.}
  \item{\code{"SupplyFunctions"}: If \code{stepFunctions} is not missing, a list of supply functions, photosynthesis functions and profit maximization functions.}
 
 }

Function \code{transp_profitMaximization} returns a list with the following elements:
\itemize{
  \item{\code{Cost}: Cost function [0-1].}
  \item{\code{Gain}: Gain function [0-1].}
  \item{\code{Profit}: Profit function [0-1].}
  \item{\code{iMaxProfit}: Index corresponding to maximum profit (starting from 0).}
}

}
\references{
Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2016. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{hydraulics_supplyFunctionNetwork}}, \code{\link{biophysics_leafTemperature}}, \code{\link{photo_photosynthesis}}, \code{\link{spwb_day}}, \code{\link{plot.spwb_day}}
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
x1 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)

# Transpiration according to Granier's model, plant water potential 
# and plant stress for a given day
t1 = transp_transpirationGranier(x1, examplesoil, examplemeteo, 1, 
                            modifyInputX = FALSE, modifyInputSoil = FALSE)

#Switch to 'Sperry' transpiration mode
control$transpirationMode="Sperry"

#Initialize soil with default soil params (2 layers)
examplesoil2 = soil(defaultSoilParams(2))

#Initialize input
x2 = forest2spwbInput(exampleforestMED,examplesoil2, SpParamsMED, control)

# Transpiration according to Sperry's model
t2 = transp_transpirationSperry(x2, examplesoil2, examplemeteo, 1, 
                           latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
                           modifyInputX = FALSE, modifyInputSoil = FALSE)
                            
# Stomatal VPD curve and chosen value for the 12th time step at day 100
transp_stomatalRegulationPlot(x2, examplesoil2, examplemeteo, day=100, timestep = 12,
                              latitude = 41.82592, elevation = 100, type="VPD")
}