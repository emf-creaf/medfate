\encoding{UTF-8}
\name{transp_modes}
\alias{transp_transpirationSperry}
\alias{transp_transpirationGranier}
\title{
Transpiration modes
}
\description{
High-level sub-models to represent transpiration, plant hydraulics and water relations within plants. The two submodels represent a very different degree of complexity, and correspond to Granier et al. (1999) or Sperry et al. (2017).
}
\usage{
transp_transpirationGranier(x, meteo, day, 
                            latitude, elevation, slope, aspect, 
                            modifyInput = TRUE)
transp_transpirationSperry(x, meteo, day,
                           latitude, elevation, slope, aspect, 
                           canopyEvaporation = 0, snowMelt = 0, soilEvaporation = 0, 
                           stepFunctions = NA, 
                           modifyInput = TRUE)
}

\arguments{
  \item{x}{An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}, built using the 'Granier' or 'Sperry' transpiration modes, depending on the function to be called.}
  \item{meteo}{A data frame with daily meteorological data series:
   \itemize{
     \item{\code{DOY}: Day of the year (Julian day).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{MinTemperature}: Minimum temperature (in degrees Celsius).}
     \item{\code{MaxTemperature}: Maximum temperature (in degrees Celsius).}
     \item{\code{MinRelativeHumidity}: Minimum relative humidity (in percent).}
     \item{\code{MaxRelativeHumidity}: Maximum relative humidity (in percent).}
     \item{\code{Radiation}: Solar radiation (in MJ/m2/day).}
     \item{\code{WindSpeed}: Wind speed (in m/s). If not available, this column can be left with \code{NA} values.}
    }
  }  
  \item{day}{An integer to identify a day within \code{meteo}.}
  \item{latitude}{Latitude (in degrees).}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North).} 
  \item{canopyEvaporation}{Canopy evaporation (from interception) for \code{day} (mm).}
  \item{soilEvaporation}{Bare soil evaporation for \code{day} (mm).}
  \item{snowMelt}{Snow melt values  for \code{day} (mm).}
  \item{stepFunctions}{An integer to indicate a simulation step for which photosynthesis and profit maximization functions are desired.}
  \item{modifyInput}{Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.}
}

\value{
Function \code{transp_transpirationGranier} and \code{transp_transpirationSperry} return a list with the following elements:
\itemize{
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
  \item{\code{"Stand"}: A vector of stand-level variables.}
  \item{\code{"Plants"}: A data frame of results for each plant cohort. When using \code{transp_transpirationGranier}, element \code{"Plants"} includes:
    \itemize{
      \item{\code{"LAI"}: Leaf area index of the plant cohort.}
      \item{\code{"LAIlive"}: Leaf area index of the plant cohort, assuming all leaves are unfolded.}
      \item{\code{"AbsorbedSWRFraction"}: Fraction of SWR absorbed by each cohort.}
      \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
      \item{\code{"GrossPhotosynthesis"}: Gross photosynthesis (in gC/m2) corresponding to each cohort.}
      \item{\code{"psi"}: Water potential (in MPa) of the plant cohort (average over soil layers).}
      \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
    }
    When using \code{transp_transpirationSperry}, element \code{"Plants"} includes:
    \itemize{
      \item{\code{"LAI"}: Leaf area index of the plant cohort.}
      \item{\code{"LAIlive"}: Leaf area index of the plant cohort, assuming all leaves are unfolded.}
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
      \item{\code{"LFMC"}: Live fuel moisture content (in percent of dry weight).}
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
            \item{\code{"Net_LWR"}: A data frame with instantaneous net long-wave radiation (LWR).} 
            \item{\code{"An"}: A data frame with instantaneous net photosynthesis (in micromol/m2/s). }
            \item{\code{"Ci"}: A data frame with instantaneous intercellular CO2 concentration (in ppm). }
            \item{\code{"GW"}: A data frame with instantaneous stomatal conductance (in mol/m2/s). }
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
  \item{\code{"CanopyTurbulence"}: Canopy turbulence (see \code{\link{wind_canopyTurbulence}}).}
  \item{\code{"SupplyFunctions"}: If \code{stepFunctions} is not missing, a list of supply functions, photosynthesis functions and profit maximization functions.}
 
 }
}
\references{
Granier A, \enc{Bréda}{Breda} N, Biron P, Villette S (1999) A lumped water balance model to evaluate duration and intensity of drought constraints in forest stands. Ecol Modell 116:269–283. https://doi.org/10.1016/S0304-3800(98)00205-1.

Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40, 816-830 (doi: 10.1111/pce.12852).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{spwb_day}}, \code{\link{plot.spwb_day}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (4 layers)
examplesoil = soil(defaultSoilParams(4))

#Initialize control parameters
control = defaultControl("Granier")

#Initialize input
x1 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)

# Transpiration according to Granier's model, plant water potential 
# and plant stress for a given day
t1 = transp_transpirationGranier(x1, examplemeteo, 1, 
                                 latitude = 41.82592, elevation = 100, slope = 0, aspect = 0, 
                                 modifyInput = FALSE)

#Switch to 'Sperry' transpiration mode
control = defaultControl("Sperry")

#Initialize input
x2 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)

# Transpiration according to Sperry's model
t2 = transp_transpirationSperry(x2, examplemeteo, 1, 
                           latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
                           modifyInput = FALSE)
}
