\encoding{UTF-8}
\name{transpiration}
\alias{transp_profitMaximization}
\alias{transp_stomatalRegulationPlot}
\alias{transp_transpirationSperry}
\alias{transp_transpirationGranier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transpiration submodel functions
}
\description{
Set of high-level functions used in the calculation of stomatal conductance and transpiration. Function \code{transp_profitMaximization} calculates gain and cost functions, as well as profit maximization from supply and photosynthesis input functions. Function \code{transp_stomatalRegulationPlot} produces a plot with the cohort supply functions against water potential and a plot with the cohort photosynthesis functions against water potential, both with the maximum profit values indicated.
}
\usage{
transp_profitMaximization(supplyFunction, photosynthesisFunction, type, 
                          Gwmin, Gwmax, kleafmax = NA)
transp_transpirationGranier(x, soil, meteo, day, modifyInput = TRUE)
transp_transpirationSperry(x, soil, meteo, day,
                           latitude, elevation, 
                           canopyEvaporation = 0, snowMelt = 0, soilEvaporation = 0, 
                           stepFunctions = NA, modifyInput = TRUE)
transp_stomatalRegulationPlot(x, soil, meteo, day, timestep, 
                              latitude, elevation)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{supplyFunction}{Water supply function (see \code{\link{hydraulics_supplyFunctionNetwork}}).}
  \item{photosynthesisFunction}{Function returned by \code{photo_photosynthesisFunction()}.}
  \item{type}{Optimization type (1,2 or 3).}
  \item{Gwmin, Gwmax}{Minimum and maximum stomatal conductance (mol·m-2·s-1).}
  \item{kleafmax}{Maximum leaf hydraulic conductance (flow per pressure drop).}
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
  \item{elevation}{Elevation above sea level (in m).}
  \item{canopyEvaporation}{Canopy evaporation (from interception) for \code{day} (mm).}
  \item{soilEvaporation}{Bare soil evaporation for \code{day} (mm).}
  \item{snowMelt}{Snow melt values  for \code{day} (mm).}
  \item{stepFunctions}{An integer to indicate a simulation step for which photosynthesis and profit maximization functions are desired.}
  \item{modifyInput}{Boolean flag to indicate that the input \code{x} and \code{soil} objects are allowed to be modified during the simulation.}
}
\details{
Details of the transpiration submodel are given in a vignette.
}

\value{
Function \code{transp_transpirationGranier} and \code{transp_transpirationSperry} return a list with the following elements:
\itemize{
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
  \item{\code{"Plants"}: A data frame of results for each plant cohort:
    \itemize{
      \item{\code{"LAI"}: Leaf area index of the plant cohort.}
      \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
      \item{\code{"psi"}: Water potential (in MPa) of the plant cohort (average over soil layers).}
      \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
    }
    When using \code{transp_transpirationSperry}, element \code{"Plants"} includes:
    \itemize{
      \item{\code{"Extraction"}: Water extracted from the soil (in mm) for each cohort.}
      \item{\code{"RootPsi"}: Minimum water potential (in MPa) at the root collar.}
      \item{\code{"StemPsi"}: Minimum water potential (in MPa) at the stem.}
      \item{\code{"LeafPsi"}: Minimum water potential (in MPa) at the leaf.}
      \item{\code{"StemPLC"}: Proportion of conductance loss in stem.}
      \item{\code{"StemRWC"}: Relative water content of symplastic stem tissue.}
      \item{\code{"LeafRWC"}: Relative water content of symplastic leaf tissue.}
      \item{\code{"dEdP"}: Overall soil-plant conductance (derivative of the supply function).}
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
   \item{\code{"ExtractionInst"}: Water extracted by each plant cohort during each time step.}
   \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort:
    \itemize{
      \item{\code{"LAIsunlit"}: Leaf area index of sunlit leaves of the plant cohort.}
      \item{\code{"LAIshade"}: Leaf area index of shade leaves of the plant cohort.}
      \item{\code{"AbsRad"}: A list with four data frames containing the instantaneous absorved radiation for each plant cohort during each time step. The data frames are combinations of short-wave radiation (SWR) vs long-wave radiation (LWR) and sunlit leaves ('SL') vs. shade leaves ('SH').}
      \item{\code{"E"}: A data frame containing the cumulative transpiration (mm) for each plant cohort during each time step. }
      \item{\code{"An"}: A data frame containing the cumulative net photosynthesis (gC/m2) for each plant cohort during each time step. }
      \item{\code{"GWsunlit"}: A data frame containing instantaneous stomatal conductance (in mmol/m2/s) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"GWshade"}: A data frame containing instantaneous stomatal conductance (in mmol/m2/s) for shade leaves each plant cohort during each time step. }
      \item{\code{"VPDsunlit"}: A data frame containing vapour pressure deficit (in kPa) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"VPDshade"}: A data frame containing vapour pressure deficit (in kPa) for shade leaves each plant cohort during each time step. }
      \item{\code{"Tempsunlit"}: A data frame containing temperature (in degrees Celsius) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"Tempshade"}: A data frame containing temperature (in degrees Celsius) for shade leaves each plant cohort during each time step. }
      \item{\code{"PsiRoot"}: A data frame containing root crown water potential (in MPa) for each plant cohort during each time step.}
      \item{\code{"PsiPlant"}: A data frame containing leaf water potential (in MPa) for each plant cohort during each time step. }
      \item{\code{"PLCstem"}: A data frame containing the proportion loss of conductance [0-1] for each plant cohort during each time step. }
      \item{\code{"RWCstem"}: A data frame containing the (average) relative water content of symplastic stem tissue [0-1] for each plant cohort during each time step. }
      \item{\code{"RWCstem"}: A data frame containing the relative water content of symplastic leaf tissue [0-1] for each plant cohort during each time step. }
    }
  }
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
\code{\link{hydraulics_supplyFunctionNetwork}}, \code{\link{biophysics_leafTemperature}}, \code{\link{photo_photosynthesis}}, \code{\link{spwb}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)


#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl()

#Initialize input
x = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)

# Transpiration according to Granier's model, plant water potential 
# and plant stress for a given day
transp_transpirationGranier(x, examplesoil, examplemeteo, 1, 
                            modifyInput = FALSE)

#Switch to 'Sperry' transpiration mode
control$transpirationMode="Sperry"

#Initialize soil with default soil params (2 layers)
examplesoil2 = soil(defaultSoilParams(2))

#Initialize input
x2 = forest2spwbInput(exampleforest,examplesoil2, SpParamsMED, control)

# Stomatal regulation for the 12th time step at day 100
transp_stomatalRegulationPlot(x2, examplesoil2, examplemeteo, day=100, timestep = 12,
                              latitude = 41.82592, elevation = 100)
}