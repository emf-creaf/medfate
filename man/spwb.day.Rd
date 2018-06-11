\encoding{UTF-8}
\name{spwb.day}
\alias{spwb.day}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil-plant water balance for a single day}
\description{
Function \code{spwb.day} performs water balance for a single day.
}
\usage{
spwb.day(x, soil, date, tmin, tmax, rhmin, rhmax, rad, wind, 
        latitude, elevation, slope, aspect, prec, runon = 0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{date}{Date as string "yyyy-mm-dd".}
  \item{tmin, tmax}{Minimum and maximum temperature (in degrees Celsius).}
  \item{rhmin, rhmax}{Minimum and maximum relative humidity (in percent).}
  \item{rad}{Solar radiation (in MJ/m2/day).}
  \item{wind}{Wind speed (in m/s).}
  \item{prec}{Precipitation (in mm).}
  \item{latitude}{Latitude (in degrees). Required when using the 'Complex' transpiration mode.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when using the 'Complex' transpiration mode.}
  \item{runon}{Surface water amount running on the target area from upslope (in mm).}
}
\details{
  Detailed model description is available in the vignettes section. The model using 'Simple' transpiration mode is described in De Caceres et al. (2015). Simulations using the 'Complex' transpiration mode are computationally much more expensive. 
}
\value{
An object (a list) of class \code{spwb.day} with the following elements:

\itemize{
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
  \item{\code{"DailyBalance"}: A list of water balance components (rain, snow, net rain, infiltration, ...) for the simulated day, equivalent to one row of 'DailyBalance' object given in \code{\link{spwb}}.}
  \item{\code{"SoilBalance"}: A list of soil water balance components (water evaporated and transpirated from each soil layer) and the current soil water potential for each soil layer.}
  \item{\code{"EnergyBalance"}: When using the 'Complex' transpiration mode, the model performs energy balance of the stand and 'EnergyBalance' is a list with the following:
    \itemize{
      \item{\code{"Temperature"}: A data frame with the temperature of the atmosphere ('Tatm'), canopy ('Tcan') and soil ('Tsoil.1', 'Tsoil.2', ...) for each time step.}
      \item{\code{"CanopyEnergyBalance"}: A data frame with the components of the canopy energy balance (in W/m2) for each time step.}
      \item{\code{"SoilEnergyBalance"}: A data frame with the components of the soil energy balance (in W/m2) for each time step.}
    }  
  }
  \item{\code{"Plants"}: A list of results for each plant cohort. When using the 'Simple' transpiration mode this includes:
    \itemize{
      \item{\code{"LAI"}: Leaf area index of the plant cohort.}
      \item{\code{"EplantCoh"}: Transpirated water (in mm) corresponding to each cohort.}
      \item{\code{"psiCoh"}: Water potential (in MPa) of the plant cohort (average over soil layers).}
      \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
    }
  When using the 'Complex' transpiration mode this includes:
    \itemize{
      \item{\code{"LAI"}: Leaf area index of the plant cohort.}
      \item{\code{"LAIsunlit"}: Leaf area index of sunlit leaves of the plant cohort.}
      \item{\code{"LAIshade"}: Leaf area index of shade leaves of the plant cohort.}
      \item{\code{"EplantCoh"}: Transpirated water (in mm) corresponding to each cohort.}
      \item{\code{"psiCoh"}: Water potential (in MPa) of the plant cohort (i.e. minimum leaf water potential).}
      \item{\code{"DDS"}: Daily drought stress [0-1] (relative whole-plant conductance).}
      \item{\code{"AbsRadinst"}: A list with four data frames containing the instantaneous absorved radiation for each plant cohort during each time step. The data frames are combinations of short-wave radiation (SWR) vs long-wave radiation (LWR) and sunlit leaves ('SL') vs. shade leaves ('SH').}
      \item{\code{"Einst"}: A data frame containing the cumulative transpiration (mm) for each plant cohort during each time step. }
      \item{\code{"Ainst"}: A data frame containing the cumulative net photosynthesis (gC/m2) for each plant cohort during each time step. }
      \item{\code{"PsiPlantinst"}: A data frame containing instantaneous leaf water potential (in MPa) for each plant cohort during each time step. }
      \item{\code{"GWsunlitinst"}: A data frame containing instantaneous stomatal conductance (in mmol/m2/s) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"GWshadeinst"}: A data frame containing instantaneous stomatal conductance (in mmol/m2/s) for shade leaves each plant cohort during each time step. }
      \item{\code{"VPDsunlitinst"}: A data frame containing instantaneous vapour pressure deficit (in kPa) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"VPDshadeinst"}: A data frame containing instantaneous vapour pressure deficit (in kPa) for shade leaves each plant cohort during each time step. }
      \item{\code{"Tempsunlitinst"}: A data frame containing instantaneous temperature (in degrees Celsius) for sunlit leaves each plant cohort during each time step. }
      \item{\code{"Tempshadeinst"}: A data frame containing instantaneous temperature (in degrees Celsius) for shade leaves each plant cohort during each time step. }
    }
  
  }
  
}
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology (doi:10.1016/j.agrformet.2015.06.012).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\note{
Objects \code{x} and \code{soil} are modified during the simulation. Daily transpiration and photosynthesis values are stored in columns \code{Transpiration} and \code{Photosynthesis} of object \code{x}. Water content relative to field capacity (vector \code{W}) of \code{soil} is also modified.  
}
\seealso{
\code{\link{spwbInput}}, \code{\link{spwb}},  \code{\link{spwbpoints}}, \code{\link{spwbgrid}}, 
\code{\link{spwb.ldrOptimization}}, \code{\link{forest}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control = defaultControl()
control$ndailysteps = 24  

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2), W=c(0.5,0.5))


#Simulate one day only
x1 = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)
d = 100
sd1<-spwb.day(x1, examplesoil, rownames(examplemeteo)[d],  
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])

#Switch to 'Complex' transpiration mode
control$transpirationMode="Complex"

#Simulate one day only
x2 = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)
d = 100
sd2<-spwb.day(x2, examplesoil, rownames(examplemeteo)[d],
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])

#Plot plant water potential (see function 'plot.swb.day()')
plot(sd2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.