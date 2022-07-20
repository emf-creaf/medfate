\encoding{UTF-8}
\name{spwb_day}
\alias{spwb_day}
\alias{growth_day}
\title{Single day simulation}
\description{
Function \code{spwb_day} performs water balance for a single day and \code{growth_day} performs water and carbon balance for a single day.
}
\usage{
spwb_day(x, date, tmin, tmax, rhmin, rhmax, rad, wind, 
         latitude, elevation, slope, aspect, prec, CO2 = NA, 
         runon = 0.0, modifyInput = TRUE)
growth_day(x, date, tmin, tmax, rhmin, rhmax, rad, wind, 
           latitude, elevation, slope, aspect, prec, CO2 = NA,
           runon = 0.0, modifyInput = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.}
  \item{date}{Date as string "yyyy-mm-dd".}
  \item{tmin, tmax}{Minimum and maximum temperature (in degrees Celsius).}
  \item{rhmin, rhmax}{Minimum and maximum relative humidity (in percent).}
  \item{rad}{Solar radiation (in MJ/m2/day).}
  \item{wind}{Wind speed (in m/s).}
  \item{prec}{Precipitation (in mm).}
  \item{latitude}{Latitude (in degrees). Required when using the 'Sperry' transpiration mode.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when using the 'Sperry' transpiration mode.}
  \item{CO2}{Atmospheric CO2 concentration (in ppm). If missing, default value is drawn from control parameter 'defaultCO2' in \code{x}.}
  \item{runon}{Surface water amount running on the target area from upslope (in mm).}
  \item{modifyInput}{Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.}
}
\details{
  Detailed model description is available in the vignettes section. The model using 'Granier' transpiration mode is described in De Caceres et al. (2015). Simulations using the 'Sperry' transpiration mode are computationally much more expensive, are described in De Cáceres et al. (2021) and are illustrated by function \code{\link{transp_transpirationSperry}}. 
}
\value{
Function \code{spwb_day()} returns a list of class \code{spwb_day} with the following elements:

\itemize{
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
  \item{\code{"topography"}: Vector with elevation, slope and aspect given as input.} 
  \item{\code{"weather"}: A vector with the input weather.}
  \item{\code{"WaterBalance"}: A vector of water balance components (rain, snow, net rain, infiltration, ...) for the simulated day, equivalent to one row of 'WaterBalance' object given in \code{\link{spwb}}.}
  \item{\code{"Soil"}: A data frame with results for each soil layer:
    \itemize{
      \item{\code{"SoilEvaporation"}: Water evaporated from the soil surface (in mm).}
      \item{\code{"HydraulicInput"}: Water entering each soil layer from other layers, transported via plant hydraulic network (in mm) (only for \code{transpirationMode = "Sperry"}).}
      \item{\code{"HydraulicOutput"}: Water leaving each soil layer (going to other layers or the transpiration stream) (in mm) (only for \code{transpirationMode = "Sperry"}).}
      \item{\code{"PlantExtraction"}: Water extracted by plants from each soil layer (in mm).}
      \item{\code{"psi"}: Soil water potential (in MPa).}
    }
  }
  \item{\code{"Stand"}: A named vector with with stand values for the simulated day, equivalent to one row of 'Stand' object returned by \code{\link{spwb}}.}
  \item{\code{"Plants"}: A data frame of results for each plant cohort (see \code{\link{transp_transpirationGranier}} or \code{\link{transp_transpirationSperry}}).}
}
The following items are only returned when \code{transpirationMode = "Sperry"}:
\itemize{
  \item{\code{"EnergyBalance"}: Energy balance of the stand (see \code{\link{transp_transpirationSperry}}).}
  \item{\code{"RhizoPsi"}: Minimum water potential (in MPa) inside roots, after crossing rhizosphere, per cohort and soil layer.}
  \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: For each leaf type, a data frame with values of LAI, Vmax298 and Jmax298 for leaves of this type in each plant cohort.}
  \item{\code{"ExtractionInst"}: Water extracted by each plant cohort during each time step.}
  \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort (see \code{\link{transp_transpirationSperry}}).}
  \item{\code{"LightExtinction"}: A list of information regarding radiation balance through the canopy, as returned by function \code{\link{light_instantaneousLightExtinctionAbsortion}}.}
  \item{\code{"CanopyTurbulence"}: Canopy turbulence (see \code{\link{wind_canopyTurbulence}}).}
  }
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology (doi:10.1016/j.agrformet.2015.06.012).

De \enc{Cáceres}{Caceres} M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L, Poyatos R, Cabon A, Granda V, Forner A, Valladares F, \enc{Martínez}{Martinez}-Vilalta J (2021) Unravelling the effect of species mixing on water use and drought stress in holm oak forests: a modelling approach. Agricultural and Forest Meteorology 296 (doi: 10.1016/j.agrformet.2020.108233).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{spwbInput}}, \code{\link{spwb}},  \code{\link{plot.spwb_day}},  
\code{\link{growthInput}}, \code{\link{growth}},  \code{\link{plot.growth_day}}  
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

# Day to be simulated
d = 100


#Simulate water balance one day only (Granier)
examplesoil = soil(defaultSoilParams(4))
x1 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
sd1<-spwb_day(x1, rownames(examplemeteo)[d],  
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])


#Simulate water balance for one day only (Sperry's mode)
control = defaultControl("Sperry")
x2 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
sd2<-spwb_day(x2, rownames(examplemeteo)[d],
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])

#Plot plant transpiration (see function 'plot.swb.day()')
plot(sd2)


#Simulate water and carbon balance for one day only
x3 = forest2growthInput(exampleforestMED,examplesoil, SpParamsMED, control)
sd3<-growth_day(x3, rownames(examplemeteo)[d],
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])

}
