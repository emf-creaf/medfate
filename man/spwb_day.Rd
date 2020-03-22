\encoding{UTF-8}
\name{spwb_day}
\alias{spwb_day}
\title{Soil-plant water balance for a single day}
\description{
Function \code{spwb_day} performs water balance for a single day.
}
\usage{
spwb_day(x, soil, date, tmin, tmax, rhmin, rhmax, rad, wind, 
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
  \item{latitude}{Latitude (in degrees). Required when using the 'Sperry' transpiration mode.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when using the 'Sperry' transpiration mode.}
  \item{runon}{Surface water amount running on the target area from upslope (in mm).}
}
\details{
  Detailed model description is available in the vignettes section. The model using 'Granier' transpiration mode is described in De Caceres et al. (2015). Simulations using the 'Sperry' transpiration mode are computationally much more expensive. 
}
\value{
An object (a list) of class \code{spwb_day} with the following elements:

\itemize{
  \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
  \item{\code{"EnergyBalance"}: When using the 'Sperry' transpiration mode, the model performs energy balance of the stand and 'EnergyBalance' is returned (see \code{\link{transp_transpirationSperry}}).}
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
  \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: For each leaf type, a data frame with values of LAI, Vmax298 and Jmax298 for leaves of this type in each plant cohort.}
  \item{\code{"RhizoPsi"}: Minimum water potential (in MPa) inside roots, after crossing rhizosphere, per cohort and soil layer.}
  \item{\code{"ExtractionInst"}: Water extracted by each plant cohort during each time step.}
  \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort (see \code{\link{transp_transpirationSperry}}).}
  \item{\code{"LightExtinction"}: A list of information regarding radiation balance through the canopy, as returned by function \code{\link{light_instantaneousLightExtinctionAbsortion}}.}
  \item{\code{"WindExtinction"}: A numeric vector with the wind speed (in m/s) for each plant cohort.}
  }
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology (doi:10.1016/j.agrformet.2015.06.012).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\note{
Objects \code{x} and \code{soil} are modified during the simulation. Daily transpiration and photosynthesis values are stored in columns \code{Transpiration} and \code{Photosynthesis} of object \code{x}. Water content relative to field capacity (vector \code{W}) of \code{soil} is also modified.  
}
\seealso{
\code{\link{spwbInput}}, \code{\link{spwb}},  \code{\link{plot.spwb_day}},  
\code{\link{spwb_ldrOptimization}}, \code{\link{forest}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control = defaultControl()

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2), W=c(0.5,0.5))


#Simulate one day only
x1 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
d = 100
sd1<-spwb_day(x1, examplesoil, rownames(examplemeteo)[d],  
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])

#Initialize soil with default soil params (2 layers)
examplesoil2 = soil(defaultSoilParams(2))

#Switch to 'Sperry' transpiration mode
control$transpirationMode="Sperry"


#Simulate one day only
x2 = forest2spwbInput(exampleforestMED,examplesoil2, SpParamsMED, control)
d = 100
sd2<-spwb_day(x2, examplesoil2, rownames(examplemeteo)[d],
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, slope=0, aspect=0,
             prec = examplemeteo$Precipitation[d])

#Plot plant transpiration (see function 'plot.swb.day()')
plot(sd2)

}