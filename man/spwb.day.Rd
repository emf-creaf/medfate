\encoding{UTF-8}
\name{spwb.day}
\alias{spwb.day}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil water balance for a single day}
\description{
Function \code{spwb.day} performs soil water balance for a single day.
}
\usage{
spwb.day(x, soil, date, doy, tmin, tmax, rhmin, rhmax, rad, wind, 
        latitude, elevation, slope, aspect, rain, runon = 0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{date}{Date as string "yyyy-mm-dd".}
  \item{doy}{Day of the year.}
  \item{tmin, tmax}{Minimum and maximum temperature (in degrees Celsius).}
  \item{rhmin, rhmax}{Minimum and maximum relative humidity (in percent).}
  \item{rad}{Solar radiation (in MJ/m2/day).}
  \item{wind}{Wind speed (in m/s).}
  \item{rain}{Precipitation (in mm).}
  \item{runon}{Surface water amount running on the target area from upslope (in mm).}
  \item{latitude}{Latitude (in degrees). Required when using the 'Complex' transpiration mode.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when using the 'Complex' transpiration mode.}
}
\details{
  Detailed model description is available in the vignettes section. The model using 'Simple' transpiration mode is described in De Caceres et al. (2015). Simulations using the 'Complex' transpiration mode are computationally much more expensive. 
}
\value{
An object of class \code{spwb.day}:
TO BE DONE
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
sd1<-spwb.day(x1, examplesoil, rownames(examplemeteo)[d], d , 
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, 
             slope= 0, aspect = 0, rain = examplemeteo$Precipitation[d])

#Switch to 'Complex' transpiration mode
control$transpirationMode="Complex"

#Simulate one day only
x2 = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)
d = 100
sd2<-spwb.day(x2, examplesoil, rownames(examplemeteo)[d], d , 
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, 
             slope= 0, aspect = 0, rain = examplemeteo$Precipitation[d])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.