\encoding{UTF-8}
\name{biophysics}
\alias{biophysics.leafTemperature}
\alias{biophysics.radiationDiurnalPattern}
\alias{biophysics.temperatureDiurnalPattern}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Physical and biophysical utility functions
}
\description{
Set of functions used in the calculation of biophysical variables. 
}
\usage{
biophysics.leafTemperature(absRad, airTemperature, u, E, leafWidth = 0.01)
biophysics.radiationDiurnalPattern(t, daylength)
biophysics.temperatureDiurnalPattern(t, tmin, tmax, daylength)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{Wind speed above the leaf boundary layer (in m/s).}
  \item{airTemperature}{Air temperature (in ºC).}
  \item{tmin, tmax}{Minimum and maximum daily temperature (ºC).}
  \item{elevation}{Elevation above sea level (in m).}
  \item{absRad}{Absorbed long- and short-wave radiation (in W·m-2).}
  \item{E}{Transpiration flow (in mmol H20·m-2·s-1) per one sided leaf area basis.}
  \item{leafWidth}{Leaf width (in m).}
  \item{t}{Time of the day (in seconds).}
  \item{daylength}{Day length (in seconds).}
}

\value{
Values returned for each function are:
\itemize{
\item{\code{biophysics.leafTemperature}: leaf temperature (in ºC)} 
\item{\code{biophysics.radiationDiurnalPattern}: the proportion of daily radiation corresponding to the input time in seconds after sunrise.} 
\item{\code{biophysics.temperatureDiurnalPattern}: diurnal pattern of temperature assuming a sinusoidal pattern with T = Tmin at sunrise and T = (Tmin+Tmax)/2 at sunset.}
}
}
\references{
Campbell, G. S., and J. M. Norman. 1998. AN INTRODUCTION TO ENVIRONMENTAL BIOPHYSICS.: 2nd edition. (eqns. 14.1 & 14.3)

 B. Y. H. Liu and R. C. Jordan, “The interrelationship and characteristic distribution of direct, diffuse and total solar radiation,” Solar Energy, vol. 4, no. 3, pp. 1–19, 1960. 
 
McMurtrie, R. E., D. A. Rook, and F. M. Kelliher. 1990. Modelling the yield of Pinus radiata on a site limited by water and nitrogen. Forest Ecology and Management 30:381–413.
 
McMahon, T. A., M. C. Peel, L. Lowe, R. Srikanthan, and T. R. McVicar. 2013. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology & Earth System Sciences 17:1331–1363. See also:  http://www.fao.org/docrep/x0490e/x0490e06.htm
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{swb}}
}