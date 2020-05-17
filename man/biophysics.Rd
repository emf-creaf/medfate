\encoding{UTF-8}
\name{biophysics}
\alias{biophysics_leafTemperature}
\alias{biophysics_radiationDiurnalPattern}
\alias{biophysics_temperatureDiurnalPattern}
\alias{biophysics_waterDynamicViscosity}

\title{
Physical and biophysical utility functions
}
\description{
Utility functions for the calculation of biophysical variables. 
}
\usage{
biophysics_leafTemperature(absRad, airTemperature, u, E, leafWidth = 1.0)
biophysics_radiationDiurnalPattern(t, daylength)
biophysics_temperatureDiurnalPattern(t, tmin, tmax, 
                                     tminPrev, tmaxPrev, tminNext, daylength)
biophysics_waterDynamicViscosity(temp)
}

\arguments{
  \item{u}{Wind speed above the leaf boundary layer (in m/s).}
  \item{airTemperature}{Air temperature (in ºC).}
  \item{tmin, tmax}{Minimum and maximum daily temperature (ºC).}
  \item{tminPrev, tmaxPrev, tminNext}{Maximum and minimum daily temperatures of the previous and following day (ºC).}
  \item{absRad}{Absorbed long- and short-wave radiation (in W·m-2).}
  \item{E}{Transpiration flow (in mmol H20·m-2·s-1) per one sided leaf area basis.}
  \item{leafWidth}{Leaf width (in cm).}
  \item{t}{Time of the day (in seconds).}
  \item{daylength}{Day length (in seconds).}
  \item{temp}{Temperature (ºC).}
}
\details{
Function \code{biophysics_leafTemperature} calculates leaf temperature according to energy balance equation given in Campbell and Norman (1988). Function \code{biophysics_radiationDiurnalPattern} follows the equations given in Liu and Jordan (1960). Function \code{biophysics_temperatureDiurnalPattern} determines diurnal temperature pattern assuming a sinusoidal pattern with T = Tmin at sunrise and T = (Tmin+Tmax)/2 at sunset and a linear change in temperature between sunset and Tmin of the day after (McMurtrie et al. 1990). Function \code{biophysics_waterDynamicViscosity} calculates water dynamic viscosity following the Vogel (1921) equation.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{biophysics_leafTemperature}: leaf temperature (in ºC)} 
\item{\code{biophysics_radiationDiurnalPattern}: the proportion of daily radiation corresponding to the input time in seconds after sunrise.} 
\item{\code{biophysics_temperatureDiurnalPattern}: diurnal pattern of temperature.}
\item{\code{biophysics_waterDynamicViscosity}: Water dynamic viscosity relative to 20ºC.} 
}
}
\references{
Campbell, G. S., and J. M. Norman. 1998. An introduction to environmental biophysics: 2nd edition. (eqns. 14.1 & 14.3)

B. Y. H. Liu and R. C. Jordan, “The interrelationship and characteristic distribution of direct, diffuse and total solar radiation,” Solar Energy, vol. 4, no. 3, pp. 1–19, 1960. 
 
McMurtrie, R. E., D. A. Rook, and F. M. Kelliher. 1990. Modelling the yield of Pinus radiata on a site limited by water and nitrogen. Forest Ecology and Management 30:381–413.

H. Vogel, "Das Temperaturabhangigkeitsgesetz der Viskositat von Flussigkeiten", Physikalische Zeitschrift, vol. 22, pp. 645–646, 1921.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{spwb}}
}