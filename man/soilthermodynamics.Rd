\encoding{UTF-8}
\name{soil thermodynamics}
\alias{soil.thermalconductivity}
\alias{soil.thermalcapacity}
\alias{soil.temperaturechange}
\alias{soil.temperaturegradient}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil thermodynamic functions}
\description{
Functions \code{soil.thermalconductivity} and \code{soil.thermalcapacity} calculate thermal conductivity and thermal capacity for each soil layer, given its texture and water content. Functions \code{soil.temperaturegradient} and \code{soil.temperaturechange} are used to calculate soil temperature gradients (in ºC/m) and temporal temperature change (in ºC/s) given soil layer texture and water content (and possibly including heat flux from above).
}
\usage{
soil.thermalconductivity(soil, model = "SX")
soil.thermalcapacity(soil, model = "SX")
soil.temperaturechange(dVec, Temp, sand, clay, W, Theta_FC, Gdown)
soil.temperaturegradient(dVec, Temp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{soil}{Soil object (returned by function \code{\link{soil}}).}
  \item{model}{Either 'SX' or 'VG' for Saxton's or Van Genuchten's pedotransfer models.}
  \item{dVec}{Width of soil layers (in mm).}
  \item{Temp}{Temperature (in ºC) for each soil layer.}
  \item{clay}{Percentage of clay (in percent weight) for each layer.}
  \item{sand}{Percentage of sand (in percent weight) for each layer.}
  \item{W}{Soil moisture (in percent of field capacity) for each layer.}
  \item{Theta_FC}{Relative water content (in percent volume) at field capacity for each layer.}
  \item{Gdown}{Downward heat flux from canopy to soil (in W·m-2).}
}
\value{
Function \code{soil.thermalconductivity} returns a vector with values of thermal conductivity (W/m/ºK) for each soil layer. Function \code{soil.thermalcapacity} returns a vector with values of heat storage capacity (J/m3/ºK) for each soil layer. Function \code{soil.temperaturegradient} returns a vector with values of temperature gradient between consecutive soil layers. Function \code{soil.temperaturechange} returns a vector with values of instantaneous temperature change (ºC/s) for each soil layer.
}
\references{
Cox, P.M., Betts, R.A., Bunton, C.B., Essery, R.L.H., Rowntree, P.R., & Smith, J. 1999. The impact of new land surface physics on the GCM simulation of climate and climate sensitivity. Climate Dynamics 15: 183–203.

Dharssi, I., Vidale, P.L., Verhoef, A., MacPherson, B., Jones, C., & Best, M. 2009. New soil physical properties implemented in the Unified Model at PS18. 9–12.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{ \code{\link{soil}}}
\examples{
examplesoil = soil(defaultSoilParams())
soil.thermalconductivity(examplesoil)
soil.thermalcapacity(examplesoil)

#Values change when altering water content (drier layers have lower conductivity and capacity)
examplesoil$W = c(0.1, 0.4, 0.7, 1.0)
soil.thermalconductivity(examplesoil)
soil.thermalcapacity(examplesoil)
}
