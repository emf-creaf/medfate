\encoding{UTF-8}
\name{soil thermodynamics}
\alias{soil_thermalConductivity}
\alias{soil_thermalCapacity}
\alias{soil_temperatureChange}
\alias{soil_temperatureGradient}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil thermodynamic functions}
\description{
Functions \code{soil_thermalConductivity} and \code{soil_thermalCapacity} calculate thermal conductivity and thermal capacity for each soil layer, given its texture and water content. Functions \code{soil_temperatureGradient} and \code{soil_temperatureChange} are used to calculate soil temperature gradients (in ºC/m) and temporal temperature change (in ºC/s) given soil layer texture and water content (and possibly including heat flux from above).
}
\usage{
soil_thermalConductivity(soil, model = "SX")
soil_thermalCapacity(soil, model = "SX")
soil_temperatureChange(dVec, Temp, sand, clay, W, Theta_FC, Gdown)
soil_temperatureGradient(dVec, Temp)
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
Function \code{soil_thermalConductivity} returns a vector with values of thermal conductivity (W/m/ºK) for each soil layer. Function \code{soil_thermalCapacity} returns a vector with values of heat storage capacity (J/m3/ºK) for each soil layer. Function \code{soil_temperatureGradient} returns a vector with values of temperature gradient between consecutive soil layers. Function \code{soil_temperatureChange} returns a vector with values of instantaneous temperature change (ºC/s) for each soil layer.
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
soil_thermalConductivity(examplesoil)
soil_thermalCapacity(examplesoil)

#Values change when altering water content (drier layers have lower conductivity and capacity)
examplesoil$W = c(0.1, 0.4, 0.7, 1.0)
soil_thermalConductivity(examplesoil)
soil_thermalCapacity(examplesoil)
}
