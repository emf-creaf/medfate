\encoding{UTF-8}
\name{spwb_maximumTranspirationRatioPlot}
\alias{spwb_maximumTranspirationRatioPlot}
\title{
Maximum transpiration vs. LAI
}
\description{
Builds a curve of maximum transpiration (Tmax) over potential evapotranspiration (ETP) for increasing leaf area index (LAI) values.
}
\usage{
spwb_maximumTranspirationRatioPlot(x, soil, meteo, latitude, elevation, slope, aspect, 
                                   ndays = 100, 
                                   LAI_seq = c(0.1, 0.25, seq(0.5, 10, by = 0.5)))
}
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{meteo}{A data frame with daily meteorological data series.}
  \item{latitude}{Latitude (in degrees). Required when using the 'Sperry' transpiration mode.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when using the 'Sperry' transpiration mode.}  \item{ndays}{Number of days to subset for curve building.}
  \item{LAI_seq}{Sequence of stand LAI values to be tested.}
}
\details{
TO BE DONE
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{spwb}}
}
\examples{

\dontrun{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params
examplesoil2 = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl()

#Switch to 'Sperry' transpiration mode
control$transpirationMode="Sperry"

#Initialize input
x2 = forest2spwbInput(exampleforestMED,examplesoil2, SpParamsMED, control)

#Run optimization
spwb_maximumTranspirationRatioPlot(x2, examplesoil2, examplemeteo, 
                                   41.82592, elevation = 100, 
                                   slope = 0, aspect = 0, ndays = 20)
}
}
