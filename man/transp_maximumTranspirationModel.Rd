\encoding{UTF-8}
\name{transp_maximumTranspirationModel}
\alias{transp_maximumTranspirationModel}
\title{
Maximum transpiration vs. LAI
}
\description{
Builds a model of maximum transpiration (Tmax) over potential evapotranspiration (ETP) for increasing leaf area index (LAI) values for each plant cohort.
}
\usage{
transp_maximumTranspirationModel(x, meteo, latitude, elevation, slope, aspect, 
                                 LAI_seq = c(0.1, 0.25, seq(0.5, 10, by = 0.5)),
                                 draw = TRUE)
}
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}, built using the 'Sperry' transpiration mode.}
  \item{meteo}{A data frame with daily meteorological data series.}
  \item{latitude}{Latitude (in degrees).}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). }   
  \item{LAI_seq}{Sequence of stand LAI values to be tested.}
  \item{draw}{Boolean flag to indicate plotting of results.}
}
\details{
This function performs a meta-modelling exercise using the Sperry transpiration model, with the aim to estimate coefficients for the equation used in the Granier transpiration model.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{spwb}}
}
\examples{

\dontrun{
#Load example daily meteorological data
data(examplemeteo)

# Load example plot plant data
data(exampleforestMED)

# Load default species parameters
data(SpParamsMED)

# Initialize soil with default soil params
examplesoil2 = soil(defaultSoilParams(2))

# Initialize control parameters for 'Sperry' transpiration mode
control = defaultControl(transpirationMode="Sperry")

# Initialize input
x2 = forest2spwbInput(exampleforestMED,examplesoil2, SpParamsMED, control)

# Estimate maximum transpiration ratio models for each cohort
m = transp_maximumTranspirationModel(x2, examplemeteo, 
                                     41.82592, elevation = 100, 
                                     slope = 0, aspect = 0)

# Inspect the model for first cohort
m[[1]]
}
}
