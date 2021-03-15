\encoding{UTF-8}
\name{transp_maximumTranspirationModel}
\alias{transp_maximumTranspirationModel}
\title{
Maximum transpiration vs. LAI
}
\description{
Builds a model of maximum transpiration (Tmax) over potential evapotranspiration (PET) for increasing leaf area index (LAI) values for each plant cohort.
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
This function performs a meta-modelling exercise using the Sperry transpiration model, with the aim to estimate coefficients for the equation used in the Granier transpiration model (Granier et al. 1999). The model to be fitted is: \code{y ~ a*LAI + b*LAI^2}, where \code{y} is the ratio between maximum transpiration (Tmax) and Penman's potential evapotranspiration (PET) and \code{LAI} is the stand LAI. Unlike the original equation of Granier et al. (1999), we fit a zero intercept model so that LAI = 0 translates into zero plant transpiration. 

The function fits the model for each cohort separately, assuming it represents the whole stand. For each stand LAI value in the input sequence, the function uses simulations with Sperry transpiration and the input weather to estimate \code{y = Tmax/PET} as a function of stand's LAI (deciduous stands include leaf phenology). Once simulations have been conducted for each stand LAI value, the function fits a Generalized Linear Model with the above equation, assuming a Gamma distribution of residuals and an identity link.

The coefficients of the model can be used to parametrize Granier's transpiration, since coefficients \code{a} and \code{b} in the equation above correspond to parameters \code{Tmax_LAI} and \code{Tmax_LAIsq}, respectively (see \code{\link{SpParamsMED}}).
}
\value{
Returns a list with as many elements as plant cohorts, each element being a \code{\link{glm}} model.
}
\references{
Granier A, \enc{Bréda}{Breda} N, Biron P, Villette S (1999) A lumped water balance model to evaluate duration and intensity of drought constraints in forest stands. Ecol Modell 116:269–283. https://doi.org/10.1016/S0304-3800(98)00205-1.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{spwb}}, \code{\link{transp_transpirationGranier}}, \code{\link{transp_transpirationSperry}}, \code{\link{SpParamsMED}}
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
