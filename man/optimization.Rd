\encoding{UTF-8}
\name{optimization}
\alias{optimization}
\alias{optimization_evaluation_function}
\alias{optimization_prediction_function}
\title{
Function factories for optimization routines
}
\description{
Function factories to generate functions to be used in model calibration, uncertainty or sensitivity analysis.
}
\usage{
optimization_evaluation_function(x, soil, parNames, 
                                 measuredData, type = "SWC", 
                                 cohort = NULL, SpParams = NULL,
                                 metric = "loglikelihood",
                                 ...)
optimization_prediction_function(x, soil, parNames, 
                                 summary_function,
                                 ...)
}
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{parNames}{A string vector of parameter names (see examples and naming rules in \code{\link{modifyParams}}).}
  \item{measuredData}{A data frame with observed/measured values. Dates should be in row names, whereas columns should be named according to the type of output to be evaluated (see details).}
  \item{type}{A string with the kind of model output to be evaluated. Accepted values are \code{"SWC"} (soil moisture content), \code{"REW"} relative extractable water, \code{"ETR"} (total evapotranspiration), \code{"E"} (transpiration per leaf area) and \code{"WP"} (plant water potentials).}
  \item{cohort}{A string of the cohort to be compared (e.g. "T1_68"). If \code{NULL} results for the first cohort will be evaluated.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}), only needed if \code{type = "FMC"}.}
  \item{metric}{An evaluation metric (see \code{\link{evaluation_metric}}).}
  \item{summary_function}{A function whose input is the result of \code{\link{spwb}} or \code{\link{growth}} and returns a scalar.}
  \item{...}{Other parameters to be passed to functions \code{\link{spwb}} or \code{\link{growth}}.}
}
\details{
See \code{\link{evaluation}} for details regarding how to specify measured data.

Functions produced by these function factories should be useful for sensitivity analyses using package 'sensitivity'.

Parameter naming should follow the rules specified in section details of \code{\link{modifyParams}}.
}
\value{
Function \code{optimization_evaluation_function} returns a function whose parameters are parameter values and whose return is an evaluation metric (e.g. loglikelihood of the data observations given model predictions).

Function \code{optimization_evaluation_function} returns a function whose parameters are parameter values and whose return is a prediction scalar (e.g. total transpiration).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{evaluation_metric}}, \code{\link{modifyParams}}, \code{\link{spwb}}, \code{\link{growth}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (4 layers)
examplesoil1 = soil(defaultSoilParams(4))

#Initialize control parameters
control = defaultControl()

#Initialize input
x1 = forest2spwbInput(exampleforestMED,examplesoil1, SpParamsMED, control)

#Load observed data (in this case the same simulation results with some added error)  
data(exampleobs)

# Generate a loglikelihood function for soil water content
# as a function of parameters "Z50" and "Z95" for cohort "T1_54"
f<-optimization_evaluation_function(x1, examplesoil1,
                                    parNames = c("T1_54/Z50", "T1_54/Z95"), 
                                    exampleobs, type = "SWC", 
                                    metric = "loglikelihood",
                                    meteo = examplemeteo, latitude = 41.82592, elevation = 100)

#Loglikelihood value for Z50 = 200, Z95 = 1000
f(c(200, 1000))

#Define a summary function as the total transpiration over the simulated period
sf<-function(x) {sum(x$WaterBalance$Transpiration, na.rm=TRUE)}

# Generate a prediction function for total transpiration over the simulated period
# as a function of parameters "Z50" and "Z95" for cohort "T1_54"
f<-optimization_prediction_function(x1, examplesoil1,
                                    parNames = c("T1_54/Z50", "T1_54/Z95"), 
                                    summary_function = sf,
                                    meteo = examplemeteo, latitude = 41.82592, elevation = 100)

# Total transpiration value for Z50 = 200, Z95 = 1000
f(c(200, 1000))
}
