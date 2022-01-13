\name{recruitment}
\alias{recruitment}
\title{
Plant recruitment
}
\description{
Annual plant recruitment observed in a forest stand
}
\usage{
recruitment(forest, SpParams, control, minMonthTemp, moistureIndex)
}
\arguments{
  \item{forest}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsDefinition}}).}
  \item{control}{A list with default control parameters (see \code{\link{defaultControl}}).}
  \item{minMonthTemp}{Minimum month temperature.}
  \item{moistureIndex}{Moisture index (annual precipitation over annual potential evapotranspiration).}
}
\details{
Species can recruit if adults (sufficiently tall individuals) are present (seed rain can also be specified in a control parameter). Minimum month temperature and moisture index values are used to determine if recruitment was successful. Species also require a minimum amount of light at the ground level.
}
\value{An object of class \code{\link{forest}} with the new plant cohorts.}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{\code{\link{fordyn}}}
\examples{
#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control = defaultControl("Granier")

#Recruitment limits
plant_parameter(exampleforestMED, SpParamsMED, "MinTempRecr")
plant_parameter(exampleforestMED, SpParamsMED, "MinMoistureRecr")

#Compare recruitment outcomes
recruitment(exampleforestMED, SpParamsMED, control, 0, 0.25)
recruitment(exampleforestMED, SpParamsMED, control, 3, 0.25)

}
