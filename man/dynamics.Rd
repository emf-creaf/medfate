\encoding{UTF-8}
\name{dynamics}
\alias{dynamics}
\title{Forest dynamics}
\description{
Function \code{dynamics} implements a forest dynamics model that simulates growth, mortality and recruitment for plant cohorts in a given forest stand during a period specified in the input climatic data. 
}
\usage{
dynamics(forest, soil, SpParams, meteo, control,
         latitude , elevation = NA, slope = NA, aspect = NA) 
}
\arguments{
  \item{forest}{An object of class \code{\link{forest}}.}
  \item{soil}{An object of class \code{\link{soil}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsMED}}).}
  \item{meteo}{A data frame with daily meteorological data series. Row names of the data frame should correspond to date strings with format "yyyy-mm-dd" (see \code{\link{Date}}).}
  \item{control}{A list with default control parameters (see \code{\link{defaultControl}}).}
  \item{latitude}{Latitude (in degrees). Required when \code{x$TranspirationMode = "Sperry"}.}
  \item{elevation, slope, aspect}{Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). Required when \code{x$TranspirationMode = "Sperry"}. Elevation is also required for 'Granier' if snowpack dynamics are simulated.}
}
\details{
  Function \code{dynamics} simulates forest dynamics for annual time steps, building on other simulation functions. For each simulated year, the function first calls function \code{\link{growth}} to simulate daily water/carbon balance, growth and mortality processes. Next, \code{dynamics} simulates recruitment for species already present in the stand and prepares the input of function \code{\link{growth}} for the next annual time step.
}
\value{
  A list of class 'dynamics' with the following elements:
  \itemize{
 \item{\code{"StandSummary"}: A data frame with stand-level summaries (leaf area index, tree basal area, tree density, shrub cover, etc.) at the beginning of the simulation and after each simulated year.}
 \item{\code{"SpeciesSummary"}: A data frame with species-level summaries (leaf area index, tree basal area, tree density, shrub cover, etc.) at the beginning of the simulation and after each simulated year.}
 \item{\code{"CohortSummary"}: A data frame with cohort-level summaries (leaf area index, tree basal area, tree density, shrub cover, etc.) at the beginning of the simulation and after each simulated year.}
 \item{\code{"TreeTable"}: A data frame with tree-cohort data (species, density, diameter, height, etc.) at the beginning of the simulation and after each simulated year.}
 \item{\code{"ShrubTable"}: A data frame with shrub-cohort data (species, density, cover, height, etc.) at the beginning of the simulation and after each simulated year.}
 \item{\code{"ForestStructures"}: A list with the structure of the stand at the beginning of the simulation and after each simulated year.}
 \item{\code{"GrowthResults"}: A list with the results of calling function \code{\link{growth}} for each simulated year.}
  }
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{growth}}, \code{\link{plot.growth}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)
#Prepare a three-year meteorological data
meteo2001 = examplemeteo
meteo2002 = examplemeteo
meteo2003 = examplemeteo
row.names(meteo2002) = seq(as.Date("2002-01-01"), 
                           as.Date("2002-12-31"), by="day")
row.names(meteo2003) = seq(as.Date("2003-01-01"), 
                           as.Date("2003-12-31"), by="day")
meteo_01_03 = rbind(meteo2001, meteo2002, meteo2003)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control = defaultControl()

#Initialize soil with default soil params (4 layers)
examplesoil = soil(defaultSoilParams(4))

#Call simulation function
fd<-dynamics(exampleforestMED, examplesoil, 
             SpParamsMED, meteo_01_03, control,
             latitude = 41.82592, elevation = 100)

#Stand-level summaries
fd$StandSummary

#Tree table by annual steps
fd$TreeTable
}
