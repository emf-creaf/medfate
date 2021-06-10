\encoding{UTF-8}
\name{modifyParams}
\alias{modifySpParams}
\alias{modifyCohortParams}
\alias{modifyInputParams}
\title{
Modify parameters
}
\description{
Routines to modify species parameter table or model input objects
}
\usage{
modifySpParams(SpParams, customParams, subsetSpecies = TRUE)
modifyCohortParams(x, customParams, verbose = TRUE)
modifyInputParams(x, customParams, verbose = TRUE)
}
\arguments{
  \item{x}{A model input object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.}
  \item{SpParams}{A species parameter data frame, typically \code{\link{SpParamsMED}}.}
  \item{customParams}{A data frame or a named vector with new parameter values (see details).}
  \item{subsetSpecies}{A flag to indicate that the output data frame should include only those species mentioned in \code{customParams}.}
  \item{verbose}{A flag to indicate that messages should be printed on the console.}
}
\details{
When calling function \code{modifySpParams}, \code{customParams} should be a data frame with as many rows as species and as many columns as parameters to modify, plus a column called 'SpIndex' to match species between the two tables.

When calling \code{modifyCohortParams}, \code{customParams} can be a data frame with as many rows as cohorts and as many columns as parameters to modify, plus a column called 'Cohort' which will be matched with the cohort names given by \code{\link{spwbInput}} or \code{\link{growthInput}}. Alternatively, \code{customParams} can be a named list or named numeric vector as for \code{modifyInputParams}.

When calling \code{modifyInputParams}, \code{customParams} must be either a named list or a named numeric vector. Cohort parameters are specified using the syntax "<cohortName>/<paramName>" for names (e.g. "T2_68/Z50" to modify parameter 'Z50' of cohort 'T2_68'). Soil layer parameters are specified using the syntax "<paramName>@#layer" for names, where #layer is the layer index (e.g. "rfc@1" will modify the rock fragment content of soil layer 1). Control parameters are specified using either "<paramName>" (e.g "phloemConductanceFactor") or "<paramName>$<subParamName>" (e.g "maximumRelativeGrowthRates$leaf"). It may seem unnecessary to modify soil or control parameters via a function, but \code{modifyInputParams} is called from optimization functions (see \code{\link{optimization}}).
}
\value{
Function \code{modifySpParams} returns a modified species parameter data frame. 

Functions \code{modifyCohortParams} and \code{modifyInputParams} return a modified \code{\link{spwbInput}} or \code{\link{growthInput}} object. Note that modifications may affect other parameters beyond those indicated in \code{customParams}, as a result of parameter dependencies (see examples).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{spwbInput}}, \code{\link{SpParamsMED}}, \code{\link{optimization}}
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
control = defaultControl("Granier")

#Initialize input
x1 = forest2spwbInput(exampleforestMED,examplesoil1, SpParamsMED, control)

# Modify Z50 and Z95 of cohort T1_54
x1m <- modifyInputParams(x1, c("T1_54/Z50" = 200, "T1_54/Z95" = 2000))

# Inspect original and modified objects 
x1$below
x1m$below

# Inspect dependencies: fine root distribution across soil layers
x1$belowLayers$V
x1m$belowLayers$V

# Modify rock fragment content and sand proportion of soil layer 1
x1s <- modifyInputParams(x1, c("rfc@1" = 5, "sand@1" = 10))

# Inspect original and modified soils 
x1$soil
x1s$soil

# When modifying growth input objects dependencies increase
x1 = forest2growthInput(exampleforestMED,examplesoil1, SpParamsMED, control)
customParams = c("T1_54/Al2As" = 2000, "T1_54/LAI_live" = 2)
x1m <- modifyInputParams(x1, customParams)


}
