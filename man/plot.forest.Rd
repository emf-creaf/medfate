\name{plot.forest}
\alias{plot.forest}
\alias{shinyplot.forest}
\title{
Plot forest attributes
}
\description{
Convenient wrappers for vertical forest profiles (see \code{\link{vprofile_leafAreaDensity}}).
}
\usage{
\method{plot}{forest}(x, SpParams, type = "LeafAreaDensity", 
                      byCohorts = FALSE, bySpecies = FALSE, ...)
\method{shinyplot}{forest}(x, SpParams, ...)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{type}{A string of the plot type: "LeafAreaDensity", "RootDistribution", "FuelBulkDensity", "PARExtinction", "SWRExtinction" or "WindExtinction". }
  \item{byCohorts}{Flag to separate profiles for each cohort.}
  \item{bySpecies}{Flag to aggregate results by species.}
  \item{\dots}{Additional parameters to vertical profiles}
}
\value{
A ggplot or a shiny application, depending on the function.
}

\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}

\seealso{
\code{\link{vprofile_leafAreaDensity}}
}
\examples{
data(exampleforestMED)
data(SpParamsMED)
plot(exampleforestMED, SpParamsMED)
}