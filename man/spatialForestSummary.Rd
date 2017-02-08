\encoding{UTF-8}
\name{spatialForestSummary}
\alias{spatialForestSummary}
\alias{spatialSoilSummary}
\alias{spatialForestSummary-methods}
\alias{spatialSoilSummary-methods}
\alias{spatialForestSummary,SpatialPointsForest-method}
\alias{spatialForestSummary,SpatialGridForest-method}
\alias{spatialSoilSummary,SpatialPointsForest-method}
\alias{spatialSoilSummary,SpatialGridForest-method}
\title{Forest and soil summaries over space}
\description{
 Functions to calculates a summary function for the forest or soil of all cells in a \code{\link{SpatialPointsForest-class}} or \code{\link{SpatialGridForest-class}}.
}
\usage{
spatialForestSummary(object, summaryFunction, ...)
spatialSoilSummary(object, summaryFunction, ...)
}
\arguments{
  \item{object}{An object of class \code{\link{SpatialPointsForest-class}} or \code{\link{SpatialGridForest-class}}.}
  \item{summaryFunction}{A function that accepts objects of class \code{\link{forest}} or \code{soil}, respectively.}
  \item{...}{Additional arguments to the summary function.}
}
\value{
An object of class \code{\link{SpatialPointsDataFrame}} or \code{\link{SpatialGridDataFrame}}, depending on the input, containing the calculated statistics.
}

\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya.
}
\seealso{
\code{\link{forest}}, \code{\link{soil}}, \code{\link{summary.forest}}
}
\examples{
#Load plot data
data(exampleSPF)

#Load species parameters
data(SpParamsMED)

#Apply summary function
y <- spatialForestSummary(exampleSPF,summary.forest, SpParamsMED)
head(y@data)

#Plot basal area
spplot(y["BA"])

}