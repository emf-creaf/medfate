\encoding{UTF-8}
\name{exampleSPL}
\alias{exampleSPL}
\docType{data}
\title{
Example of spatial points with forest data
}
\description{
An example of an object of \code{\link{SpatialPointsLandscape-class}} with data for 30 plots, taken from the Spanish Forest Inventory (DGCN 2005).
}
\usage{data("exampleSPL")}
\format{
  The data format is that of an object \code{\link{SpatialPointsLandscape-class}} 
}
\source{
DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya. Dirección General de Conservación de la Naturaleza, Ministerio de Medio Ambiente, Madrid.
}
\seealso{ \code{\link{forest}}, \code{\link{exampleforest}}, \code{\link{SpatialPointsLandscape-class}}}
\examples{
data(exampleSPL)

#Plot forest coordinates
plot(exampleSPL)

#Inspect forest object corresponding to the first point
exampleSPL@forestlist[[1]]

}
\keyword{datasets}
