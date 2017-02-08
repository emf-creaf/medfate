\encoding{UTF-8}
\name{exampleSPF}
\alias{exampleSPF}
\docType{data}
\title{
Example of spatial points with forest data
}
\description{
An example of an object of \code{\link{SpatialPointsForest-class}} with data for 30 plots, taken from the Spanish Forest Inventory (DGCN 2005).
}
\usage{data("exampleSPF")}
\format{
  The data format is that of an object \code{\link{SpatialPointsForest-class}} 
}
\source{
DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya. Dirección General de Conservación de la Naturaleza, Ministerio de Medio Ambiente, Madrid.
}
\seealso{ \code{\link{forest}}, \code{\link{exampleforest}}, \code{\link{SpatialPointsForest-class}}}
\examples{
data(exampleSPF)

#Plot forest coordinates
plot(exampleSPF)

#Inspect forest object corresponding to the first point
exampleSPF@forestlist[[1]]

}
\keyword{datasets}
