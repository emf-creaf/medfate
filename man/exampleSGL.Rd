\name{exampleSGL}
\alias{exampleSGL}
\encoding{UTF-8}
\docType{data}
\title{
Example of spatial grid with forest data
}
\description{
An example of an object of \code{\link{SpatialGridLandscape-class}}, with data taken from the Spanish Forest Inventory (DGCN 2005).
}
\usage{data("exampleSGL")}
\format{
  The data format is that of an object \code{\link{SpatialGridLandscape-class}} 
}

\source{
DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya. Dirección General de Conservación de la Naturaleza, Ministerio de Medio Ambiente, Madrid.
}
\seealso{ \code{\link{forest}}, \code{\link{exampleforest}}, \code{\link{SpatialGridLandscape-class}}}
\examples{
data(exampleSGL)

#Inspect forest object corresponding to the first pixel
exampleSGL@forestlist[[1]]

}
\keyword{datasets}