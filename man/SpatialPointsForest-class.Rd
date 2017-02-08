\encoding{UTF-8}
\name{SpatialPointsForest-class}
\Rdversion{1.1}
\docType{class}
\alias{SpatialPointsForest-class}
\alias{[,SpatialPointsForest,ANY,ANY,ANY-method}

\title{Class \code{"SpatialPointsForest"}}
\description{
An S4 class that represents a set of forest stands along with their spatial location
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SpatialPointsForest", ...)} or by calls to function \code{\link{SpatialPointsForest}}.
}
\section{Slots}{
  \describe{
    \item{\code{forestlist}:}{Object of class \code{"list"} containing a set of \code{\link{forest}} objects.  }
    \item{\code{soillist}:}{Object of class \code{"list"} containing a set of \code{\link{soil}} objects.}
    \item{\code{data}:}{Object of class "data.frame" containing the elevation (in m), slope (in degrees) and aspect (in degrees from North) of every cell.}
    \item{\code{coords}:}{Object of class \code{"matrix"} with spatial coordinates. }
    \item{\code{bbox}:}{Object of class \code{"matrix"} with the boundary box. }
    \item{\code{proj4string}:}{Object of class \code{"CRS"} with the projection string. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SpatialPointsTopography}"}, directly.
Class \code{"\linkS4class{SpatialPointsDataFrame}"}, distance 2.
Class \code{"\linkS4class{SpatialPoints}"}, distance 3.
Class \code{"\linkS4class{Spatial}"}, by class "SpatialPoints", distance 4.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SpatialPointsTopography", i = "ANY", j = "ANY", drop = "ANY")}: subsets the points, the corresponding topography data, and the forest and soil lists; only rows (points can be subsetted.}
    \item{spatialForestSummary}{\code{signature(object = "SpatialPointsForest")}: calculates a summary function for all forest stands and returns an object of class \code{\link{SpatialPointsDataFrame-class}}. }
    \item{spatialSoilSummary}{\code{signature(object = "SpatialPointsForest")}: calculates a summary function for the soil of all stands and returns an object of class \code{\link{SpatialPointsDataFrame-class}}. }
 }
}

\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpatialPointsTopography}}, \code{\link{SpatialPointsForest}}, \code{\link{SpatialGridForest}}, \code{\link{spatialForestSummary}}, \code{\link{forest}}
}
\examples{
showClass("SpatialPointsForest")
}
\keyword{classes}
