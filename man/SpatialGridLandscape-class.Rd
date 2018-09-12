\encoding{UTF-8}
\name{SpatialGridLandscape-class}
\Rdversion{1.1}
\docType{class}
\alias{SpatialGridLandscape-class}

\title{Class \code{"SpatialGridLandscape"}}
\description{
An S4 class that represents a landscape configuration over a grid of coordinates.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SpatialGridLandscape", ...)}, or by calls to the function \code{\link{SpatialGridLandscape}}.
}
\section{Slots}{
  \describe{
    \item{\code{grid}:}{Object of class \code{"GridTopology"}. }
    \item{\code{bbox}:}{Object of class \code{"matrix"} with the boundary box. }
    \item{\code{proj4string}:}{Object of class \code{"CRS"} with the projection string.}   
    \item{\code{data}:}{Object of class "data.frame" containing the elevation (in m), slope (in degrees) and aspect (in degrees from North) of every cell.}
    \item{\code{lct}:}{A character vector of land cover type for each cell. Values allowed are: 'wildland', 'agriculture', 'rock' and 'static'.  }
    \item{\code{forestlist}:}{Object of class \code{"list"} containing a set of \code{\link{forest}} objects.  }
    \item{\code{soillist}:}{Object of class \code{"list"} containing a set of \code{\link{soil}} objects.}
    \item{\code{waterOrder}:}{A numeric vector of cell processing order.}
    \item{\code{waterQ}:}{A list of water discharge values to neighbours.}
    \item{\code{queenNeigh}:}{A list of neighbours for each cell.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SpatialGridTopography}"}, directly.
Class \code{"\linkS4class{SpatialGridDataFrame}"}, distance 2.
Class \code{"\linkS4class{SpatialGrid}"}, distance 3.
Class \code{"\linkS4class{Spatial}"}, by class "SpatialGrid", distance 4.
}
\section{Methods}{
  \describe{
    \item{getLCTs}{\code{signature(object = "SpatialGridLandscape")}: returns a SpatialGridDataFrame with the land cover types of the landscape cells.}
    \item{spatialForestSummary}{\code{signature(object = "SpatialGridLandscape")}: calculates a summary function for all forest stands and returns an object of class \code{\link{SpatialGridDataFrame-class}}. }
    \item{spatialSoilSummary}{\code{signature(object = "SpatialGridLandscape")}: calculates a summary function for the soil of all stands and returns an object of class \code{\link{SpatialGridDataFrame-class}}. }
    \item{spplot}{\code{signature(object = "SpatialGridLandscape")}: allows plotting maps of the landscape state.}
	 }
}

\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}

\seealso{
\code{\link{SpatialGridTopography-class}}, \code{\link{SpatialPointsLandscape-class}}, \code{\link{spatialForestSummary}}, \code{\link{forest}}, \code{\link{soil}}
}
\examples{
#Structure of the S4 object
showClass("SpatialGridLandscape")
}
\keyword{classes}
