\encoding{UTF-8}
\name{SpatialPointsLandscape}
\alias{SpatialPointsLandscape}
\alias{SpatialPixelsLandscape}
\alias{SpatialGridLandscape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates Spatial Landscape Classes
}
\description{
Functions to initialize spatial landscape classes.
}
\usage{
SpatialPointsLandscape(spt, forestlist, soillist)
SpatialPixelsLandscape(spxt, lct, forestlist, soillist, verbose=TRUE)
SpatialGridLandscape(sgt, lct, forestlist, soillist, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spt}{An object of class \code{\link{SpatialPointsTopography}}.}
  \item{spxt}{An object of class \code{\link{SpatialPixelsTopography}}.}
  \item{sgt}{An object of class \code{\link{SpatialGridTopography}}.}
  \item{lct}{A character vector with the land cover type of each grid cell (values should be 'wildland', 'agriculture', 'rock' or 'static').}
  \item{forestlist}{A list of objects of class 'forest' with the same number of elements as spatial points/pixels.}
  \item{soillist}{A list of objects of class 'forest' with the same number of elements as spatial points/pixels.}
  \item{verbose}{A flag to specify console output during spatial object initialization.}
}
\value{
An object of class \code{\link{SpatialPointsLandscape-class}}, \code{\link{SpatialPixelsLandscape-class}} or \code{\link{SpatialGridLandscape-class}} depending on the function.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}

\seealso{\code{\link{forest}}, \code{\link{soil}}}
