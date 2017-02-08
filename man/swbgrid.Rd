\encoding{UTF-8}
\name{swbgrid}
\alias{swbgrid}
\alias{plot.swbgrid}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil water balance and lateral water discharge}
\description{
Function \code{swbgrid} conducts daily soil water balance over a grid of cells while incorporating water runoff from upperslope cells into the current cell.
}
\usage{
swbgrid(y, SpParams, meteo, dates, 
           summaryFreq = "years", trackSpecies = numeric(), 
           control = defaultControl())
\method{plot}{swbgrid}(x, type = "Runon", summaryIndex = 1, spIndex = NULL, ...)
}
\arguments{
  \item{y}{An object of class \code{\link{SpatialGridForest-class}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{meteo}{A \code{\link{SpatialGridMeteorology-class}} object or a data frame with two columns: 'dir' and 'filename, to indicate the path to the meteorological data.}
  \item{dates}{A \code{\link{Date}} object describing the days of the period to be modeled.}
  \item{summaryFreq}{Frequency in which summary layers will be produced (e.g. "years", "months", ...) (see \code{\link{cut.Date}}).}
  \item{trackSpecies}{An integer vector containing the indices of species for which transpiration and drought stress is to be tracked.}
  \item{control}{A list of control parameters (see \code{\link{defaultControl}}).}
  \item{x}{An object of class \code{swbgrid}.}
  \item{type}{Type of information to be drawn.}
  \item{summaryIndex}{The index of the summary to be ploted.}
  \item{spIndex}{The index of the species to be ploted (for some types).}
  \item{...}{Additional parameters to function \code{\link{spplot}}.}
}
\details{
Function \code{swbgrid} requires daily meteorological data over a grid. The user may supply an object of class \code{\link{SpatialGridMeteorology}} (see package 'meteoland') or a data frame with information regarding where to read meteorological data.
}
\value{
A list of class 'swbgrid' with the following elements:
\itemize{
 \item{\code{grid}: The \code{\link{GridTopology}} object corresponding to the simulated area.}
 \item{\code{LandscapeBalance}: A data frame with as many rows as summary points and where columns are components of the water balance at the landscape level (i.e., precipitation, interception, soil evaporation, plant transpiration, ...).}
}
Then, the following matrices are included (each with as many rows as cells and as many columns as summary points):
\itemize{
   \item{\code{NetPrec}: Net precipitation, after accounting for interception (in mm).}
   \item{\code{Infiltration}: The amount of water infiltrating into the soil (in mm).}
   \item{\code{Runon}: The amount of water imported from other cells via surface runoff (in mm).}
   \item{\code{Runoff}: The amount of water exported via surface runoff (in mm).}
   \item{\code{DeepDrainage}: The amount of water exported via deep drainage (in mm).}
   \item{\code{Esoil}: Bare soil evaporation (in mm).}
   \item{\code{Eplant}: Plant transpiration (in mm).}
}
The same list contains two three-dimensional arrays (each with dimensions number of cells, number of summary layers and number of tracked species):
\itemize{
   \item{\code{Transpiration}: Total transpiration (in mm) of the tracked species for each summary period.}
   \item{\code{DI}: Drought intensity (from 0 to 1) of the tracked species for each summary period.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya.
}
\seealso{
\code{\link{swb}}, \code{\link{swbpoints}}, \code{\link{SpatialGridForest-class}}
}