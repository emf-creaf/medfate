\encoding{UTF-8}
\name{spwbgrid}
\alias{spwbgrid}
\alias{plot.spwbgrid}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil-plant water balance and lateral water discharge}
\description{
Function \code{spwbgrid} conducts daily soil and plant water balance over a grid of cells while incorporating water runoff from upperslope cells into the current cell.
}
\usage{
spwbgrid(y, SpParams, meteo, dates, 
           summaryFreq = "years", trackSpecies = numeric(), 
           control = defaultControl())
\method{plot}{spwbgrid}(x, type = "Runon", summaryIndex = 1, spIndex = NULL, ...)
}
\arguments{
  \item{y}{An object of class \code{\link{SpatialGridLandscape-class}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{meteo}{A \code{\link{SpatialGridMeteorology-class}} object, a data frame with two columns: 'dir' and 'filename, to indicate the path to the meteorological data, or a data frame with meteorological data (the same for all cells).}
  \item{dates}{A \code{\link{Date}} object describing the days of the period to be modeled.}
  \item{summaryFreq}{Frequency in which summary layers will be produced (e.g. "years", "months", ...) (see \code{\link{cut.Date}}).}
  \item{trackSpecies}{An integer vector containing the indices of species for which transpiration and drought stress is to be tracked.}
  \item{control}{A list of control parameters (see \code{\link{defaultControl}}).}
  \item{x}{An object of class \code{spwbgrid}.}
  \item{type}{Type of information to be drawn.}
  \item{summaryIndex}{The index of the summary to be ploted.}
  \item{spIndex}{The index of the species to be ploted (for some types).}
  \item{...}{Additional parameters to function \code{\link{spplot}}.}
}
\details{
Function \code{spwbgrid} requires daily meteorological data over a grid. The user may supply: 
\enumerate{
\item{an object of class \code{\link{SpatialGridMeteorology}}.}
\item{a data frame with information regarding where to read meteorological data.}
\item{a data frame with meteorological data common for all cells of the grid.}
}
}
\value{
A list of class 'spwbgrid' with the following elements:
\itemize{
 \item{\code{grid}: The \code{\link{GridTopology}} object corresponding to the simulated area.}
 \item{\code{LandscapeBalance}: A data frame with as many rows as summary points and where columns are components of the water balance at the landscape level (i.e., rain, snow, interception, infiltration, soil evaporation, plant transpiration, ...).}
}
Then, the following matrices are included (each with as many rows as cells and as many columns as summary points):
\itemize{
   \item{\code{Rain}: Rainfall (in mm).}
   \item{\code{Snow}: Snowfall (in mm).}
   \item{\code{Interception}: Rainfall interception (in mm).}
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
\code{\link{spwb}}, \code{\link{spwbpoints}}, \code{\link{SpatialGridLandscape-class}}
}
\examples{
\dontrun{
# Load spatial grid landscape object
data("exampleSGL")

#Load meteorological data (the same will be used for all points)
data("examplemeteo")

#Load species parameters
data("SpParamsMED")

# Deactivate extra console output
control = defaultControl()
control$verbose = FALSE

# Run water balance simulation over the grid for five days
# (you can test longer periods)
res <- spwbgrid(exampleSGL, SpParamsMED, examplemeteo[1:5,], 
                control = control, summaryFreq = "months")
}
}