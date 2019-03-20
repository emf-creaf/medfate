\encoding{UTF-8}
\name{spwbpoints}
\alias{spwbpoints}
\alias{growthpoints}

\title{Model simulations for spatially-distributed forest stands}
\description{
Functions \code{spwbpoints} and \code{growthpoints} allow calling local models \code{\link{spwb}} and \code{\link{growth}}, respectively, for a set of forest stands distributed in specific locations. No spatial processes are simulated.
}
\usage{
spwbpoints(y, SpParams, meteo, control = defaultControl(),
       dates = NULL, summaryFunction = NULL, args=NULL)
growthpoints(y, SpParams, meteo, control = defaultControl(),
       dates = NULL, summaryFunction = NULL, args=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An object of class \code{\link{SpatialPointsLandscape-class}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{meteo}{Meteorology data (see details).}
  \item{control}{A list of control parameters (see \code{\link{defaultControl}}).}
  \item{dates}{A \code{\link{Date}} object with the days of the period to be modeled. If \code{NULL}, then the whole period of \code{meteo} is used.}
  \item{summaryFunction}{An appropriate function to calculate summaries (e.g., \code{\link{summary.spwb}}).}
  \item{args}{List with additional arguments for the summary function.}
}
\details{
Functions \code{spwbpoints} and \code{growthpoints}  accept different formats for meteorological input (parameter \code{meteo}). If a \code{data.frame} is supplied (as in \code{\link{spwb}} or \code{\link{growth}}) then the same meteorology is used for all points (not recommended). To specify different meteorology for different points, the user can use an object of \code{\link{SpatialPointsMeteorology-class}}. Alternatively, the user can supply an object of class \code{\link{SpatialPointsDataFrame-class}} containing the meta data (columns dir and filename) of meteorological files that will be read from the disk.
}
\value{
Functions \code{spwbpoints} and \code{growthpoints} return a list with the following elements:
\itemize{
\item{\code{sp}: An object of class \code{\link{SpatialPoints}} with the spatial coordinates of forest stands.}
\item{\code{input}: A list of objects of class \code{\link{spwbInput}} or \code{\link{growthInput}} (one per forest stand). In the case of \code{growthpoints}, the input variables will have been modified by the call to the \code{\link{growth}} function.}
\item{\code{result}: A list the result of calling \code{\link{spwb}} or \code{\link{growth}} on each forest stand. If \code{summaryFunction} is not null, then each element of the list will contain the result of the summary function.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{spwb}}, \code{\link{growth}}, \code{\link{SpatialPointsLandscape-class}}
}
\examples{
\dontrun{
# Load spatial points forest object
data("exampleSPL")

#Load meteorological data (the same will be used for all points)
data("examplemeteo")

#Load species parameters
data("SpParamsMED")

# Deactivate extra console output
control = defaultControl()
control$verbose = FALSE

# Run simulation for the first two forest plots
# (you can run the code without subsetting 'exampleSPL')
res <- spwbpoints(exampleSPL, SpParamsMED, examplemeteo, control = control)

#Extract summaries for a given forest plot
summary(res$result$`80001`, freq="months",FUN=sum, 
        output="PlantStress", bySpecies = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.