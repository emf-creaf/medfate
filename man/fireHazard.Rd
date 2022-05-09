\encoding{UTF-8}
\name{fireHazard}
\alias{fireHazard}
\title{
Fire hazard
}
\description{
Estimates potential fire behaviour at each daily step of a simulation
}
\usage{
fireHazard(x, SpParams, forest = NULL, standardConditions = FALSE, 
           freq="days", fun = "max")
}
\arguments{
  \item{x}{An object of class \code{\link{spwb}}, \code{\link{pwb}}, \code{\link{growth}} or \code{\link{fordyn}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsDefinition}} and \code{\link{SpParamsMED}}).}
  \item{forest}{An object of class \code{\link{forest}} (needed if \code{x} is of class \code{\link{spwb}}, \code{\link{pwb}} or \code{\link{growth}}).}
  \item{standardConditions}{A flag to indicate that standard fire weather conditions are to be used (instead of deriving fuel moisture and windspeed from \code{x}).}
  \item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
  \item{fun}{Summary function (by default, maximum values).}
}
\details{
Live fuel moisture of shrub and canopy layers is estimated from plant water status. Dead fuel moisture is estimated following Resco-de-Dios et al. (2015).
}
\value{
A matrix with fire behaviour variables (columns) for each simulated day (rows) or coarser time steps if summaries are requested.
}
\references{
Resco de Dios, V., A. W. Fellows, R. H. Nolan, M. M. Boer, R. A. Bradstock, F. Domingo, and M. L. Goulden. 2015. A semi-mechanistic model for predicting the moisture content of fine litter. Agricultural and Forest Meteorology 203:64–73.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{spwb}}, \code{\link{fuel_FCCS}}, \code{\link{fire_FCCS}}
}