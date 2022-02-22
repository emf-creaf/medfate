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
fireHazard(x, SpParams, forest = NULL, standardConditions = FALSE)
}
\arguments{
  \item{x}{An object of class \code{\link{spwb}}, \code{\link{pwb}}, \code{\link{growth}} or \code{\link{fordyn}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsDefinition}} and \code{\link{SpParamsMED}}).}
  \item{forest}{An object of class \code{\link{forest}} (needed if \code{x} is of class \code{\link{spwb}}, \code{\link{pwb}} or \code{\link{growth}}).}
  \item{standardConditions}{A flag to indicate that standard fire weather conditions are to be used (instead of deriving fuel moisture and windspeed from \code{x}).}
}
\value{
A matrix with fire behaviour variables (columns) for each simulated day (rows)
}
\references{
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{spwb}}, \code{\link{fire_FCCS}}
}