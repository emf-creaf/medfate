\encoding{UTF-8}
\name{resetInputs}
\alias{resetInputs}
\title{Reset simulation inputs}
\description{
Function \code{resetInputs()} allows resetting state variables in \code{x} and \code{soil} to their defaults, or to copy values of state variables from a previous simulation results stored in \code{from}.
}
\usage{
resetInputs(x, soil, from = NULL, day = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{from}{An object of class \code{spwb} or \code{growth} storing the results of a previous simulation, including values of state variables. If \code{from = NULL}, state variables are set to their defaults (i.e. soil moisture set to field capacity and cumulative growth degree days set to zero).}
  \item{day}{An integer with the day from which state variable values stored in \code{from} should be copied. If missing, values are copied from the first day of stored values.}
}
\value{
Does not return any value. Instead, it modifies input objects \code{x} and \code{soil}.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{spwbInput}}, \code{\link{growthInput}}
}