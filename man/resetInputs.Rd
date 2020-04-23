\encoding{UTF-8}
\name{resetInputs}
\alias{resetInputs}
\title{Reset simulation inputs}
\description{
Function \code{resetInputs()} allows resetting state variables in \code{x} and \code{soil} to their defaults.
}
\usage{
resetInputs(x, soil)
}
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
}
\value{
Does not return any value. Instead, it modifies input objects \code{x} and \code{soil}.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{spwbInput}}, \code{\link{growthInput}}, \code{\link{spwb}}
}