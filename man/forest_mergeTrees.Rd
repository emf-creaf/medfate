\encoding{UTF-8}
\name{Forest manipulation}
\alias{forest_mergeTrees}
\alias{forest_mergeShrubs}

\title{Forest utility functions}
\description{
Functions to manipulate a \code{\link{forest}} object.
}
\usage{
forest_mergeTrees(x, byDBHclass = TRUE)
forest_mergeShrubs(x, byHeightclass = TRUE)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{byDBHclass}{Boolean flag to indicate that 5-cm tree DBH classes should be kept separated.}
  \item{byHeightclass}{Boolean flag to indicate that 10-cm shrub height classes should be kept separated.}
}
\value{
Another \code{\link{forest}} object with merged trees or shrubs, depending on the function.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{spwb}}, \code{\link{forest}},  \code{\link{fordyn}}, \code{\link{summary.forest}}
}
