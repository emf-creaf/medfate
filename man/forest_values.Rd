\encoding{UTF-8}
\name{Forest values}
\alias{forest_basalArea}
\alias{forest_mergeTrees}

\title{Forest utility functions}
\description{
Functions to manipulate or calculate overall attributes of a \code{\link{forest}} object.
}
\usage{
forest_basalArea(x)
forest_mergeTrees(x)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
}
\value{
Depending on the function:
\itemize{
\item{\code{forest_basalArea}: Total basal area (m2/ha).}
\item{\code{forest_mergeTrees}: Another \code{\link{forest}} object with merged trees.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{spwb}}, \code{\link{forest}}, \code{\link{plant_basalArea}}, \code{\link{summary.forest}}
}
\examples{
#Default species parameterization
data(SpParamsMED)

#Load example plot
data(exampleforest)

#A short way to obtain total basal area
forest_basalArea(exampleforest)
}