\encoding{UTF-8}
\name{Forest values}
\alias{forest_basalArea}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Forest description functions}
\description{
Functions to calculate overall attributes of a \code{\link{forest}} object.
}
\usage{
forest_basalArea(x)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
}
\value{
A vector with values for each cohort of the input \code{\link{forest}} object:
\itemize{
\item{\code{forest_basalArea}: Total basal area (m2/ha).}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
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