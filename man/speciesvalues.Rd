\encoding{UTF-8}
\name{Species values}
\alias{species.BasalArea}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Species description functions}
\description{
Functions to calculate attributes of a \code{\link{forest}} object by species.
}
\usage{
species.BasalArea(x, SpParams)
}
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
}
\value{
A vector with values for each species in \code{SpParams}:
\itemize{
\item{\code{species.BasalArea}: Species basal area (m2/ha).}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swb}}, \code{\link{forest}}, \code{\link{plant.BasalArea}}, \code{\link{summary.forest}}
}
\examples{
#Default species parameterization
data(SpParamsMED)

#Load example plot
data(exampleforest)

#Species basal area
species.BasalArea(exampleforest, SpParamsMED)
}