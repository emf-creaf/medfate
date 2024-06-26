% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/soil_redefineLayers.R
\name{soil_redefineLayers}
\alias{soil_redefineLayers}
\title{Redefine soil layer widths}
\usage{
soil_redefineLayers(x, widths = c(300, 700, 1000, 2000))
}
\arguments{
\item{x}{A data frame of soil parameters (see an example in \code{\link{defaultSoilParams}}) or an object of class \code{\link{soil}}.}

\item{widths}{A numeric vector indicating the desired layer widths, in mm.}
}
\value{
A data frame or \code{\link{soil}} object with soil parameters, depending on the class of \code{x}.
}
\description{
Allows redefining soil layer widths of an input data frame of soil parameters.
}
\details{
If an initialized \code{\link{soil}} is supplied, its hydraulic parameters will be recalculated and the value of state variables will be lost.
}
\examples{
# Define initial soil with 5 layers
spar <- defaultSoilParams(5)
spar

# Redefine to four layers
soil_redefineLayers(spar)

# Same but after soil parameter initialization
examplesoil <- soil(spar)
examplesoil

soil_redefineLayers(examplesoil)

}
\seealso{
\code{\link{soil}}, \code{\link{defaultSoilParams}}
}
\author{
\enc{Víctor}{Victor} Granda, EMF-CREAF

Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
}
