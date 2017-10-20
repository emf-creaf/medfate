\encoding{UTF-8}
\name{soildgridsBuilder}
\alias{soildgridsBuilder}

\title{SoilGrids soil description fetcher}

\description{\code{soilgridsBuilder} takes a vector of depths and returns a
list of soil characteristics ready to use with \code{\link{soil}} function.}

\usage{
soilgridsBuilder(lat = , long = , depths = c(300, 500, 1200))
}

\arguments{
  \item{lat}{Latitude in decimal format}
  \item{long}{Longitude in decimal format}
  \item{depths}{A numeric vector indicating the desired depths, in \emph{cm}}
}

\details{
This function connects with the SoilGrids REST API (https://rest.soilgrids.org)
to retrieve the soil physical and chemical characteristics for a site, selected
by its coordinates. Also, in case the depths are not the default ones in the
SoilGrids API, the function uses the trapezoidal rule to calculate the values
for the desired depths (as described in Hengl \emph{et al}. 2007).  
To do this, it uses internally the \code{rjson} and \code{GSIF} R packages.
}

\value{
A list, containing the soil characteristics ready to use with the
\code{\link{soil}} function.
}

\author{
\enc{Víctor}{Victor} Granda, Centre Tecnologic Forestal de Catalunya
}

\references{
Hengl T, Mendes de Jesus J, Heuvelink GBM, Ruiperez Gonzalez M, Kilibarda M, Blagotić
A, et al. (2017) SoilGrids250m: Global gridded soil information based on machine
learning. PLoS ONE 12(2): e0169748. doi:10.1371/journal.pone.0169748
}

\examples{
  \dontrun{
  foo <- soilgridsBuilder(lat = 42.6667, long = -5.6333, depths = c(300, 600, 1200))
  foo_soil <- soil(foo)
  foo_soil
  }
}