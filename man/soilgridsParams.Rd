\encoding{UTF-8}
\name{soilgridsParams}
\alias{soilgridsParams}

\title{SoilGrids soil description fetcher}

\description{\code{soilgridsParams} takes a vector of depths and returns a
list of soil characteristics ready to use with \code{\link{soil}} function.}

\usage{
soilgridsParams(lat, long, depths = c(300, 500, 1200))
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
To do this, it uses internally the \code{rjson} and \code{GSIF} R packages. Macroporosity is calculated from bulk density and sand percentage following Stolf et al. (2011).
}

\value{
A list, containing the soil characteristics ready to use with the
\code{\link{soil}} function.
}

\author{
\enc{Víctor}{Victor} Granda, Centre Tecnologic Forestal de Catalunya
}

\references{
Hengl T, Mendes de Jesus J, Heuvelink GBM, Ruiperez Gonzalez M, Kilibarda M, \enc{Blagotić}{Blagotic}
A, et al. (2017) SoilGrids250m: Global gridded soil information based on machine
learning. PLoS ONE 12(2): e0169748. doi:10.1371/journal.pone.0169748

Stolf, R., Thurler, Á., Oliveira, O., Bacchi, S., & Reichardt, K. 2011. Method to estimate soil macroporosity and microporosity based on sand content and bulk density. Revista Brasileira de Ciencias do Solo 35: 447–459.
}
\seealso{
 \code{\link{soil}}, \code{\link{defaultSoilParams}}
}
\examples{
  \dontrun{
  foo <- soilgridsParams(lat = 42.6667, long = -5.6333, depths = c(300, 600, 1200))
  foo_soil <- soil(foo)
  foo_soil
  }
}