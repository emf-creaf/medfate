\encoding{UTF-8}
\name{soilgridsParams}
\alias{soilgridsParams}

\title{SoilGrids soil description fetcher}

\description{\code{soilgridsParams} takes a vector of depths and returns a
list of soil characteristics ready to use with \code{\link{soil}} function.}

\usage{
soilgridsParams(points, widths = c(300, 700, 1000, 2000), verbose = FALSE)
}

\arguments{
  \item{points}{An object of class (or subclass) \code{SpatialPoints} with a valid CRS definition. }
  \item{widths}{A numeric vector indicating the desired layer widths, in \emph{mm}. If \code{NULL} the default soil grids layer definition is returned.}
  \item{verbose}{A logical flag to include a progress bar while processing the output of the query to the SoilGrids REST API.}
}

\details{
This function connects with the SoilGrids REST API (https://rest.soilgrids.org)
to retrieve the soil physical and chemical characteristics for a site (Hengl \emph{et al}. 2007), selected
by its coordinates. Also, in case the depths are not the default ones in the SoilGrids API, the function uses averages the values of soil grid layers depending on the overlap between soil layer definitions.  Input coordinates are transformed to longitude/latitude within the function.
}

\value{
If only one point is supplied, a data frame containing the soil characteristics ready to use with the
\code{\link{soil}} function. If more than one point is supplied, the function returns a list with as many elements as points, each one containing the mentioned list.
}

\author{
\enc{Víctor}{Victor} Granda, CREAF/CREAF

Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}

\references{
Hengl T, Mendes de Jesus J, Heuvelink GBM, Ruiperez Gonzalez M, Kilibarda M, \enc{Blagotić}{Blagotic}
A, et al. (2017) SoilGrids250m: Global gridded soil information based on machine
learning. PLoS ONE 12(2): e0169748. doi:10.1371/journal.pone.0169748.
}
\seealso{
 \code{\link{soil}}, \code{\link{defaultSoilParams}}
}
\examples{
\dontrun{
  coords_sp <- SpatialPoints(cbind(long = -5.6333, lat = 42.6667), 
                             CRS("+proj=longlat +datum=WGS84"))
  foo <- soilgridsParams(coords_sp, widths = c(300, 700, 1000))
  foo_soil <- soil(foo)
  foo_soil
}
}
