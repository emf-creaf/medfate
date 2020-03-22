\encoding{UTF-8}
\name{soilgridsParams}
\alias{soilgridsParams}

\title{SoilGrids soil description fetcher}

\description{\code{soilgridsParams} takes a vector of depths and returns a
list of soil characteristics ready to use with \code{\link{soil}} function (requires package GSIF).}

\usage{
soilgridsParams(points, depths = c(300, 500, 1200), verbose = FALSE)
}

\arguments{
  \item{points}{An object of class (or subclass) \code{SpatialPoints} with a valid CRS definition. }
  \item{depths}{A numeric vector indicating the desired depths, in \emph{mm}}
  \item{verbose}{A logical flag to include a progress bar while processing the output of the query to the SoilGrids REST API.}
}

\details{
This function connects with the SoilGrids REST API (https://rest.soilgrids.org)
to retrieve the soil physical and chemical characteristics for a site, selected
by its coordinates. Also, in case the depths are not the default ones in the
SoilGrids API, the function uses the trapezoidal rule to calculate the values
for the desired depths (as described in Hengl \emph{et al}. 2007).  
To do this, it uses internally the \code{rjson} and \code{GSIF} R packages. Although soil layer definition is taken from the input, the function includes absolute depth to bedrock and depth to the R horizon from SoilGrids, in case the user wants to redefine soil definition. Input coordinates are transformed to longitude/latitude within the function.
}

\value{
If only one point is supplied, the function returns a list with the following elements:
\itemize{
  \item{\code{soilparams}: A data frame containing the soil characteristics ready to use with the
\code{\link{soil}} function.}
  \item{\code{soilgrids_Rhorizondepth}: Depth to the R horizon (in mm).}
  \item{\code{soilgrids_Rhorizonprob}: Probability of the R horizon.}
  \item{\code{soilgrids_absolutesoildepth}: Absolute depth to bedrock (in mm).}
}
If more than one point is supplied, the function returns a list with as many elements as points, each one containing the mentioned list.
}

\author{
\enc{Víctor}{Victor} Granda, CREAF/CTFC

Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
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
  library(GSIF)
  library(rjson)
  coords_sp <- SpatialPoints(cbind(long = -5.6333, lat = 42.6667), 
                             CRS("+proj=longlat +datum=WGS84"))
  foo <- soilgridsParams(coords_sp, depths = c(300, 600, 1200))
  foo_soil <- soil(foo$soilparams)
  foo_soil
}
}