\encoding{UTF-8}
\name{pheno_updateLeaves}
\alias{pheno_updateLeaves}
\alias{pheno_leafDevelopmentStatus}
\title{
Leaf phenology
}
\description{
Function \code{pheno_leafDevelopmentStatus} returns the fenological status ([0-1]) of species according to the growth degree days required for leaf flush. Function \code{pheno_updateLeaves} updates the status of expanded leaves and dead leaves of object \code{x} given the temperature and wind of a given day.
}
\usage{
pheno_leafDevelopmentStatus(Sgdd, gdd)
pheno_updateLeaves(x, doy, tmean, wind, Tbase = 5)
}
\arguments{
  \item{Sgdd}{A numeric vector of growth degree days for leaf budburst (in Celsius).}
  \item{gdd}{Growth degree days (in Celsius)}
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{doy}{Day of the year.}
  \item{tmean}{Average day temperature (in Celsius).}
  \item{wind}{Average day wind speed (in m/s).}
  \item{Tbase}{Base temperature to update growth degree days (in Celsius).}
}
\value{
Function \code{pheno_leafDevelopmentStatus} returns a vector of values between 0 and 1.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\note{
Function \code{pheno_updateLeaves} modifies the object \code{x}.
}
\seealso{
\code{\link{spwb}}
}
