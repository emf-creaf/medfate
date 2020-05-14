\encoding{UTF-8}
\name{pheno_updateLeaves}
\alias{pheno_updateLeaves}
\alias{pheno_leafDevelopmentStatus}
\alias{pheno_leafSenescenceStatus}
\title{
Leaf phenology
}
\description{
Function \code{pheno_leafDevelopmentStatus} returns the expanded status (0 to 1) of leaves according to the growth degree days required to start bud burst and leaf unfolding, as dictated by a simple ecodormancy (one-phase) model (Chuine et al. 2013). Function \code{pheno_leafSenescenceStatus} returns the 0/1 senescence status of leaves according to the one-phase senescence model of Delpierre et al. (2009) on the basis of photoperiod and temperature.  
Function \code{pheno_updateLeaves} updates the status of expanded leaves and dead leaves of object \code{x} given the photoperiod, temperature and wind of a given day. It applies the development model for 1 < doy < 180 and the senescence model for 181 > doy > 365.
}
\usage{
pheno_leafDevelopmentStatus(Sgdd, gdd, unfoldingDD = 300)
pheno_leafSenescenceStatus(Ssen, sen)
pheno_updateLeaves(x, doy, photoperiod, tmean, wind)
}
\arguments{
  \item{Sgdd}{Degree days required for leaf budburst (in Celsius).}
  \item{gdd}{Cumulative degree days (in Celsius)}
  \item{unfoldingDD}{Degree-days for complete leaf unfolding after budburst has occurred.}
  \item{Ssen}{Threshold to start leaf senescence.}
  \item{sen}{Cumulative senescence variable.}
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{doy}{Day of the year.}
  \item{photoperiod}{Day length (in hours).}
  \item{tmean}{Average day temperature (in Celsius).}
  \item{wind}{Average day wind speed (in m/s).}
}
\value{
Function \code{pheno_leafDevelopmentStatus} returns a vector of values between 0 and 1, whereas function \code{pheno_leafSenescenceStatus} returns a vector of 0 (senescent) and 1 (expanded) values.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\note{
Function \code{pheno_updateLeaves} modifies the object \code{x}.
}
\references{
Chuine, I., De Cortazar-Atauri, I.G., Kramer, K., \enc{Hänninen}{Hanninen}, H., 2013. Plant development
models. Phenology: An Integrative Environmental Science. Springer, pp. 275–293.

Delpierre N, Dufrêne E, Soudani K et al (2009a) Modelling interannual and spatial variability of leaf senescence for three deciduous tree species in France. Agric For Meteorol 149:938–948. doi:10.1016/j.agrformet.2008.11.014
}
\seealso{
\code{\link{spwb}}
}
