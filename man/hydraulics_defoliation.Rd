% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{hydraulics_defoliation}
\alias{hydraulics_defoliation}
\alias{hydraulics_proportionDefoliationSigmoid}
\alias{hydraulics_proportionDefoliationWeibull}
\title{Hydraulic-related defoliation}
\usage{
hydraulics_proportionDefoliationSigmoid(
  psiLeaf,
  P50,
  slope,
  PLC_crit = 0.88,
  P50_cv = 10
)

hydraulics_proportionDefoliationWeibull(
  psiLeaf,
  c,
  d,
  PLC_crit = 0.88,
  P50_cv = 10
)
}
\arguments{
\item{psiLeaf}{Leaf water potential (in MPa).}

\item{P50, slope}{Parameters of the Sigmoid function.}

\item{PLC_crit}{Critical leaf PLC corresponding to defoliation}

\item{P50_cv}{Coefficient of variation (in percent) of leaf P50, to describe the
variability in hydraulic vulnerability across crown leaves.}

\item{c, d}{Parameters of the Weibull function.}
}
\value{
The proportion of crown defoliation.
}
\description{
Functions to calculate the proportion of crown defoliation due to hydraulic disconnection.
}
\details{
The functions assume that crowns are made of a population of leaves whose
hydraulic vulnerability (i.e. the water potential corresponding to 50\% loss of conductance)
follows a Gaussian distribution centered on the input P50 and with a known coefficient of variation (\code{P50_cv}).
The slope parameter (or the c exponent in the case of a Weibull function) is considered constant.
Leaves are hydraulically disconnected, and shedded, when their embolism rate exceeds a critical value (\code{PLC_crit}).
}
\seealso{
\code{\link{hydraulics_conductancefunctions}}
}
\author{
Hervé Cochard, INRAE

Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\keyword{internal}
