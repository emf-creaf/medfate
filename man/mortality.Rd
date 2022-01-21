\encoding{UTF-8}
\name{Mortality}
\alias{mortality_dailyProbability}

\title{Mortality}
\description{
A simple function to determine a daily mortality likelihood according to the value of a stress variable.
}
\usage{
mortality_dailyProbability(basalMortalityRate, stressValue, 
                           stressThreshold, minValue = 0.0,
                           exponent = 10.0)
}
\arguments{
  \item{basalMortalityRate}{Basal daily mortality rate per individual.}
  \item{stressValue}{Current value of the stress variable (lower values indicate stronger stress).}
  \item{stressThreshold}{Threshold to indicate that lower values increase mortality likelihood.}
  \item{minValue}{Minimum value of the stress variable (i.e. maximum stress), corresponding to probability of mortality equal to 1.}
  \item{exponent}{Coefficient modulating the steepness of the relationship.}
}
\value{
Returns a probability (between 0 and 1).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{growth}}
}
