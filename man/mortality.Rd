\encoding{UTF-8}
\name{Mortality}
\alias{mortality_dailyProbability}

\title{Mortality}
\description{
A simple function to determine a daily mortality likelihood according to the value of a stress variable.
}
\usage{
mortality_dailyProbability(mortalityBaselineRate, 
                           stressValue, stressThreshold, 
                           allowStress = TRUE,
                           minValue = 0.0, slope = 1.0)
}
\arguments{
  \item{mortalityBaselineRate}{Baseline mortality rate at the annual scale.}
  \item{stressValue}{Current value of the stress variable (lower values indicate stronger stress).}
  \item{stressThreshold}{Threshold to indicate that lower values increase mortality likelihood.}
  \item{allowStress}{A boolean function to activate stress effects.}
  \item{minValue}{Minimum value of the stress variable (i.e. maximum stress), corresponding to probability of mortality equal to 1.}
  \item{slope}{Slope coefficient modulating how fast probability increases after the stress threshold has been reached.}
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
