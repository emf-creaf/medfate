\name{spwb_sensitivity}
\alias{spwb_sensitivity}
\title{
Sensitivity analysis for soil plant water balance simulations
}
\description{
Performs a set of calls to \code{\link{spwb}} with the aim to determine the sensitivity to particular parameters.
}
\usage{
spwb_sensitivity(x, soil, meteo, 
                 paramType = "above", paramName = "LAI_live", cohort = NA,
                 p_change = c(-80,-40,-20,0,20,40,80), summary.fun = NULL, simplify=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{meteo}{A data frame with daily meteorological data series (see \code{\link{spwb}}).}
  \item{paramType}{Data frame of \code{x} to modify.}
  \item{paramName}{Name of the parameter to modify.}
  \item{cohort}{Integer with the cohort to modify (if \code{NA} parameter values of all cohort are modified).}
  \item{p_change}{Numerical vector with percentages of change.}
  \item{summary.fun}{Summary function to be applied to the results of each simulation. }
  \item{simplify}{Whether the result of \code{summary.fun} should be simplified (see \code{\link{sapply}}). }
  \item{...}{Additional parameters to function \code{\link{spwb}}.}
}
\details{
Due to parameter dependence, modifying some parameters affects others:
\itemize{
  \item{Modifying \code{LAI_live} also affects \code{LAI_expanded}.}
  \item{Modifying \code{VCroot_kmax} from \code{paramsTransp} affects both \code{paramsTransp$VCroot_kmax} and \code{below$VCroot_kmax}.}
}

}
\value{
If \code{summary.fun = NULL} the function returns a list whose elements are the result of calling \code{\link{spwb}}. Otherwise, the function applies \code{summary.fun} to each simulation result and returns these summaries (actually, a call to \code{\link{sapply}} is done).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}

\seealso{
\code{\link{spwb}}, \code{\link{summary.spwb}}
}

