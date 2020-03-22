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
                 p_change = c(-80,-40,-20,0,20,40,80), 
                 summary.fun = NULL, simplify=TRUE,...)
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
  \item{Setting \code{paramName = "Z50/Z95"} affects \code{below$V}, \code{below$VCroot_kmax} and \code{below$VGrhizo_kmax}.}
  \item{Modifying \code{LAI_live} also affects \code{LAI_expanded}.}
  \item{Modifying \code{VCroot_kmax} from \code{paramsTransp} affects both \code{paramsTransp$VCroot_kmax} and \code{below$VCroot_kmax}.}
  \item{Modifying \code{WaterStorage} affects simultaneously \code{paramsWaterStorage$Vleaf} and \code{paramsWaterStorage$Vsapwood}.}
  \item{Modifying \code{c} affects simultaneously \code{paramsTransp$VCleaf_c}, \code{paramsTransp$VCstem_c} and \code{paramsTransp$VCroot_c}.}
  \item{Modifying \code{d} affects simultaneously \code{paramsTransp$VCleaf_d}, \code{paramsTransp$VCstem_d} and \code{paramsTransp$VCroot_d}.}
  \item{Modifying \code{Plant_kmax} from \code{paramsTransp} affects \code{paramsTransp$VCleaf_kmax}, \code{paramsTransp$VCstem_kmax}, \code{paramsTransp$VCroot_kmax} and \code{below$VCroot_kmax}.}
  \item{Modifying \code{Al2As} from \code{paramsAnatomy} affects \code{paramsWaterStorage$Vsapwood}, \code{paramsTransp$VCstem_kmax}, \code{paramsTransp$VCroot_kmax} and \code{below$VCroot_kmax}.}
  \item{Setting \code{paramName = "Vmax298/Jmax298"} affects both \code{paramsTransp$Vmax298} and \code{paramsTransp$Jmax298}.}
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

