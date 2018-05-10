\encoding{UTF-8}
\name{conductancefunctions}
\alias{hydraulics.psi2K}
\alias{hydraulics.K2Psi}
\alias{hydraulics.averagePsi}
\alias{hydraulics.psiCrit}
\alias{hydraulics.vulnerabilityCurvePlot}
\alias{hydraulics.vanGenuchtenConductance}
\alias{hydraulics.xylemConductance}
\alias{hydraulics.xylemPsi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydraulic confuctance functions
}
\description{
Set of functions used in the calculation of soil and plant hydraulic conductance.
}
\usage{
hydraulics.psi2K(psi, Psi_extract, ws = 3.0)
hydraulics.K2Psi(K, Psi_extract, ws = 3.0)
hydraulics.averagePsi(psi, v, c, d)
hydraulics.vulnerabilityCurvePlot(x, soil = NULL, type="leaf",  
                                  psiVec =  seq(-0.1, -8.0, by =-0.01), 
                                  relative = FALSE)   
hydraulics.psiCrit(c, d)
hydraulics.vanGenuchtenConductance(psi, krhizomax, n, alpha)
hydraulics.xylemConductance(psi, kxylemmax, c, d)
hydraulics.xylemPsi(kxylem, kxylemmax, c, d)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi}{A scalar (or a vector, depending on the function) with water potential (in MPa).}
  \item{K}{Whole-plant relative conductance (0-1).}
  \item{Psi_extract}{Soil water potential (in MPa) corresponding to 50\% whole-plant relative conductance.}
  \item{ws}{Exponent of the whole-plant relative conductance Weibull function.}
  \item{v}{Proportion of fine roots within each soil layer.}
  \item{krhizomax}{Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kxylemmax}{Maximum xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kleafmax}{Maximum leaf hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kstemmax}{Maximum stem xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{krootmax}{Maximum root xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{E}{Flow per surface unit.}
  \item{Emax}{Maximum flow per surface unit.}
  \item{psiUpstream}{Water potential upstream (in MPa). In a one-component model corresponds to soil potential. In a two-component model corresponds to the potential inside the roots.}
  \item{psiCav}{Minimum water potential (in MPa) experienced (for irreversible cavitation).}
  \item{minFlow}{Minimum flow in supply function.}
  \item{psiPlant}{Plant water potential (in MPa).}
  \item{psiSoil}{Soil water potential (in MPa). A scalar or a vector depending on the function.}
  \item{psiRhizo}{Soil water potential (in MPa) in the rhizosphere (root surface).}
  \item{psiStep}{Water potential precision (in MPa).}
  \item{psiTol}{Precision for water potential estimates (in MPa).}
  \item{psiIni}{Vector of initial water potential values (in MPa).}
  \item{psiMax}{Minimum (maximum in absolute value) water potential to be considered (in MPa).}
  \item{dE}{Increment of flow per surface unit.}
  \item{ETol}{Precision for water flow per surface unit.}
  \item{c, d}{Parameters of the Weibull function (generic xylem vulnerability curve).}
  \item{rootc, rootd}{Parameters of the Weibull function for roots (root xylem vulnerability curve).}
  \item{stemc, stemd}{Parameters of the Weibull function for stems (stem xylem vulnerability curve).}
  \item{leafc, leafd}{Parameters of the Weibull function for leaves (leaf vulnerability curve).}
  \item{n, alpha}{Parameters of the Van Genuchten function (rhizosphere vulnerability curve).}
  \item{nsoil, alphasoil}{Parameter vectors of the Van Genuchten function (rhizosphere vulnerability curve) with one value for each soil layer.}
  \item{allowNegativeFlux}{A boolean to indicate wether negative flux (i.e. from plant to soil) is allowed.}
  \item{maxNsteps}{Maximum number of steps in the construction of supply functions.}
  \item{ntrial}{Maximum number of steps in Newton-Raphson optimization.}
  \item{kxylem}{Xylem hydraulic conductance (defined as flow per surface unit and per pressure drop).}
  \item{x}{An object of class \code{\link{swbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{type}{Plot type. For \code{hydraulics.supplyFunctionPlot}, either \code{"E"}, \code{"Elayers"}, \code{"PsiStem"}, \code{"PsiRoot"}, \code{"PsiRhizo"} or \code{"dEdP"}). For \code{hydraulics.vulnerabilityCurvePlot}, either \code{"leaf"}, \code{"stem"}, \code{"root"}, or \code{"rhizosphere"}).}
  \item{psiVec}{Vector of water potential values to evaluate for the vulnerability curve.}
  \item{relative}{A flag to relativize vulnerability curves to the [0-1] interval.}
}
\details{
Details of the hydraulic model are given in a vignette. Function \code{hydraulics.vulnerabilityCurvePlot} draws a plot of the vulnerability curves for the given \code{soil} object and network properties of each plant cohort in \code{x}.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{hydraulics.psi2K}: Whole-plant relative conductance (0-1).}
\item{\code{hydraulics.K2Psi}: Soil water potential (in MPa) corresponding to the given whole-plant relative conductance value (inverse of \code{hydraulics.psi2K()}).}
\item{\code{hydraulics.averagePsi}: The average water potential (in MPa) across soil layers.}
\item{\code{hydraulics.vanGenuchtenConductance}: Rhizosphere conductance corresponding to an input water potential (soil vulnerability curve).}
\item{\code{hydraulics.xylemConductance}: Xylem conductance (flow rate per pressure drop) corresponding to an input water potential (plant vulnerability curve).}
\item{\code{hydraulics.xylemPsi}: Xylem water potential (in MPa) corresponding to an input xylem conductance (flow rate per pressure drop).}
}
}
\references{
Sperry, J. S., F. R. Adler, G. S. Campbell, and J. P. Comstock. 1998. Limitation of plant water use by rhizosphere and xylem conductance: results from a model. Plant, Cell \& Environment 21:347–359.

Sperry, J. S., and D. M. Love. 2015. What plant hydraulics can tell us about responses to climate-change droughts. New Phytologist 207:14–27.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{hydraulics.supplyFunctionPlot}}, \code{\link{hydraulics.maximumStemHydraulicConductance}}, \code{\link{swb}}, \code{\link{soil}}
}
\examples{
kstemmax = 4 # in mmol·m-2·s-1·MPa-1
stemc = 3 
stemd = -4 # in MPa

psiVec = seq(-0.1, -7.0, by =-0.01)

#Vulnerability curve
kstem = unlist(lapply(psiVec, hydraulics.xylemConductance, kstemmax, stemc, stemd))
plot(-psiVec, kstem, type="l",ylab="Xylem conductance (mmol·m-2·s-1·MPa-1)", 
xlab="Canopy pressure (-MPa)", lwd=1.5,ylim=c(0,kstemmax))
}
