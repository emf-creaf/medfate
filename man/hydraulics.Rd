\encoding{UTF-8}
\name{hydraulics}
\alias{hydraulics.psi2K}
\alias{hydraulics.K2Psi}
\alias{hydraulics.averagePsi}
\alias{hydraulics.averageRhizosphereResistancePercent}
\alias{hydraulics.E2psiXylem}
\alias{hydraulics.E2psiTwoElements}
\alias{hydraulics.E2psiVanGenuchten}
\alias{hydraulics.E2psiNetwork}
\alias{hydraulics.Ecrit}
\alias{hydraulics.EVanGenuchten}
\alias{hydraulics.EXylem}
\alias{hydraulics.findRhizosphereMaximumConductance}
\alias{hydraulics.maximumStemHydraulicConductance}
\alias{hydraulics.psiCrit}
\alias{hydraulics.regulatedPsiXylem}
\alias{hydraulics.regulatedPsiTwoElements}
\alias{hydraulics.supplyFunctionOneXylem}
\alias{hydraulics.supplyFunctionTwoElements}
\alias{hydraulics.supplyFunctionNetwork}
\alias{hydraulics.supplyFunctionPlot}
\alias{hydraulics.taperFactor}
\alias{hydraulics.vanGenuchtenConductance}
\alias{hydraulics.xylemConductance}
\alias{hydraulics.xylemPsi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydraulic submodel functions
}
\description{
Set of functions used in the implementation of Sperry and Love's (2015) hydraulic model.
}
\usage{
hydraulics.psi2K(psi, Psi_extract, ws = 3.0)
hydraulics.K2Psi(K, Psi_extract, ws = 3.0)
hydraulics.averagePsi(psi, v, c, d)
hydraulics.averageRhizosphereResistancePercent(krhizomax, n, alpha,
                                               krootmax, rootc, rootd,
                                               kstemmax, stemc, stemd, psiStep = -0.01)
hydraulics.E2psiXylem(E, psiUpstream, kxylemmax, c, d, psiCav = 0, psiStep = -0.01, 
                 psiMax = -10.0)
hydraulics.E2psiVanGenuchten(E, psiSoil, krhizomax, n, alpha, 
                             psiStep = -0.01, psiMax = -10.0)
hydraulics.E2psiTwoElements(E, psiSoil, krhizomax, kxylemmax, n, alpha, c, d, 
                            psiCav = 0, psiStep = -0.001, psiMax = -10.0)
hydraulics.E2psiNetwork(E, psiSoil, 
                        krhizomax, nsoil, alphasoil,
                        krootmax, rootc, rootd, 
                        kstemmax, stemc, stemd,
                        psiIni = as.numeric(c(0)), psiCav = 0,
                        psiStep = -0.001, psiMax = -10.0, ntrial = 10, 
                        psiTol = 0.0001, ETol = 1e-04)
hydraulics.Ecrit(psiUpstream, kxylemmax, c, d)
hydraulics.EVanGenuchten(psiRhizo, psiSoil, krhizomax, n, alpha, 
                         psiStep = -0.001, psiTol = 0.0001, allowNegativeFlux = TRUE)
hydraulics.EXylem(psiPlant, psiUpstream, kxylemmax, c, d, allowNegativeFlux = TRUE, 
                  psiCav = 0)
hydraulics.findRhizosphereMaximumConductance(averageResistancePercent, n, alpha,
                                             krootmax, rootc, rootd,
                                             kstemmax, stemc, stemd)        
hydraulics.maximumStemHydraulicConductance(xylemConductivity, Al2As, height, 
                                           taper = FALSE)
hydraulics.supplyFunctionOneXylem(psiSoil, v, 
                          kstemmax, stemc, stemd, psiCav = 0,
                          maxNsteps=200, psiStep = -0.001, psiMax = -10.0, dE=0.01)
hydraulics.supplyFunctionTwoElements(Emax, psiSoil, 
                          krhizomax, kxylemmax, n, alpha, c, d, 
                          psiCav = 0, dE = 0.1, psiMax = -10.0)
hydraulics.supplyFunctionNetwork(psiSoil, 
                          krhizomax, nsoil, alphasoil,
                          krootmax, rootc, rootd, 
                          kstemmax, stemc, stemd, psiCav = 0,
                          maxNsteps = 400, psiStep = -0.001, psiMax = -10.0, 
                          ntrial = 10, psiTol = 0.0001, ETol = 1e-04)
hydraulics.supplyFunctionPlot(soil, x, type = "E")   
hydraulics.taperFactor(height)
hydraulics.psiCrit(c, d)
hydraulics.regulatedPsiXylem(E, psiUpstream, kxylemmax, c, d, psiStep = -0.01)
hydraulics.regulatedPsiTwoElements(Emax, psiSoil, krhizomax, kxylemmax, n, alpha, 
                                   c, d, dE = 0.1, psiMax = -10.0)
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
  \item{krhizomax}{Maximum rhizosphere hydraulic conductance (defined as flow per surface unit and per pressure drop).}
  \item{kxylemmax}{Maximum xylem hydraulic conductance (defined as flow per surface unit and per pressure drop).}
  \item{kstemmax}{Maximum stem xylem hydraulic conductance (defined as flow per surface unit and per pressure drop).}
  \item{krootmax}{Maximum root xylem hydraulic conductance (defined as flow per surface unit and per pressure drop).}
  \item{E}{Flow per surface unit.}
  \item{Emax}{Maximum flow per surface unit.}
  \item{psiUpstream}{Water potential upstream (in MPa). In a one-component model corresponds to soil potential. In a two-component model corresponds to the potential inside the roots.}
  \item{psiCav}{Minimum water potential (in MPa) experienced (for irreversible cavitation).}
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
  \item{rootc, rootd}{Parameters of the Weibull function (root xylem vulnerability curve).}
  \item{stemc, stemd}{Parameters of the Weibull function (stem xylem vulnerability curve).}
  \item{n, alpha}{Parameters of the Van Genuchten function (rhizosphere vulnerability curve).}
  \item{nsoil, alphasoil}{Parameter vectors of the Van Genuchten function (rhizosphere vulnerability curve) with one value for each soil layer.}
  \item{allowNegativeFlux}{A boolean to indicate wether negative flux (i.e. from plant to soil) is allowed.}
  \item{averageResistancePercent}{Average (across water potential values) resistance percent of the rhizosphere, with respect to total resistance (rhizosphere + root xylem + stem xylem).}
  \item{xylemConductivity}{Xylem conductivity as flow per length of conduit and pressure drop (in kg·m-1·s-1·MPa-1).}
  \item{Al2As}{Leaf area to sapwood area (in m2·m-2).}
  \item{height}{Plant height (in cm).}
  \item{taper}{A boolean flag to indicate correction by taper of xylem conduits (Christoffersen et al. 2017)}
  \item{maxNsteps}{Maximum number of steps in the construction of supply functions.}
  \item{ntrial}{Maximum number of steps in Newton-Raphson optimization.}
  \item{kxylem}{Xylem hydraulic conductance (defined as flow per surface unit and per pressure drop).}
  \item{x}{An object of class \code{\link{swbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{type}{Plot type (either \code{"E"}, \code{"Elayers"}, \code{"PsiRhizo"} or \code{"dEdP"}).}
}
\details{
Details of the hydraulic model are given in a vignette. Function \code{hydraulics.supplyFunctionPlot} draws a plot of the supply function for the given \code{soil} object and network properties of each plant cohort in \code{x}.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{hydraulics.psi2K}: Whole-plant relative conductance (0-1).}
\item{\code{hydraulics.K2Psi}: Soil water potential (in MPa) corresponding to the given whole-plant relative conductance value (inverse of \code{hydraulics.psi2K()}).}
\item{\code{hydraulics.averagePsi}: The average water potential (in MPa) across soil layers.}
\item{\code{hydraulics.averageRhizosphereResistancePercent}: The average percentage of resistance due to the rhizosphere, calculated across water potential values.}
\item{\code{hydraulics.E2psiXylem}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the xylem supply function and given an upstream (soil or root) water potential.}
\item{\code{hydraulics.E2psiVanGenuchten}: The root water potential (in MPa) corresponding to the input flow, according to the rhizosphere supply function and given a soil water potential.}
\item{\code{hydraulics.E2psiTwoElements}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the rhizosphere and plant supply functions and given an input soil water potential.}
\item{\code{hydraulics.E2psiNetwork}: The rhizosphere, root crown and plant (leaf water potential (in MPa) corresponding to the input flow, according to the vulnerability curves of rhizosphere, root and stem elements in a network.}
\item{\code{hydraulics.Ecrit}: The critical flow according to the xylem supply function and given an input soil water potential.}
\item{\code{hydraulics.EVanGenuchten}: The flow (integral of the vulnerability curve) according to the rhizosphere supply function and given an input drop in water potential (soil and rhizosphere).}
\item{\code{hydraulics.EXylem}: The flow (integral of the vulnerability curve) according to the xylem supply function and given an input drop in water potential (rhizosphere and plant).}
\item{\code{hydraulics.findRhizosphereMaximumConductance}: The maximum rhizosphere conductance value given an average rhizosphere resistance and the vulnerability curves of rhizosphere, root and stem elements.}
\item{\code{hydraulics.supplyFunctionOneXylem}, \code{hydraulics.supplyFunctionTwoElements} and
\code{hydraulics.supplyFunctionNetwork}: A list with different numeric vectors with information of the two-element supply function:
\itemize{
\item{\code{E}: Flow values (supply values).}
\item{\code{FittedE}: Fitted flow values (for \code{hydraulics.supplyFunctionTwoElements}).}
\item{\code{Elayers}: Flow values across the roots of each soil layer (only for \code{hydraulics.supplyFunctionNetwork}).}
\item{\code{PsiRhizo}: Water potential values at the root surface (only for \code{hydraulics.supplyFunctionNetwork}).}
\item{\code{PsiRoot}: Water potential values inside the root crown (not for \code{hydraulics.supplyFunctionOneXylem}).}
\item{\code{PsiPlant}: Water potential values at the canopy (leaf).}
\item{\code{dEdP}: Derivatives of the supply function.}
}
}
\item{\code{hydraulics.supplyFunctionPlot}: A (hidden) list with the result of calling \code{hydraulics.supplyFunctionNetwork} for each cohort. }
\item{\code{hydraulics.psiCrit}: Critical water potential (in MPa) according to the xylem vulnerability curve. }
\item{\code{hydraulics.regulatedPsiXylem}: Plant water potential after regulation (one-element loss function) given an input water potential.}
\item{\code{hydraulics.regulatedPsiTwoElements}: Plant water potential after regulation (two-element loss function) given an input soil water potential.}
\item{\code{hydraulics.vanGenuchtenConductance}: Rhizosphere conductance corresponding to an input water potential (soil vulnerability curve).}
\item{\code{hydraulics.xylemConductance}: Xylem conductance (flow rate per pressure drop) corresponding to an input water potential (plant vulnerability curve).}
\item{\code{hydraulics.xylemPsi}: Xylem water potential (in MPa) corresponding to an input xylem conductance (flow rate per pressure drop).}
\item{\code{hydraulics.taperFunction}: Taper factor according to Savage et al. (2010).}
}
}
\references{
Sperry, J. S., and D. M. Love. 2015. What plant hydraulics can tell us about responses to climate-change droughts. New Phytologist 207:14–27.

Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R. Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A. Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and P. Meir. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 9: 4227–4255.

Savage, V. M., L. P. Bentley, B. J. Enquist, J. S. Sperry, D. D. Smith, P. B. Reich, and E. I. von Allmen. 2010. Hydraulic trade-offs and space filling enable better predictions of vascular structure and function in plants. Proceedings of the National Academy of Sciences of the United States of America 107:22722–7.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{swb}}, \code{\link{soil}}
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
