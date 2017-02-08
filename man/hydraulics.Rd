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
\alias{hydraulics.vanGenuchtenConductance}
\alias{hydraulics.xylemConductance}
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
hydraulics.E2psiXylem(E, psiUpstream, kxylemmax, c, d, psiStep = -0.01, 
                 psiMax = -10.0)
hydraulics.E2psiVanGenuchten(E, psiSoil, krhizomax, n, alpha, 
                             psiStep = -0.01, psiMax = -10.0)
hydraulics.E2psiTwoElements(E, psiSoil, krhizomax, kxylemmax, n, alpha, c, d, 
                            psiStep = -0.001, psiMax = -10.0)
hydraulics.E2psiNetwork(E, psiSoil, 
                        krhizomax, nsoil, alphasoil,
                        krootmax, rootc, rootd, 
                        kstemmax, stemc, stemd,
                        psiIni = as.numeric(c(0)),
                        psiStep = -0.001, psiMax = -10.0, ntrial = 10, 
                        psiTol = 0.0001, ETol = 0.001)
hydraulics.Ecrit(psiUpstream, kxylemmax, c, d)
hydraulics.EVanGenuchten(psiRhizo, psiSoil, krhizomax, n, alpha, 
                         psiStep = -0.001, psiTol = 0.0001, allowNegativeFlux = TRUE)
hydraulics.EXylem(psiPlant, psiUpstream, kxylemmax, c, d, allowNegativeFlux = TRUE)
hydraulics.findRhizosphereMaximumConductance(averageResistancePercent, n, alpha,
                                             krootmax, rootc, rootd,
                                             kstemmax, stemc, stemd)        
hydraulics.maximumStemHydraulicConductance(xylemConductivity, Al2As, height)
hydraulics.supplyFunctionOneXylem(psiSoil, v, 
                          kstemmax, stemc, stemd,
                          maxNsteps=200, psiStep = -0.001, psiMax = -10.0, dE=0.01)
hydraulics.supplyFunctionTwoElements(Emax, psiSoil, 
                          krhizomax, kxylemmax, n, alpha, c, d, 
                          dE = 0.1, psiMax = -10.0)
hydraulics.supplyFunctionNetwork(psiSoil, 
                          krhizomax, nsoil, alphasoil,
                          krootmax, rootc, rootd, 
                          kstemmax, stemc, stemd,
                          maxNsteps=200, psiStep = -0.001, psiMax = -10.0, 
                          ntrial = 10, psiTol = 0.0001, ETol = 0.001)
hydraulics.psiCrit(c, d)
hydraulics.regulatedPsiXylem(E, psiUpstream, kxylemmax, c, d, psiStep = -0.01)
hydraulics.regulatedPsiTwoElements(Emax, psiSoil, krhizomax, kxylemmax, n, alpha, 
                                   c, d, dE = 0.1, psiMax = -10.0)
hydraulics.vanGenuchtenConductance(psi, krhizomax, n, alpha)
hydraulics.xylemConductance(psi, kxylemmax, c, d)
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
  \item{maxNsteps}{Maximum number of steps in the construction of supply functions.}
  \item{ntrial}{Maximum number of steps in Newton-Raphson optimization.}
}
\details{
Details of the hydraulic model are given in a vignette.
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
\item{\code{hydraulics.psiCrit}: Critical water potential (in MPa) according to the xylem vulnerability curve. }
\item{\code{hydraulics.regulatedPsiXylem}: Plant water potential after regulation (one-element loss function) given an input water potential.}
\item{\code{hydraulics.regulatedPsiTwoElements}: Plant water potential after regulation (two-element loss function) given an input soil water potential.}
\item{\code{hydraulics.vanGenuchtenConductance}: Rhizosphere conductance corresponding to an input water potential (soil vulnerability curve).}
\item{\code{hydraulics.xylemConductance}: Xylem conductance (flow rate per pressure drop) corresponding to an input water potential (plant vulnerability curve).}
}
}
\references{
Sperry, J. S., and D. M. Love. 2015. What plant hydraulics can tell us about responses to climate-change droughts. New Phytologist 207:14–27.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{swb}}
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
