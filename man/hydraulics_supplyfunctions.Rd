\encoding{UTF-8}
\name{supplyfunctions}
\alias{hydraulics_ECrit}
\alias{hydraulics_EXylem}
\alias{hydraulics_EVanGenuchten}
\alias{hydraulics_ECapacitance}
\alias{hydraulics_E2psiXylem}
\alias{hydraulics_E2psiXylemUp}
\alias{hydraulics_E2psiVanGenuchten}
\alias{hydraulics_E2psiTwoElements}
\alias{hydraulics_E2psiFineRootLeaf}
\alias{hydraulics_E2psiBelowground}
\alias{hydraulics_E2psiAboveground}
\alias{hydraulics_E2psiNetwork}
\alias{hydraulics_E2psiNetworkStem1}
\alias{hydraulics_supplyFunctionOneXylem}
\alias{hydraulics_supplyFunctionTwoElements}
\alias{hydraulics_supplyFunctionThreeElements}
\alias{hydraulics_supplyFunctionFineRootLeaf}
\alias{hydraulics_supplyFunctionBelowground}
\alias{hydraulics_supplyFunctionAboveground}
\alias{hydraulics_supplyFunctionNetwork}
\alias{hydraulics_supplyFunctionNetworkStem1}
\alias{hydraulics_supplyFunctionPlot}
\alias{hydraulics_regulatedPsiXylem}
\alias{hydraulics_regulatedPsiTwoElements}

\title{
Hydraulic supply functions
}
\description{
Set of functions used in the implementation of hydraulic supply functions (Sperry \& Love 2015).
}
\usage{
hydraulics_EXylem(psiPlant, psiUpstream, 
                  kxylemmax, c, d, allowNegativeFlux = TRUE, 
                  psiCav = 0)
hydraulics_EVanGenuchten(psiRhizo, psiSoil, krhizomax, 
                      n, alpha, l = 0.5)
hydraulics_ECapacitance(psi, psiPrev, PLCprev,
                        V, fapo, c, d, 
                        pi0, eps,
                        timestep)
hydraulics_ECrit(psiUpstream, kxylemmax, c, d, pCrit = 0.001)
hydraulics_E2psiXylem(E, psiUpstream, 
                      kxylemmax, c, d, psiCav = 0)
hydraulics_E2psiXylemUp(E, psiDownstream, 
                      kxylemmax, c, d, psiCav = 0)
hydraulics_E2psiVanGenuchten(E, psiSoil, krhizomax, n, alpha, 
                      psiStep = -0.0001, psiMax = -10.0)
hydraulics_E2psiTwoElements(E, psiSoil, krhizomax, kxylemmax, n, alpha, c, d, 
                      psiCav = 0, psiStep = -1e-04, psiMax = -10.0)
hydraulics_E2psiFineRootLeaf(E, psiFineRoot, 
                      krootmax, rootc, rootd,
                      kstemmax, stemc, stemd,
                      kleafmax, leafc, leafd,
                      PLCstem)
hydraulics_E2psiBelowground(E, psiSoil, 
                      krhizomax, nsoil, alphasoil,
                      krootmax, rootc, rootd, 
                      psiIni = as.numeric(c(0)),
                      ntrial = 10, psiTol = 0.0001, ETol = 0.0001)
hydraulics_E2psiAboveground(E, psiRootCrown, 
                      kstemmax, stemc, stemd,
                      kleafmax, leafc, leafd,
                      PLCstem)  
hydraulics_E2psiNetwork(E, psiSoil, 
                  krhizomax, nsoil, alphasoil,
                  krootmax, rootc, rootd, 
                  kstemmax, stemc, stemd,
                  kleafmax, leafc, leafd,
                  PLCstem,
                  psiIni = as.numeric(c(0)),
                  ntrial = 10, 
                  psiTol = 0.0001, ETol = 0.0001)
hydraulics_E2psiNetworkStem1(E, psiSoil, 
                  krhizomax, nsoil, alphasoil,
                  krootmax, rootc, rootd, 
                  kstemmax, stemc, stemd,
                  PLCstem,
                  psiIni = as.numeric(c(0)),
                  ntrial = 10, 
                  psiTol = 0.0001, ETol = 0.0001)
hydraulics_supplyFunctionOneXylem(psiSoil, v, 
                  kstemmax, stemc, stemd, psiCav = 0,
                  maxNsteps=200, dE=0.01)
hydraulics_supplyFunctionTwoElements(Emax, psiSoil, 
                  krhizomax, kxylemmax, n, alpha, c, d, 
                  psiCav = 0, dE = 0.1, psiMax = -10.0)
hydraulics_supplyFunctionThreeElements(Emax, psiSoil, 
                  krhizomax, kxylemmax, kleafmax,
                  n, alpha, stemc, stemd, leafc, leafd,
                  psiCav = 0, dE = 0.1, psiMax = -10.0)
                          
                          
hydraulics_supplyFunctionBelowground(psiSoil, 
                  krhizomax, nsoil, alphasoil,
                  krootmax, rootc, rootd, 
                  minFlow = 0.0, maxNsteps=400,  
                  ntrial = 10, psiTol = 0.0001, ETol = 0.0001, pCrit = 0.001)

hydraulics_supplyFunctionAboveground(Erootcrown, psiRootCrown, 
                  kstemmax, stemc, stemd,
                  kleafmax, leafc, leafd,
                  PLCstem)

hydraulics_supplyFunctionFineRootLeaf(psiFineRoot,
                    krootmax, rootc, rootd,
                    kstemmax, stemc, stemd,
                    kleafmax, leafc, leafd,
                    PLCstem,
                    minFlow = 0.0, maxNsteps=400, 
                    ETol = 0.0001, pCrit = 0.001)
hydraulics_supplyFunctionNetworkStem1(psiSoil, 
                    krhizomax, nsoil, alphasoil,
                    krootmax, rootc, rootd, 
                    kstemmax, stemc, stemd,
                    PLCstem,
                    minFlow = 0.0, maxNsteps=400, 
                    ntrial = 200, psiTol = 0.0001, ETol = 0.0001, pCrit = 0.001)
hydraulics_supplyFunctionNetwork(psiSoil, 
                  krhizomax, nsoil, alphasoil,
                  krootmax, rootc, rootd, 
                  kstemmax, stemc, stemd,
                  kleafmax, leafc, leafd,
                  PLCstem, minFlow = 0.0, maxNsteps=400,
                  ntrial = 200, psiTol = 0.0001, ETol = 0.0001, pCrit = 0.001)
                           
hydraulics_supplyFunctionPlot(x, soil, draw = TRUE, type = "E",
                              speciesNames = FALSE, ylim=NULL)   
hydraulics_regulatedPsiXylem(E, psiUpstream, kxylemmax, c, d, psiStep = -0.01)
hydraulics_regulatedPsiTwoElements(Emax, psiSoil, krhizomax, kxylemmax, n, alpha, 
                                   c, d, dE = 0.1, psiMax = -10.0)
}

\arguments{
  \item{v}{Proportion of fine roots within each soil layer.}
  \item{krhizomax}{Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kxylemmax}{Maximum xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kleafmax}{Maximum leaf hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kstemmax}{Maximum stem xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{krootmax}{Maximum root xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{E}{Flow per surface unit.}
  \item{Emax}{Maximum flow per surface unit.}
  \item{Erootcrown}{Flow per surface unit at the root crown.}
  \item{psi}{Water potential (in MPa).}
  \item{psiPrev}{Water potential (in MPa) in the previous time step.}
  \item{psiDownstream}{Water potential upstream (in MPa).}
  \item{psiUpstream}{Water potential upstream (in MPa). In a one-component model corresponds to soil potential. In a two-component model corresponds to the potential inside the roots.}
  \item{psiCav}{Minimum water potential (in MPa) experienced (for irreversible cavitation).}
  \item{minFlow}{Minimum flow in supply function.}
  \item{psiPlant}{Plant water potential (in MPa).}
  \item{psiFineRoot}{Water potential (in MPa) inside fine roots.}
  \item{psiSoil}{Soil water potential (in MPa). A scalar or a vector depending on the function.}
  \item{psiRhizo}{Soil water potential (in MPa) in the rhizosphere (root surface).}
  \item{psiRootCrown}{Soil water potential (in MPa) at the root crown.}
  \item{psiStep}{Water potential precision (in MPa).}
  \item{psiTol}{Precision for water potential estimates (in MPa).}
  \item{psiIni}{Vector of initial water potential values (in MPa).}
  \item{psiMax}{Minimum (maximum in absolute value) water potential to be considered (in MPa).}
  \item{pCrit}{Critical water potential (in MPa).}
  \item{PLCstem}{Proportion of loss conductance in the stem [0-1].}
  \item{PLCprev}{Previous proportion of loss conductance [0-1].}
  \item{V}{Capacity of the compartment per leaf area (in L/m2).}
  \item{fapo}{Apoplastic fraction (proportion) in the segment.}
  \item{pi0}{Full turgor osmotic potential (MPa).}
  \item{eps}{Bulk modulus of elasticity (MPa).}
  \item{dE}{Increment of flow per surface unit.}
  \item{ETol}{Precision for water flow per surface unit.}
  \item{c, d}{Parameters of the Weibull function (generic xylem vulnerability curve).}
  \item{rootc, rootd}{Parameters of the Weibull function for roots (root xylem vulnerability curve).}
  \item{stemc, stemd}{Parameters of the Weibull function for stems (stem xylem vulnerability curve).}
  \item{leafc, leafd}{Parameters of the Weibull function for leaves (leaf vulnerability curve).}
  \item{n, alpha, l}{Parameters of the Van Genuchten function (rhizosphere vulnerability curve).}
  \item{nsoil, alphasoil}{Parameter vectors of the Van Genuchten function (rhizosphere vulnerability curve) with one value for each soil layer.}
  \item{allowNegativeFlux}{A boolean to indicate wether negative flux (i.e. from plant to soil) is allowed.}
  \item{maxNsteps}{Maximum number of steps in the construction of supply functions.}
  \item{ntrial}{Maximum number of steps in Newton-Raphson optimization.}
  \item{x}{An object of class \code{\link{spwbInput}}.}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{type}{Plot type for \code{hydraulics_supplyFunctionPlot}, either \code{"E"}, \code{"ERhizo"}, \code{"psiStem"}, \code{"psiRoot"}, \code{"psiRhizo"} or \code{"dEdP"}).}
  \item{draw}{A flag to indicate whether the supply function should be drawn or just returned.}
  \item{speciesNames}{A flag to indicate the use of species names instead of cohort names in plots.}
  \item{ylim}{Graphical parameter to override function defaults.}
  \item{timestep}{Time step in seconds.}
}
\details{
Function \code{hydraulics_supplyFunctionPlot} draws a plot of the supply function for the given \code{soil} object and network properties of each plant cohort in \code{x}. Function \code{hydraulics_vulnerabilityCurvePlot} draws a plot of the vulnerability curves for the given \code{soil} object and network properties of each plant cohort in \code{x}.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{hydraulics_E2psiXylem}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the xylem supply function and given an upstream (soil or root) water potential.}
\item{\code{hydraulics_E2psiVanGenuchten}: The root water potential (in MPa) corresponding to the input flow, according to the rhizosphere supply function and given a soil water potential.}
\item{\code{hydraulics_E2psiTwoElements}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the rhizosphere and plant supply functions and given an input soil water potential.}
\item{\code{hydraulics_E2psiNetwork}: The rhizosphere, root crown and plant (leaf water potential (in MPa) corresponding to the input flow, according to the vulnerability curves of rhizosphere, root and stem elements in a network.}
\item{\code{hydraulics_Ecrit}: The critical flow according to the xylem supply function and given an input soil water potential.}
\item{\code{hydraulics_EVanGenuchten}: The flow (integral of the vulnerability curve) according to the rhizosphere supply function and given an input drop in water potential (soil and rhizosphere).}
\item{\code{hydraulics_EXylem}: The flow (integral of the vulnerability curve) according to the xylem supply function and given an input drop in water potential (rhizosphere and plant).}
\item{\code{hydraulics_supplyFunctionOneXylem}, \code{hydraulics_supplyFunctionTwoElements} and
\code{hydraulics_supplyFunctionNetwork}: A list with different numeric vectors with information of the two-element supply function:
\itemize{
\item{\code{E}: Flow values (supply values).}
\item{\code{FittedE}: Fitted flow values (for \code{hydraulics_supplyFunctionTwoElements}).}
\item{\code{Elayers}: Flow values across the roots of each soil layer (only for \code{hydraulics_supplyFunctionNetwork}).}
\item{\code{PsiRhizo}: Water potential values at the root surface (only for \code{hydraulics_supplyFunctionNetwork}).}
\item{\code{PsiRoot}: Water potential values inside the root crown (not for \code{hydraulics_supplyFunctionOneXylem}).}
\item{\code{PsiPlant}: Water potential values at the canopy (leaf).}
\item{\code{dEdP}: Derivatives of the supply function.}
}
}
\item{\code{hydraulics_supplyFunctionPlot}: If \code{draw = FALSE} a list with the result of calling \code{hydraulics_supplyFunctionNetwork} for each cohort. }
\item{\code{hydraulics_regulatedPsiXylem}: Plant water potential after regulation (one-element loss function) given an input water potential.}
\item{\code{hydraulics_regulatedPsiTwoElements}: Plant water potential after regulation (two-element loss function) given an input soil water potential.}
}
}
\references{
Sperry, J. S., F. R. Adler, G. S. Campbell, and J. P. Comstock. 1998. Limitation of plant water use by rhizosphere and xylem conductance: results from a model. Plant, Cell \& Environment 21:347–359.

Sperry, J. S., and D. M. Love. 2015. What plant hydraulics can tell us about responses to climate-change droughts. New Phytologist 207:14–27.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF.
}
\seealso{
\code{\link{hydraulics_psi2K}}, \code{\link{hydraulics_maximumStemHydraulicConductance}}, \code{\link{spwb}}, \code{\link{soil}}
}
\examples{
kstemmax = 4 # in mmol·m-2·s-1·MPa-1
stemc = 3 
stemd = -4 # in MPa

psiVec = seq(-0.1, -7.0, by =-0.01)

#Vulnerability curve
kstem = unlist(lapply(psiVec, hydraulics_xylemConductance, kstemmax, stemc, stemd))
plot(-psiVec, kstem, type="l",ylab="Xylem conductance (mmol·m-2·s-1·MPa-1)", 
xlab="Canopy pressure (-MPa)", lwd=1.5,ylim=c(0,kstemmax))
}
