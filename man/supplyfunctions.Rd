\encoding{UTF-8}
\name{supplyfunctions}
\alias{hydraulics.ECrit}
\alias{hydraulics.EXylem}
\alias{hydraulics.EVanGenuchten}
\alias{hydraulics.ECapacitance}
\alias{hydraulics.E2psiXylem}
\alias{hydraulics.E2psiVanGenuchten}
\alias{hydraulics.E2psiTwoElements}
\alias{hydraulics.E2psiBelowground}
\alias{hydraulics.E2psiAboveground}
\alias{hydraulics.E2psiAbovegroundCapacitance}
\alias{hydraulics.E2psiAbovegroundCapacitanceDisconnected}
\alias{hydraulics.E2psiNetwork}
\alias{hydraulics.E2psiNetworkCapacitance}
\alias{hydraulics.supplyFunctionOneXylem}
\alias{hydraulics.supplyFunctionTwoElements}
\alias{hydraulics.supplyFunctionThreeElements}
\alias{hydraulics.supplyFunctionBelowground}
\alias{hydraulics.supplyFunctionAboveground}
\alias{hydraulics.supplyFunctionAbovegroundCapacitance}
\alias{hydraulics.supplyFunctionNetwork}
\alias{hydraulics.supplyFunctionNetworkCapacitance}
\alias{hydraulics.supplyFunctionPlot}
\alias{hydraulics.regulatedPsiXylem}
\alias{hydraulics.regulatedPsiTwoElements}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydraulic supply functions
}
\description{
Set of functions used in the implementation of hydraulic supply functions (Sperry \& Love 2015).
}
\usage{
hydraulics.EXylem(psiPlant, psiUpstream, 
                  kxylemmax, c, d, allowNegativeFlux = TRUE, 
                  psiCav = 0)
hydraulics.EVanGenuchten(psiRhizo, psiSoil, krhizomax, 
                      n, alpha, l = 0.5)
hydraulics.ECapacitance(psi, psiPrev, PLCprev,
                      V, fapo, c,  d, 
                      pi0, eps,timestep)                         
hydraulics.ECrit(psiUpstream, kxylemmax, c, d, pCrit = 0.001)
hydraulics.E2psiXylem(E, psiUpstream, 
                      kxylemmax, c, d, psiCav = 0)
hydraulics.E2psiVanGenuchten(E, psiSoil, krhizomax, n, alpha, 
                      psiStep = -0.0001, psiMax = -10.0)
hydraulics.E2psiTwoElements(E, psiSoil, krhizomax, kxylemmax, n, alpha, c, d, 
                      psiCav = 0, psiStep = -1e-04, psiMax = -10.0)
hydraulics.E2psiBelowground(E, psiSoil, 
                      krhizomax, nsoil, alphasoil,
                      krootmax, rootc, rootd, 
                      psiIni = as.numeric(c(0)),
                      ntrial = 10, psiTol = 0.0001, ETol = 0.0001)
hydraulics.E2psiAboveground(E, psiRootCrown, 
                      kstemmax, stemc, stemd,
                      kleafmax, leafc, leafd,
                      PLCstem)  
hydraulics.E2psiAbovegroundCapacitance(E, psiRootCrown, 
                      psiStemPrev, PLCstem,
                      psiLeafPrev, 
                      kstemmax, stemc, stemd,
                      kleafmax, leafc, leafd,
                      Vsapwood, stemfapo, stempi0, stemeps,
                      Vleaf, leaffapo, leafpi0, leafeps,
                      tstep)
hydraulics.E2psiAbovegroundCapacitanceDisconnected(E,                           
                      psiStemPrev, PLCstem, RWCsympstemPrev, 
                      psiLeafPrev, RWCsympleafPrev,
                      kstemmax, stemc, stemd,
                      kleafmax, leafc, leafd,
                      Vsapwood, stemfapo, stempi0, stemeps,
                      Vleaf, leaffapo, leafpi0, leafeps,
                      klat,
                      tstep = 3600.0)                                       
hydraulics.E2psiNetwork(E, psiSoil, 
                  krhizomax, nsoil, alphasoil,
                  krootmax, rootc, rootd, 
                  kstemmax, stemc, stemd,
                  kleafmax, leafc, leafd,
                  PLCstem,
                  psiIni = as.numeric(c(0)),
                  ntrial = 10, 
                  psiTol = 0.0001, ETol = 0.0001)
hydraulics.E2psiNetworkCapacitance(E, psiSoil, 
                  psiStemPrev, PLCstem,
                  psiLeafPrev, 
                  krhizomax, nsoil, alphasoil,
                  krootmax, rootc, rootd, 
                  kstemmax, stemc, stemd,
                  kleafmax, leafc, leafd,
                  Vsapwood, stemfapo, stempi0, stemeps,
                  Vleaf, leaffapo, leafpi0, leafeps,
                  tstep = 3600.0,
                  psiIni = as.numeric(c(0)),
                  ntrial = 10, psiTol = 0.0001, ETol = 0.0001)

hydraulics.supplyFunctionOneXylem(psiSoil, v, 
                  kstemmax, stemc, stemd, psiCav = 0,
                  maxNsteps=200, dE=0.01)
hydraulics.supplyFunctionTwoElements(Emax, psiSoil, 
                  krhizomax, kxylemmax, n, alpha, c, d, 
                  psiCav = 0, dE = 0.1, psiMax = -10.0)
hydraulics.supplyFunctionThreeElements(Emax, psiSoil, 
                  krhizomax, kxylemmax, kleafmax,
                  n, alpha, stemc, stemd, leafc, leafd,
                  psiCav = 0, dE = 0.1, psiMax = -10.0)
                          
                          
hydraulics.supplyFunctionBelowground(psiSoil, 
                  krhizomax, nsoil, alphasoil,
                  krootmax, rootc, rootd, 
                  minFlow = 0.0, maxNsteps=400,  
                  ntrial = 10, psiTol = 0.0001, ETol = 0.0001, pCrit = 0.001)

hydraulics.supplyFunctionAboveground(Erootcrown, psiRootCrown, 
                  kstemmax, stemc, stemd,
                  kleafmax, leafc, leafd,
                  PLCstem)

hydraulics.supplyFunctionAbovegroundCapacitance(Erootcrown, psiRootCrown,
                  psiStemPrev, PLCstemPrev, 
                  psiLeafPrev, 
                  kstemmax, stemc, stemd,
                  kleafmax, leafc, leafd,
                  Vsapwood, stemfapo, stempi0, stemeps,
                  Vleaf, leaffapo, leafpi0, leafeps,
                  tstep = 3600.0)
                               
hydraulics.supplyFunctionNetwork(psiSoil, 
                  krhizomax, nsoil, alphasoil,
                  krootmax, rootc, rootd, 
                  kstemmax, stemc, stemd,
                  kleafmax, leafc, leafd,
                  PLCstem, minFlow = 0.0, maxNsteps=400,
                  ntrial = 200, psiTol = 0.0001, ETol = 0.0001, pCrit = 0.001)
                           
hydraulics.supplyFunctionPlot(x, soil, draw = TRUE, type = "E")   
hydraulics.regulatedPsiXylem(E, psiUpstream, kxylemmax, c, d, psiStep = -0.01)
hydraulics.regulatedPsiTwoElements(Emax, psiSoil, krhizomax, kxylemmax, n, alpha, 
                                   c, d, dE = 0.1, psiMax = -10.0)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{Proportion of fine roots within each soil layer.}
  \item{krhizomax}{Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kxylemmax}{Maximum xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kleafmax}{Maximum leaf hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kstemmax}{Maximum stem xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{krootmax}{Maximum root xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{klat}{Lateral hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{E}{Flow per surface unit.}
  \item{Emax}{Maximum flow per surface unit.}
  \item{Erootcrown}{Flow per surface unit at the root crown.}
  \item{psi}{Water potential (in MPa).}
  \item{psiPrev}{Water potential (in MPa) in the previous time step.}
  \item{psiStemPrev, psiLeafPrev}{Stem or leaf water potential (in MPa) in the previous time step.}
  \item{psiUpstream}{Water potential upstream (in MPa). In a one-component model corresponds to soil potential. In a two-component model corresponds to the potential inside the roots.}
  \item{psiCav}{Minimum water potential (in MPa) experienced (for irreversible cavitation).}
  \item{minFlow}{Minimum flow in supply function.}
  \item{psiPlant}{Plant water potential (in MPa).}
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
  \item{PLCstemPrev}{Previous proportion of loss conductance [0-1] in the stem.}
  \item{RWCsympstemPrev, RWCsympleafPrev}{Previous relative water content in the stem or leaf [0-1].}
  \item{V}{Capacity of the compartment per leaf area (in L/m2).}
  \item{Vleaf, Vsapwood}{Capacity of the compartment per leaf area (in L/m2) of leaves and sapwood.}
  \item{fapo}{Apoplastic fraction (proportion) in the segment.}
  \item{stemfapo, leaffapo}{Apoplastic fraction (proportion) in the leaf or stem.}
  \item{pi0}{Full turgor osmotic potential (MPa).}
  \item{stempi0,leafpi0}{Full turgor osmotic potential in the stem or the leaf (MPa).}
  \item{eps}{Bulk modulus of elasticity (MPa).}
  \item{stemeps,leafeps}{Bulk modulus of elasticity (MPa) in the stem or the leaf.}
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
  \item{type}{Plot type. For \code{hydraulics.supplyFunctionPlot}, either \code{"E"}, \code{"Elayers"}, \code{"PsiStem"}, \code{"PsiRoot"}, \code{"PsiRhizo"} or \code{"dEdP"}). For \code{hydraulics.vulnerabilityCurvePlot}, either \code{"leaf"}, \code{"stem"}, \code{"root"}, or \code{"rhizosphere"}).}
  \item{draw}{A flag to indicate whether the supply function should be drawn or just returned.}
  \item{tstep, timestep}{Time step in seconds.}
}
\details{
Details of the hydraulic model are given in a vignette. Function \code{hydraulics.supplyFunctionPlot} draws a plot of the supply function for the given \code{soil} object and network properties of each plant cohort in \code{x}. Function \code{hydraulics.vulnerabilityCurvePlot} draws a plot of the vulnerability curves for the given \code{soil} object and network properties of each plant cohort in \code{x}.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{hydraulics.E2psiXylem}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the xylem supply function and given an upstream (soil or root) water potential.}
\item{\code{hydraulics.E2psiVanGenuchten}: The root water potential (in MPa) corresponding to the input flow, according to the rhizosphere supply function and given a soil water potential.}
\item{\code{hydraulics.E2psiTwoElements}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the rhizosphere and plant supply functions and given an input soil water potential.}
\item{\code{hydraulics.E2psiNetwork}: The rhizosphere, root crown and plant (leaf water potential (in MPa) corresponding to the input flow, according to the vulnerability curves of rhizosphere, root and stem elements in a network.}
\item{\code{hydraulics.Ecrit}: The critical flow according to the xylem supply function and given an input soil water potential.}
\item{\code{hydraulics.EVanGenuchten}: The flow (integral of the vulnerability curve) according to the rhizosphere supply function and given an input drop in water potential (soil and rhizosphere).}
\item{\code{hydraulics.EXylem}: The flow (integral of the vulnerability curve) according to the xylem supply function and given an input drop in water potential (rhizosphere and plant).}
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
\item{\code{hydraulics.regulatedPsiXylem}: Plant water potential after regulation (one-element loss function) given an input water potential.}
\item{\code{hydraulics.regulatedPsiTwoElements}: Plant water potential after regulation (two-element loss function) given an input soil water potential.}
}
}
\references{
Sperry, J. S., F. R. Adler, G. S. Campbell, and J. P. Comstock. 1998. Limitation of plant water use by rhizosphere and xylem conductance: results from a model. Plant, Cell \& Environment 21:347–359.

Sperry, J. S., and D. M. Love. 2015. What plant hydraulics can tell us about responses to climate-change droughts. New Phytologist 207:14–27.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC, Catalonia, Spain
}
\seealso{
\code{\link{hydraulics.psi2K}}, \code{\link{hydraulics.maximumStemHydraulicConductance}}, \code{\link{spwb}}, \code{\link{soil}}
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
