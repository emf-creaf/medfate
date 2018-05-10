\encoding{UTF-8}
\name{scalingconductance}
\alias{hydraulics.averageRhizosphereResistancePercent}
\alias{hydraulics.findRhizosphereMaximumConductance}
\alias{hydraulics.maximumStemHydraulicConductance}
\alias{hydraulics.maximumRootHydraulicConductance}
\alias{hydraulics.referenceConductivityHeightFactor}
\alias{hydraulics.taperFactorSavage}
\alias{hydraulics.terminalConduitRadius}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scaling from conductivity to conductance
}
\description{
Functions used to scale from tissue conductivity to conductance of different elements of the continuum.
}
\usage{
hydraulics.averageRhizosphereResistancePercent(krhizomax, n, alpha,
                                               krootmax, rootc, rootd,
                                               kstemmax, stemc, stemd, 
                                               kleafmax, leafc, leafd, psiStep = -0.01)
hydraulics.findRhizosphereMaximumConductance(averageResistancePercent, n, alpha,
                                             krootmax, rootc, rootd,
                                             kstemmax, stemc, stemd,
                                             kleafmax, leafc, leafd) 
hydraulics.maximumRootHydraulicConductance(xylemConductivity, Al2As, 
                                           v, widths, depthWidthRatio = 1.0)
hydraulics.maximumStemHydraulicConductance(xylemConductivity, refheight, Al2As, height, 
                                           angiosperm = TRUE, taper = FALSE)
hydraulics.referenceConductivityHeightFactor(refheight, height)
hydraulics.terminalConduitRadius(height)
hydraulics.taperFactorSavage(height)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{Proportion of fine roots within each soil layer.}
  \item{krhizomax}{Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kxylemmax}{Maximum xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kleafmax}{Maximum leaf hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{kstemmax}{Maximum stem xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{krootmax}{Maximum root xylem hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
  \item{psiStep}{Water potential precision (in MPa).}
  \item{rootc, rootd}{Parameters of the Weibull function for roots (root xylem vulnerability curve).}
  \item{stemc, stemd}{Parameters of the Weibull function for stems (stem xylem vulnerability curve).}
  \item{leafc, leafd}{Parameters of the Weibull function for leaves (leaf vulnerability curve).}
  \item{n, alpha}{Parameters of the Van Genuchten function (rhizosphere vulnerability curve).}
  \item{averageResistancePercent}{Average (across water potential values) resistance percent of the rhizosphere, with respect to total resistance (rhizosphere + root xylem + stem xylem).}
  \item{xylemConductivity}{Xylem conductivity as flow per length of conduit and pressure drop (in kg·m-1·s-1·MPa-1).}
  \item{Al2As}{Leaf area to sapwood area (in m2·m-2).}
  \item{height}{Plant height (in cm).}
  \item{refheight}{Reference plant height (in cm).}
  \item{angiosperm}{A boolean flag to indicate an angiosperm species.}
  \item{taper}{A boolean flag to indicate correction by taper of xylem conduits (Christoffersen et al. 2017).}
  \item{widths}{Soil layer depths (in mm).}
}
\details{
Details of the hydraulic model are given in a vignette.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{hydraulics.averageRhizosphereResistancePercent}: The average percentage of resistance due to the rhizosphere, calculated across water potential values.}
\item{\code{hydraulics.findRhizosphereMaximumConductance}: The maximum rhizosphere conductance value given an average rhizosphere resistance and the vulnerability curves of rhizosphere, root and stem elements.}
\item{\code{hydraulics.taperFactorSavage}: Taper factor according to Savage et al. (2010).}
}
}
\references{
Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R. Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A. Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and P. Meir. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 9: 4227–4255.

Savage, V. M., L. P. Bentley, B. J. Enquist, J. S. Sperry, D. D. Smith, P. B. Reich, and E. I. von Allmen. 2010. Hydraulic trade-offs and space filling enable better predictions of vascular structure and function in plants. Proceedings of the National Academy of Sciences of the United States of America 107:22722–7.

Olson, M.E., Anfodillo, T., Rosell, J.A., Petit, G., Crivellaro, A., Isnard, S., León-Gómez, C., Alvarado-Cárdenas, L.O., & Castorena, M. 2014. Universal hydraulics of the flowering plants: Vessel diameter scales with stem length across angiosperm lineages, habits and climates. Ecology Letters 17: 988–997.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{hydraulics.psi2K}}, \code{\link{hydraulics.supplyFunctionPlot}}, \code{\link{swb}}, \code{\link{soil}}
}
\examples{
kstemmax = 4 # in mmol·m-2·s-1·MPa-1
stemc = 3 
stemd = -4 # in MPa

}
