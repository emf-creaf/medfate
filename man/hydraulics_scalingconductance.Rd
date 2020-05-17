\encoding{UTF-8}
\name{scalingconductance}
\alias{hydraulics_maximumSoilPlantConductance}
\alias{hydraulics_soilPlantResistances}
\alias{hydraulics_averageRhizosphereResistancePercent}
\alias{hydraulics_findRhizosphereMaximumConductance}
\alias{hydraulics_maximumStemHydraulicConductance}
\alias{hydraulics_maximumRootHydraulicConductance}
\alias{hydraulics_referenceConductivityHeightFactor}
\alias{hydraulics_taperFactorSavage}
\alias{hydraulics_terminalConduitRadius}
\alias{hydraulics_stemWaterCapacity}
\alias{hydraulics_leafWaterCapacity}
\title{
Scaling from conductivity to conductance
}
\description{
Functions used to scale from tissue conductivity to conductance of different elements of the continuum.
}
\usage{
hydraulics_maximumSoilPlantConductance(krhizomax, krootmax, 
                                       kstemmax, kleafmax)
hydraulics_soilPlantResistances(psiSoil, psiRhizo, 
                                psiStem, PLCstem, psiLeaf, 
                                krhizomax, n, alpha,
                                krootmax, rootc, rootd, 
                                kstemmax, stemc, stemd,
                                kleafmax, leafc, leafd)
hydraulics_averageRhizosphereResistancePercent(krhizomax, n, alpha,
                                               krootmax, rootc, rootd,
                                               kstemmax, stemc, stemd, 
                                               kleafmax, leafc, leafd, psiStep = -0.01)
hydraulics_findRhizosphereMaximumConductance(averageResistancePercent, n, alpha,
                                             krootmax, rootc, rootd,
                                             kstemmax, stemc, stemd,
                                             kleafmax, leafc, leafd) 
hydraulics_maximumRootHydraulicConductance(xylemConductivity, Al2As, 
                                           v, widths, depthWidthRatio = 1.0)
hydraulics_maximumStemHydraulicConductance(xylemConductivity, refheight, Al2As, height, 
                                           taper = FALSE)
hydraulics_referenceConductivityHeightFactor(refheight, height)
hydraulics_terminalConduitRadius(height)
hydraulics_taperFactorSavage(height)
hydraulics_stemWaterCapacity(Al2As, height, wd)
hydraulics_leafWaterCapacity(SLA, ld)
}

\arguments{
  \item{psiSoil}{Soil water potential (in MPa). A scalar or a vector depending on the function.}
  \item{psiRhizo}{Water potential (in MPa) in the rhizosphere (root surface).}
  \item{psiStem}{Water potential (in MPa) in the stem.}
  \item{psiLeaf}{Water potential (in MPa) in the leaf.}
  \item{PLCstem}{Percent loss of conductance (in \%) in the stem.}
  \item{v}{Proportion of fine roots within each soil layer.}
  \item{krhizomax}{Maximum rhizosphere hydraulic conductance (defined as flow per leaf surface unit and per pressure drop).}
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
  \item{refheight}{Reference plant height of measurement of xylem conductivity (in cm).}
  \item{taper}{A boolean flag to indicate correction by taper of xylem conduits (Christoffersen et al. 2017).}
  \item{widths}{Soil layer depths (in mm).}
  \item{depthWidthRatio}{Ratio between radius of the soil layer with the largest radius and maximum rooting depth.}
  \item{SLA}{Specific leaf area (mm2·mg-1).}
  \item{wd}{Wood density (g·cm-3).}
  \item{ld}{Leaf tissue density (g·cm-3).}
}
\details{
Details of the hydraulic model are given in a vignette.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{hydraulics_maximumSoilPlantConductance}: The maximum soil-plant conductance, in the same units as the input segment conductances.}
\item{\code{hydraulics_averageRhizosphereResistancePercent}: The average percentage of resistance due to the rhizosphere, calculated across water potential values.}
\item{\code{hydraulics_findRhizosphereMaximumConductance}: The maximum rhizosphere conductance value given an average rhizosphere resistance and the vulnerability curves of rhizosphere, root and stem elements.}
\item{\code{hydraulics_taperFactorSavage}: Taper factor according to Savage et al. (2010).}
}
}
\references{
Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R. Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A. Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and P. Meir. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 9: 4227–4255.

Savage, V. M., L. P. Bentley, B. J. Enquist, J. S. Sperry, D. D. Smith, P. B. Reich, and E. I. von Allmen. 2010. Hydraulic trade-offs and space filling enable better predictions of vascular structure and function in plants. Proceedings of the National Academy of Sciences of the United States of America 107:22722–7.

Olson, M.E., Anfodillo, T., Rosell, J.A., Petit, G., Crivellaro, A., Isnard, S., León-Gómez, C., Alvarado-Cárdenas, L.O., & Castorena, M. 2014. Universal hydraulics of the flowering plants: Vessel diameter scales with stem length across angiosperm lineages, habits and climates. Ecology Letters 17: 988–997.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
\code{\link{hydraulics_psi2K}}, \code{\link{hydraulics_supplyFunctionPlot}}, \code{\link{spwb}}, \code{\link{soil}}
}
\examples{
kstemmax = 4 # in mmol·m-2·s-1·MPa-1
stemc = 3 
stemd = -4 # in MPa

}
