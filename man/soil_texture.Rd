\encoding{UTF-8}
\name{soil texture and hydraulics}
\alias{soil_psi2thetaSX}
\alias{soil_psi2thetaVG}
\alias{soil_psi}
\alias{soil_theta2psiSX}
\alias{soil_theta2psiVG}
\alias{soil_thetaFC}
\alias{soil_thetaWP}
\alias{soil_thetaSAT}
\alias{soil_thetaSATSX}
\alias{soil_theta}
\alias{soil_water}
\alias{soil_waterFC}
\alias{soil_waterWP}
\alias{soil_waterExtractable}
\alias{soil_waterSAT}
\alias{soil_USDAType}
\alias{soil_saturatedConductivitySX}
\alias{soil_vanGenuchtenParamsCarsel}
\alias{soil_vanGenuchtenParamsToth}
\alias{soil_waterTableDepth}
\alias{soil_retentionCurvePlot}
\alias{soil_rockWeight2Volume}
\title{Soil texture and hydraulics}
\description{
Collection of functions relating soil texture with soil hydraulics and soil water content (see description of each function in details).
}
\usage{
soil_psi2thetaSX(clay, sand, psi, om = NA)
soil_psi2thetaVG(n, alpha, theta_res, theta_sat, psi)
soil_theta2psiSX(clay, sand, theta, om = NA)
soil_theta2psiVG(n, alpha, theta_res, theta_sat, theta)
soil_saturatedConductivitySX(clay, sand, om = NA, mmol = TRUE)
soil_USDAType(clay, sand)
soil_vanGenuchtenParamsCarsel(soilType)
soil_vanGenuchtenParamsToth(clay, sand, om, bd, topsoil)
soil_psi(soil, model="SX")
soil_theta(soil, model="SX")
soil_water(soil, model="SX")
soil_waterFC(soil, model="SX")
soil_waterWP(soil, model="SX")
soil_waterSAT(soil, model="SX")
soil_waterExtractable(soil, model="SX", minPsi = -5.0)
soil_thetaFC(soil, model="SX")
soil_thetaWP(soil, model="SX")
soil_thetaSAT(soil, model="SX")
soil_thetaSATSX(clay, sand, om = NA)
soil_waterTableDepth(soil, model="SX")
soil_rockWeight2Volume(pWeight, bulkDensity, rockDensity = 2.3)
soil_retentionCurvePlot(soil, model="SX", layer = 1, 
                        psi = seq(0, -6.0, by=-0.01),
                        relative = TRUE, to = "SAT")
}
\arguments{
  \item{clay}{Percentage of clay (in percent weight).}
  \item{sand}{Percentage of sand (in percent weight).}
  \item{n, alpha, theta_res, theta_sat}{Parameters of the Van Genuchten-Mualem model (m = 1 - 1/n).}
  \item{psi}{Water potential (in MPa).}
  \item{theta}{Relative water content (in percent volume).}
  \item{om}{Percentage of organic matter (optional, in percent weight).}
  \item{mmol}{Boolean flag to indicate that saturated conductivity units should be returned in mmol/m/s/MPa. If \code{mmol = FALSE} then units are cm/day.}
  \item{bd}{Bulk density (in g/cm3).}
  \item{topsoil}{A boolean flag to indicate topsoil layer.}
  \item{soilType}{A string indicating the soil type.}
  \item{soil}{Soil object (returned by function \code{\link{soil}}).}
  \item{model}{Either 'SX' or 'VG' for Saxton's or Van Genuchten's water retention models; or 'both' to plot both retention models.}
  \item{layer}{Soil layer to be plotted.}
  \item{relative}{Boolean flag to indicate that retention curve should be relative to field capacity or saturation.}
  \item{to}{Either 'SAT' (saturation) or 'FC' (field capacity).}
  \item{minPsi}{Minimum water potential (in MPa) to calculate the amount of extractable water.}
  \item{pWeight}{Percentage of corresponding to rocks, in weight.}
  \item{bulkDensity}{Bulk density of the soil fraction (g/cm3).}
  \item{rockDensity}{Rock density (g/cm3).}
}
\details{
\itemize{
\item{\code{soil_psi2thetaSX()} and \code{soil_theta2psiSX()} calculate water potentials (MPa) and water contents (theta) using texture data the formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available.}
\item{\code{soil_psi2thetaVG()} and \code{soil_theta2psiVG()} to the same calculations as before, but using the Van Genuchten - Mualem equations (\enc{Wösten}{Wosten} & van Genuchten 1988). }
\item{\code{soil_saturatedConductivitySX()} returns the saturated conductivity of the soil (in cm/day or mmol/m/s/MPa), estimated from formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available.}
\item{\code{soil_USDAType()} returns the USDA type (a string) for a given texture.}
\item{\code{soil_vanGenuchtenParamsCarsel()} gives parameters for van Genuchten-Mualem equations (alpha, n, theta_res and theta_sat, where alpha is in MPa-1) for a given texture type (Leij et al. 1996) }
\item{\code{soil_vanGenuchtenParamsToth()} gives parameters for van Genuchten-Mualem equations (alpha, n, theta_res and theta_sat, where alpha is in MPa-1) for a given texture, organic matter and bulk density (Toth et al. 2015).}
\item{\code{soil_psi()} returns the water potential (MPa) of each soil layer, according to its water retention model.}
\item{\code{soil_theta()} returns the moisture content (as percent of soil volume) of each soil layer, according to its water retention model.}
\item{\code{soil_water()} returns the water volume (mm) of each soil layer, according to its water retention model.}
\item{\code{soil_waterExtractable()} returns the water volume (mm) extractable from the soil according to its water retention curves and up to a given soil water potential.}
\item{\code{soil_waterFC()} and \code{soil_thetaFC()} calculate the water volume (in mm) and moisture content (as percent of soil volume) of each soil layer at field capacity, respectively.}
\item{\code{soil_waterWP()} and \code{soil_thetaWP()} calculate the water volume (in mm) and moisture content (as percent of soil volume) of each soil layer at wilting point (-1.5 MPa), respectively. }
\item{\code{soil_waterSAT()}, \code{soil_thetaSATSX()} and \code{soil_thetaSAT()} calculate the saturated water volume (in mm) and moisture content (as percent of soil volume) of each soil layer.}
\item{\code{soil_waterTableDepth()} returns water table depth in mm from surface.}
\item{\code{soil_rockWeight2Volume()} transforms rock percentage from weight to volume basis.}
\item{\code{soil_retentionCurvePlot()} allows ploting the water retention curve of a given soil layer.}
}
}
\value{
Depends on the function (see details).
}
\references{
Leij, F.J., Alves, W.J., Genuchten, M.T. Van, Williams, J.R., 1996. The UNSODA Unsaturated Soil Hydraulic Database User’s Manual Version 1.0.

Saxton, K.E., Rawls, W.J., Romberger, J.S., Papendick, R.I., 1986. Estimating generalized soil-water characteristics from texture. Soil Sci. Soc. Am. J. 50, 1031–1036.

Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic estimates by texture and organic matter for hydrologic solutions. Soil Sci. Soc. Am. J. 70, 1569. doi:10.2136/sssaj2005.0117

\enc{Wösten}{Wosten}, J.H.M., & van Genuchten, M.T. 1988. Using texture and other soil properties to predict the unsaturated soil hydraulic functions. Soil Science Society of America Journal 52: 1762–1770.

\enc{Tóth}{Toth}, B., Weynants, M., Nemes, A., \enc{Makó}{Mako}, A., Bilas, G., & \enc{Tóth}{Toth}, G. 2015. New generation of hydraulic pedotransfer functions for Europe. European Journal of Soil Science 66: 226–238.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{ \code{\link{soil}}}
\examples{

#Determine USDA soil texture type
type = soil_USDAType(clay=40, sand=10)
type

#Van Genuchten's params (bulk density = 1.3 g/cm)
vg = soil_vanGenuchtenParamsToth(40,10,1,1.3,TRUE)
vg

# Initialize soil object with default params
s = soil(defaultSoilParams())

# Plot Saxton's and Van Genuchten's water retention curves
soil_retentionCurvePlot(s, model="both")

}
