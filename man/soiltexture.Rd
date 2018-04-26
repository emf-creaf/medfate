\encoding{UTF-8}
\name{soil texture and hydraulics}
\alias{soil.psi2thetaSX}
\alias{soil.psi2thetaVG}
\alias{soil.psi}
\alias{soil.theta2psiSX}
\alias{soil.theta2psiVG}
\alias{soil.thetaFC}
\alias{soil.theta}
\alias{soil.waterFC}
\alias{soil.USDAType}
\alias{soil.vanGenuchtenParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil texture and hydraulics}
\description{
Functions \code{soil.psi2thetaSX} and \code{soil.theta2psiSX} calculate water potentials (psi) and water contents (theta) using texture data the formulae of Saxton et al. (1986) or Saxton & Rawls (2006) depending on whether organic matter is available. Functions code{soil.psi2thetaVG} and \code{soil.theta2psiVG} to the same calculations as before, but using the Van Genuchten - Mualem equations (\enc{Wösten}{Wosten} & van Genuchten 1988). Function \code{soil.USDAType} returns the USDA type for a given texture. Function \code{soil.vanGenuchtenParams} gives parameters for van Genuchten's conductance function for a given texture type (Leij et al. 1996). Functions \code{soil.waterFC} and \code{soil.thetaFC} calculate the water volume (in mm or as percent volume) of each soil layer at field capacity, according to a given pedotransfer model. Functions \code{soil.psi} and \code{soil.theta} return the current water potential and water content of the soil object, according to a given pedotransfer model.
}
\usage{
soil.psi2thetaSX(clay, sand, psi, om = NA)
soil.psi2thetaVG(n, alpha, theta_res, theta_sat, psi)
soil.theta2psiSX(clay, sand, theta, om = NA)
soil.theta2psiVG(n, alpha, theta_res, theta_sat, theta)
soil.USDAType(clay, sand)
soil.vanGenuchtenParams(soilType)
soil.psi(soil, model="SX")
soil.theta(soil, model="SX")
soil.waterFC(soil, model="SX")
soil.thetaFC(soil, model="SX")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clay}{Percentage of clay (in percent weight).}
  \item{sand}{Percentage of sand (in percent weight).}
  \item{n, alpha, theta_res, theta_sat}{Parameters of the Van Genuchten-Mualem model (m = 1 - 1/n).}
  \item{psi}{Water potential (in MPa).}
  \item{theta}{Relative water content (in percent volume).}
  \item{om}{Percentage of organic matter (optional, in percent weight).}
  \item{soilType}{A string indicating the soil type.}
  \item{soil}{Soil object (returned by function \code{\link{soil}}).}
  \item{model}{Either 'SX' or 'VG' for Saxton's or Van Genuchten's pedotransfer models.}
}
\value{
Functions \code{soil.psi2thetaSX} and \code{soil.psi2thetaVG} return the soil water potential (in MPa) from soil volumetric water content, and functions \code{soil.theta2psiSX} and \code{soil.theta2psiVG} do the reverse calculation returning water potential in MPa. Function \code{soil.USDAType} returns a string. Function \code{soil.vanGenuchtenParams} returns a vector with four parameter values.
}
\references{
Leij, F.J., Alves, W.J., Genuchten, M.T. Van, Williams, J.R., 1996. The UNSODA Unsaturated Soil Hydraulic Database User’s Manual Version 1.0.

Saxton, K.E., Rawls, W.J., Romberger, J.S., Papendick, R.I., 1986. Estimating generalized soil-water characteristics from texture. Soil Sci. Soc. Am. J. 50, 1031–1036.

Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic estimates by texture and organic matter for hydrologic solutions. Soil Sci. Soc. Am. J. 70, 1569. doi:10.2136/sssaj2005.0117

\enc{Wösten}{Wosten}, J.H.M., & van Genuchten, M.T. 1988. Using texture and other soil properties to predict the unsaturated soil hydraulic functions. Soil Science Society of America Journal 52: 1762–1770.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{ \code{\link{soil}}}
\examples{

# Plot Saxton's pedotransfer function
psi = seq(0, -6.0, by=-0.01)
plot(-psi, lapply(as.list(psi), FUN=soil.psi2thetaSX, clay=40, sand=10), 
     type="l", ylim=c(0,0.6),ylab="Water content (prop. volume)", xlab = "Soil water potential (-MPa)")

#Determine USDA soil texture type
type = soil.USDAType(clay=40, sand=10)
type

#Van Genuchten's params (according to the USDA type and UNSODA database)
vg = soil.vanGenuchtenParams(type)
vg

#Add Van Genuchten pedotransfer function
lines(-psi, lapply(as.list(psi), FUN=soil.psi2thetaVG, alpha=vg[1], n = vg[2], 
      theta_res = vg[3], theta_sat = vg[4]), lty=2) 

legend("topright", legend=c("Saxton", "Van Genuchten"), lty=c(1,2), bty="n")
}