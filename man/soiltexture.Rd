\encoding{UTF-8}
\name{soil texture}
\alias{soil.psi2theta}
\alias{soil.theta2psi}
\alias{soil.USDAType}
\alias{soil.vanGenuchtenParams}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil texture functions}
\description{
Functions \code{soil.psi2theta} and \code{soil.theta2psi} calculate water potentials (psi) and relative water contents (theta) using texture data the formulae of Saxton (1986). Function \code{soil.USDAType} returns the USDA type for a given texture. Function \code{soil.vanGenuchtenParams} gives parameters for van Genuchten's conductance function for a given texture type (Leij et al. 1996).
}
\usage{
soil.psi2theta(clay, sand, psi)
soil.theta2psi(clay, sand, theta)
soil.USDAType(clay, sand)
soil.vanGenuchtenParams(soilType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clay}{Percentage of clay.}
  \item{sand}{Percentage of sand.}
  \item{psi}{Water potential (in MPa).}
  \item{theta}{Relative water content (in percent)}
  \item{soilType}{A string indicating the soil type.}
}
\value{
Function \code{soil.psi2theta} returns the soil water potential (in MPa) from soil relative water content, and the function \code{soil.theta2psi} does the reverse calculation. Function \code{soil.USDAType} returns a string. Function \code{soil.vanGenuchtenParams} returns a vector with two parameter values.
}
\references{
Saxton, K.E., Rawls, W.J., Romberger, J.S., Papendick, R.I., 1986. Estimating generalized soil-water characteristics from texture. Soil Sci. Soc. Am. J. 50, 1031–1036.

Leij, F.J., Alves, W.J., Genuchten, M.T. Van, Williams, J.R., 1996. The UNSODA Unsaturated Soil Hydraulic Database User’s Manual Version 1.0.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{ \code{\link{soil}}}
\examples{
type = soil.USDAType(clay=10, sand=70)
type
soil.vanGenuchtenParams(type)

psi = seq(0, -6.0, by=-0.01)
plot(-psi, lapply(as.list(psi), FUN=soil.psi2theta, clay=70, sand=10), 
     type="l", ylim=c(0,0.6),ylab="Theta")
lines(-psi, lapply(as.list(psi), FUN=soil.psi2theta, clay=10, sand=70), lty=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.