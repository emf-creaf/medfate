\encoding{UTF-8}
\name{tissuemoisture}
\alias{moisture.symplasticRWC}
\alias{moisture.apoplasticRWC}
\alias{moisture.leafRWC}
\alias{moisture.branchRWC}
\alias{moisture.fineFuelRWC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tissue moisture functions
}
\description{
Set of functions used to calculate tissue moisture from water potential.
}
\usage{
moisture.symplasticRWC(psi, pi0, epsilon)
moisture.apoplasticRWC(psi, c, d, cellWallFraction = 0.07)
moisture.leafRWC(psi, pi0, epsilon, af)
moisture.branchRWC(psi, wd, c, d, af = 0.80)
moisture.fineFuelRWC(psi, leaf_pi0, leaf_eps, leaf_af, 
                     wd, c, d, r635)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi}{Tissue water potential (MPa).}
  \item{pi0}{Full turgor osmotic potential (MPa).}
  \item{epsilon}{Bulk modulus of elasticity (MPa).}
  \item{c,d}{Parameters of the xylem vulnerability curve.}
  \item{cellWallFraction}{Fraction of relative water content in the cell wall.}
  \item{af}{Apoplastic fraction (proportion).}
  \item{wd}{Wood density (g/cm3).}
  \item{leaf_pi0, leaf_eps, leaf_af}{Full turgor osmotic potential, elasticity modulus and apoplastic fraction of leaf tissue.}
  \item{r635}{Ratio of biomass of foliar + small branches (<6.35) to foliar biomass.}
}
\details{
Details of the tissue moisture calculations are given in a vignette. 
}
\value{
Values returned for each function are:
\itemize{
\item{\code{moisture.symplasticRWC}: Relative water content of the symplastic fraction.}
\item{\code{moisture.apoplasticRWC}: Relative water content of the apoplastic fraction.}
\item{\code{moisture.leafRWC}: Leaf relative water content.}
\item{\code{moisture.branchRWC}: Branch relative water content.}
\item{\code{moisture.fineFuelRWC}: Fine fuel (i.e. leaves + small branches) relative water content.}
}
}
\references{
Bartlett, M.K., Scoffoni, C., Sack, L. 2012. The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecology Letters 15: 393–405.

\enc{Hölttä}{Holtta}, T., Cochard, H., Nikinmaa, E., Mencuccini, M. 2009. Capacitive effect of cavitation in xylem conduits: Results from a dynamic model. Plant, Cell and Environment 32: 10–21.

Martin-StPaul, N., Delzon, S., Cochard, H. 2017. Plant resistance to drought depends on timely stomatal closure. Ecology Letters 20: 1437–1447.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{hydraulics.psi2K}}, \code{\link{hydraulics.supplyFunctionPlot}}, \code{\link{swb}}, \code{\link{soil}}
}
\examples{
psi = seq(-10,0, by=0.1)
rwc_s = rep(NA, length(psi))
for(i in 1:length(psi)) rwc_s[i] = moisture.symplasticRWC(psi[i],-3,12)
plot(psi, rwc_s, type="l", xlab="Water potential (MPa)", ylab = "Symplasmic RWC")
}
