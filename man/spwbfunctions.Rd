\encoding{UTF-8}
\name{spwbprocesses}
\alias{spwb.SoilInfiltration}
\alias{spwb.SoilEvaporation}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil infiltration and bare soil evaporation}
\description{
Function \code{spwb.SoilInfiltration} calculates the amount of water that infiltrates into the topsoil, according to the USDA SCS curve number method (Boughton 1989). The remaining is assumed to be lost as surface runoff. Function \code{spwb.SoilEvaporation} calculates the amount of evaporation from bare soil, following Ritchie (1972).
}
\usage{
spwb.SoilInfiltration(input, Ssoil) 
spwb.SoilEvaporation(DEF,PETs, Gsoil)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{A numeric vector of (daily) water input (in mm of water).}
  \item{Ssoil}{Soil water storage capacity (can be referred to topsoil) (in mm of water).}
  \item{DEF}{Water deficit in the (topsoil) layer.}
  \item{PETs}{Potential evapotranspiration at the soil surface.}
  \item{Gsoil}{Gamma parameter (maximum daily evaporation).}
}
\details{
See description of infiltration and soil evaporation processes in De Caceres et al. (submitted).
}
\value{
Function \code{spwb.SoilInfiltration} a vector of the same length as \code{input} containing the daily amount of water that infiltrates into the soil (in mm of water). Function \code{spwb.SoilEvaporation} returns the amount of water evaporated from the soil.
}
\references{
Boughton (1989). A review of the USDA SCS curve number method. - Australian Journal of Soil Research 27: 511-523.

De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (submitted) Coupling a water balance model with forest inventory data to evaluate plant drought stress at the regional level. Agricultural and Forest Meteorology.

Ritchie (1972). Model for predicting evaporation from a row crop with incomplete cover. - Water resources research.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{spwb}}
}
\examples{
SoilDepth = c(200,400,800,1200,1500)

#TOPSOIL LAYERS
d1 = pmin(SoilDepth, 300) #<300
#SUBSOIL LAYERS
d2 = pmax(0, pmin(SoilDepth-300,1200)) #300-1500 mm
#ROCK LAYER
d3 = 4000-(d1+d2) #From SoilDepth down to 4.0 m

TS_clay = 15
TS_sand = 25
SS_clay = 15
SS_sand = 25
RL_clay = 15
RL_sand = 25
TS_gravel = 20
SS_gravel = 40
RL_gravel = 95

Theta_FC1=soil.psi2theta(TS_clay, TS_sand, -33) #in m3/m3
Theta_FC2=soil.psi2theta(SS_clay, SS_sand, -33) #in m3/m3
Theta_FC3=soil.psi2theta(RL_clay, RL_sand, -33) #in m3/m3
pcTS_gravel = 1-(TS_gravel/100)
pcSS_gravel = 1-(SS_gravel/100)
pcRL_gravel = 1-(RL_gravel/100)
MaxVol1 = (d1*Theta_FC1*pcTS_gravel)
MaxVol2 = (d2*Theta_FC2*pcSS_gravel)
MaxVol3 = (d3*Theta_FC3*pcRL_gravel)
V = MaxVol1+MaxVol2+MaxVol3

par(mar=c(5,5,1,1), mfrow=c(1,2))
NP = seq(0,60, by=1)
plot(NP,spwb.SoilInfiltration(NP, V[1]), type="l", xlim=c(0,60), ylim=c(0,60), 
     ylab="Infiltration (mm)", xlab="Net rainfall (mm)", frame=FALSE)
lines(NP,spwb.SoilInfiltration(NP, V[2]), lty=2)
lines(NP,spwb.SoilInfiltration(NP, V[3]), lty=3)
lines(NP,spwb.SoilInfiltration(NP, V[4]), lty=4)
lines(NP,spwb.SoilInfiltration(NP, V[5]), lty=5)
legend("topleft", bty="n", lty=1:5, legend=c(paste("d =", SoilDepth, "Vsoil =",round(V),"mm")))
plot(NP,NP-spwb.SoilInfiltration(NP, V[1]), type="l", xlim=c(0,60), ylim=c(0,60), 
     ylab="Runoff (mm)", xlab="Net rainfall (mm)", frame=FALSE)
lines(NP,NP-spwb.SoilInfiltration(NP, V[2]), lty=2)
lines(NP,NP-spwb.SoilInfiltration(NP, V[3]), lty=3)
lines(NP,NP-spwb.SoilInfiltration(NP, V[4]), lty=4)
lines(NP,NP-spwb.SoilInfiltration(NP, V[5]), lty=5)
legend("topleft", bty="n", lty=1:5, legend=c(paste("d =", SoilDepth,"Vsoil =",round(V),"mm")))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.