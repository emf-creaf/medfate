\encoding{UTF-8}
\name{soil hydrology}
\alias{hydrology_infiltrationAmount}
\alias{hydrology_soilEvaporationAmount}
\alias{hydrology_soilEvaporation}
\alias{hydrology_infiltrationRepartition}
\alias{hydrology_snowMelt}

\title{Soil infiltration and bare soil evaporation}
\description{
Function \code{hydrology_infiltrationAmount} calculates the amount of water that infiltrates into the topsoil, according to the USDA SCS curve number method (Boughton 1989). The remaining is assumed to be lost as surface runoff. Function \code{hydrology_soilEvaporationAmount} calculates the amount of evaporation from bare soil, following Ritchie (1972). Function \code{hydrology_snowMelt} calculates the maximum amount of snowmelt according to Kergoat (1998).
}
\usage{
hydrology_infiltrationAmount(input, Ssoil) 
hydrology_infiltrationRepartition(I, dVec, macro, a = -0.005, b = 3)
hydrology_soilEvaporationAmount(DEF,PETs, Gsoil)
hydrology_soilEvaporation(soil, soilFunctions, pet, LgroundSWR, modifySoil = TRUE)
hydrology_snowMelt(tday, rad, LgroundSWR, elevation)
}
\arguments{
  \item{input}{A numeric vector of (daily) water input (in mm of water).}
  \item{Ssoil}{Soil water storage capacity (can be referred to topsoil) (in mm of water).}
  \item{DEF}{Water deficit in the (topsoil) layer.}
  \item{PETs}{Potential evapotranspiration at the soil surface.}
  \item{Gsoil}{Gamma parameter (maximum daily evaporation).}
  \item{I}{Soil infiltration (in mm of water).}
  \item{dVec}{Width of soil layers (in mm).}
  \item{macro}{Macroporosity of soil layers (in \%).}
  \item{soil}{An object of class \code{\link{soil}}.}
  \item{soilFunctions}{Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).}
  \item{pet}{Potential evapotranspiration for a given day (mm)}
  \item{LgroundSWR}{Proportion of short-wave radiation (SWR) reaching the ground.}
  \item{modifySoil}{Boolean flag to indicate that the input \code{soil} object should be modified during the simulation.}
  \item{tday}{Average day temperature (ºC).}
  \item{rad}{Solar radiation (in MJ/m2/day).}
  \item{elevation}{Altitude above sea level (m).}
  \item{a,b}{Parameters of the extinction function used for water infitration.}
}
\details{
See description of infiltration and soil evaporation processes in De Caceres et al. (2015).
}
\value{
Function \code{hydrology_infiltrationAmount} a vector of the same length as \code{input} containing the daily amount of water that infiltrates into the soil (in mm of water). Function \code{hydrology_infiltrationRepartition} estimates the amount of infiltrated water that reaches each soil layer. Function \code{hydrology_soilEvaporationAmount} returns the amount of water evaporated from the soil. Function \code{hydrology_soilEvaporation} returns a vector of water evaporated from each soil layer.
}
\references{
Boughton (1989). A review of the USDA SCS curve number method. - Australian Journal of Soil Research 27: 511-523.

De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to evaluate plant drought stress at the regional level. Agricultural and Forest Meteorology.

Kergoat L. (1998). A model for hydrological equilibrium of leaf area index on a global scale. Journal of Hydrology 212–213: 268–286.

Ritchie (1972). Model for predicting evaporation from a row crop with incomplete cover. - Water resources research.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{spwb}}, \code{\link{hydrology_soilWaterInputs}}
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

Theta_FC1=soil_psi2thetaSX(TS_clay, TS_sand, -33) #in m3/m3
Theta_FC2=soil_psi2thetaSX(SS_clay, SS_sand, -33) #in m3/m3
Theta_FC3=soil_psi2thetaSX(RL_clay, RL_sand, -33) #in m3/m3
pcTS_gravel = 1-(TS_gravel/100)
pcSS_gravel = 1-(SS_gravel/100)
pcRL_gravel = 1-(RL_gravel/100)
MaxVol1 = (d1*Theta_FC1*pcTS_gravel)
MaxVol2 = (d2*Theta_FC2*pcSS_gravel)
MaxVol3 = (d3*Theta_FC3*pcRL_gravel)
V = MaxVol1+MaxVol2+MaxVol3

par(mar=c(5,5,1,1), mfrow=c(1,2))
NP = seq(0,60, by=1)
plot(NP,hydrology_infiltrationAmount(NP, V[1]), type="l", xlim=c(0,60), ylim=c(0,60), 
     ylab="Infiltration (mm)", xlab="Net rainfall (mm)", frame=FALSE)
lines(NP,hydrology_infiltrationAmount(NP, V[2]), lty=2)
lines(NP,hydrology_infiltrationAmount(NP, V[3]), lty=3)
lines(NP,hydrology_infiltrationAmount(NP, V[4]), lty=4)
lines(NP,hydrology_infiltrationAmount(NP, V[5]), lty=5)
legend("topleft", bty="n", lty=1:5, 
       legend=c(paste("d =", SoilDepth, "Vsoil =",round(V),"mm")))
plot(NP,NP-hydrology_infiltrationAmount(NP, V[1]), type="l", xlim=c(0,60), ylim=c(0,60), 
     ylab="Runoff (mm)", xlab="Net rainfall (mm)", frame=FALSE)
lines(NP,NP-hydrology_infiltrationAmount(NP, V[2]), lty=2)
lines(NP,NP-hydrology_infiltrationAmount(NP, V[3]), lty=3)
lines(NP,NP-hydrology_infiltrationAmount(NP, V[4]), lty=4)
lines(NP,NP-hydrology_infiltrationAmount(NP, V[5]), lty=5)
legend("topleft", bty="n", lty=1:5, 
       legend=c(paste("d =", SoilDepth,"Vsoil =",round(V),"mm")))


}