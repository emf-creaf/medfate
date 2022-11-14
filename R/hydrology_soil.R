#' Soil infiltration, percolation and bare soil evaporation
#'
#' Function \code{hydrology_infiltrationAmount} calculates the amount of water that infiltrates 
#' into the topsoil, according to the USDA SCS curve number method (Boughton 1989). 
#' The remaining is assumed to be lost as surface runoff. 
#' Function \code{hydrology_soilEvaporationAmount} calculates the amount of evaporation from bare soil, following Ritchie (1972). 
#' Function \code{hydrology_snowMelt} calculates the maximum amount of snowmelt according to Kergoat (1998).
#' 
#' @param input A numeric vector of (daily) water input (in mm of water).
#' @param Ssoil Soil water storage capacity (can be referred to topsoil) (in mm of water).
#' 
#' @details See description of infiltration and soil evaporation processes in De Caceres et al. (2015).
#' 
#' @return 
#' Function \code{hydrology_infiltrationAmount} a vector of the same length as \code{input} containing the daily amount of water that infiltrates into the soil (in mm of water). 
#' 
#' Function \code{hydrology_infiltrationRepartition} estimates the amount of infiltrated water that reaches each soil layer. 
#' 
#' Function \code{hydrology_soilEvaporationAmount} returns the amount of water evaporated from the soil. 
#' 
#' Function \code{hydrology_soilEvaporation} returns a vector of water evaporated from each soil layer.
#' 
#' @references 
#' Boughton (1989). A review of the USDA SCS curve number method. - Australian Journal of Soil Research 27: 511-523.
#' 
#' De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to evaluate plant drought stress at the regional level. Agricultural and Forest Meteorology.
#' 
#' Kergoat L. (1998). A model for hydrological equilibrium of leaf area index on a global scale. Journal of Hydrology 212–213: 268–286.
#' 
#' Ritchie (1972). Model for predicting evaporation from a row crop with incomplete cover. - Water resources research.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso  \code{\link{spwb}}, \code{\link{hydrology_soilWaterInputs}}
#' 
#' @examples 
#' SoilDepth = c(200,400,800,1200,1500)
#' 
#' #TOPSOIL LAYERS
#' d1 = pmin(SoilDepth, 300) #<300
#' #SUBSOIL LAYERS
#' d2 = pmax(0, pmin(SoilDepth-300,1200)) #300-1500 mm
#' #ROCK LAYER
#' d3 = 4000-(d1+d2) #From SoilDepth down to 4.0 m
#' 
#' TS_clay = 15
#' TS_sand = 25
#' SS_clay = 15
#' SS_sand = 25
#' RL_clay = 15
#' RL_sand = 25
#' TS_gravel = 20
#' SS_gravel = 40
#' RL_gravel = 95
#' 
#' Theta_FC1=soil_psi2thetaSX(TS_clay, TS_sand, -33) #in m3/m3
#' Theta_FC2=soil_psi2thetaSX(SS_clay, SS_sand, -33) #in m3/m3
#' Theta_FC3=soil_psi2thetaSX(RL_clay, RL_sand, -33) #in m3/m3
#' pcTS_gravel = 1-(TS_gravel/100)
#' pcSS_gravel = 1-(SS_gravel/100)
#' pcRL_gravel = 1-(RL_gravel/100)
#' MaxVol1 = (d1*Theta_FC1*pcTS_gravel)
#' MaxVol2 = (d2*Theta_FC2*pcSS_gravel)
#' MaxVol3 = (d3*Theta_FC3*pcRL_gravel)
#' V = MaxVol1+MaxVol2+MaxVol3
#' 
#' par(mar=c(5,5,1,1), mfrow=c(1,2))
#' NP = seq(0,60, by=1)
#' plot(NP,hydrology_infiltrationAmount(NP, V[1]), type="l", xlim=c(0,60), ylim=c(0,60), 
#'      ylab="Infiltration (mm)", xlab="Net rainfall (mm)", frame=FALSE)
#' lines(NP,hydrology_infiltrationAmount(NP, V[2]), lty=2)
#' lines(NP,hydrology_infiltrationAmount(NP, V[3]), lty=3)
#' lines(NP,hydrology_infiltrationAmount(NP, V[4]), lty=4)
#' lines(NP,hydrology_infiltrationAmount(NP, V[5]), lty=5)
#' legend("topleft", bty="n", lty=1:5, 
#'        legend=c(paste("d =", SoilDepth, "Vsoil =",round(V),"mm")))
#' plot(NP,NP-hydrology_infiltrationAmount(NP, V[1]), type="l", xlim=c(0,60), ylim=c(0,60), 
#'      ylab="Runoff (mm)", xlab="Net rainfall (mm)", frame=FALSE)
#' lines(NP,NP-hydrology_infiltrationAmount(NP, V[2]), lty=2)
#' lines(NP,NP-hydrology_infiltrationAmount(NP, V[3]), lty=3)
#' lines(NP,NP-hydrology_infiltrationAmount(NP, V[4]), lty=4)
#' lines(NP,NP-hydrology_infiltrationAmount(NP, V[5]), lty=5)
#' legend("topleft", bty="n", lty=1:5, 
#'        legend=c(paste("d =", SoilDepth,"Vsoil =",round(V),"mm")))
#' 
#' @name hydrology_soil
hydrology_infiltrationAmount<-function(input, Ssoil){
  I = rep(0, length(input))
  sel = input>(0.2*Ssoil)
  I[sel] = input[sel]-((input[sel]-0.2*Ssoil)^2/(input[sel]+0.8*Ssoil))
  I[!sel] = input[!sel]
  return(I)
}

