print.soil<-function(x, model="SX",...) {
  #Depth
  cat(paste("Soil depth (mm):", round(x$SoilDepth, digits=0),"\n"))
  #Soil parameters related to texture
  nlayers = length(x$dVec)
  dini = 0;
  dfin = 0;
  ##Water content at field capacity
  Water_WP = soil_waterWP(x, model)
  Theta_WP = soil_thetaWP(x, model)
  Water_FC = soil_waterFC(x, model)
  Theta_FC = soil_thetaFC(x, model)
  Water_SAT = soil_waterSAT(x, model)
  Theta_SAT = soil_thetaSAT(x, model)
  Water_EXTR = soil_waterExtractable(x,model)
  
  for(l in 1:nlayers) {
    dfin = dfin+x$dVec[l]
    silt = 100-x$sand[l]-x$clay[l]
    if(!is.na(x$om[l])) silt = silt - x$om[l]
    cat(paste("\nLayer ",l," [",dini," to ",dfin,"mm ]",
              "\n    clay (%):", round(x$clay[l]),"silt (%):", round(silt), "sand (%):", round(x$sand[l]), "organic matter (%):", round(x$om[l]),
              "[", x$usda_Type[l],"]\n    Rock fragment content (%):", round(x$rfc[l]),"Macroporosity (%):", round(x$macro[l]*100),  
              "\n    Theta WP (%):", round(100*Theta_WP[l]),"Theta FC (%):", round(100*Theta_FC[l]), "Theta SAT (%):", round(100*Theta_SAT[l]), "Theta current (%)", round(100*x$W[l]*Theta_FC[l]),
              "\n    Vol. WP (mm):", round(Water_WP[l]),"Vol. FC (mm):", round(Water_FC[l]),"Vol. SAT (mm):", round(Water_SAT[l]), "Vol. current (mm):", round(x$W[l]*Water_FC[l]), 
              "\n    Temperature (Celsius):", round(x$Temp[l],1),
              "\n"))
    dini = dini+x$dVec[l]
  }
  cat(paste("\nTotal soil saturated capacity (mm):", round(sum(Water_SAT), digits=0),"\n"))  
  cat(paste("Total soil water holding capacity (mm):", round(sum(Water_FC), digits=0),"\n"))  
  cat(paste("Total soil extractable water (mm):", round(sum(Water_EXTR), digits=0),"\n"))  
  cat(paste("Total soil current Volume (mm):",round(sum(x$W*Water_FC), digits=0),"\n"))
  cat(paste("Water table depth (mm):",round(soil_waterTableDepth(x, model), digits=0),"\n"))
  cat(paste("\nSnow pack water equivalent (mm):",round(x$SWE, digits=0),"\n"))
}