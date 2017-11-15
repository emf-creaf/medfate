print.soil<-function(x,...) {
  #Depth
  cat(paste("Soil depth (mm):", round(x$SoilDepth, digits=0),"\n"))
  #Soil parameters related to texture
  nlayers = length(x$dVec)
  dini = 0;
  dfin = 0;
  ##Water content at field capacity
  Water_FC = soil.waterFC(x)
  
  for(l in 1:nlayers) {
    dfin = dfin+x$dVec[l]
    silt = 100-x$sand[l]-x$clay[l]
    if(!is.na(x$om[l])) silt = silt - x$om[l]
    cat(paste("\nLayer ",l," [",dini," to ",dfin,"mm ]",
              "\n    clay (%):", round(x$clay[l]),"silt (%):", round(silt), "sand (%):", round(x$sand[l]), "organic matter (%):", round(x$om[l]),
              "[", x$usda_Type[l],"]\n    Rock fragment content (%):", round(x$rfc[l]),"Macroporosity (%):", round(x$macro[l]*100),  
              "\n    Theta FC (%):", round(100*x$Theta_FC[l]),"Vol. FC (mm):", round(Water_FC[l]), "Vol. current (mm):", round(x$W[l]*Water_FC[l]), 
              "\n    Temperature (ºC):", round(x$Temp[l],1),
              "\n"))
    dini = dini+x$dVec[l]
  }
  cat(paste("\nTotal soil water holding capacity (mm):", round(sum(Water_FC), digits=0),"\n"))  
  cat(paste("\nTotal current Volume (mm):",round(sum(x$W*Water_FC), digits=0),"\n"))
}