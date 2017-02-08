print.soil<-function(x,...) {
  #Depth
  cat(paste("Soil depth (mm):", round(x$SoilDepth, digits=0),"\n"))
  #Soil parameters related to texture
  nlayers = length(x$dVec)
  dini = 0;
  dfin = 0;
  for(l in 1:nlayers) {
    dfin = dfin+x$dVec[l]
    cat(paste("\nLayer ",l," [",dini," to ",dfin,"mm ]",
              "\n    clay (%):", round(x$clay[l]),"silt (%):", 100-round(x$sand[l])-round(x$clay[l]), "sand (%):", round(x$sand[l]),
              "[", x$usda_Type[l],"]\n    Rock fragment content (%):", round(x$rfc[l]),"Macroporosity (%):", round(x$macro[l]*100),  
              "\n    Theta FC (%):", round(100*x$Theta_FC[l]),"Vol. FC (mm):", round(x$Water_FC[l]), "Vol. current (mm):", round(x$W[l]*x$Water_FC[l]), "\n"))
    dini = dini+x$dVec[l]
  }
  cat(paste("\nTotal soil water holding capacity (mm):", round(sum(x$Water_FC), digits=0),"\n"))  
  cat(paste("\nTotal current Volume (mm):",round(sum(x$W*x$Water_FC), digits=0),"\n"))
}