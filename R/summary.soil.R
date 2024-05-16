#' @rdname soil
#' 
#' @param object An object of class \code{soil}.
#' @param model Either 'SX' or 'VG' for Saxton or Van Genuchten pedotransfer models.
#' @param ... Additional parameters to \code{summary}.
#' 
summary.soil<-function(object, model="SX",...) {
  #Depth
  cat(paste("Soil depth (mm):", round(sum(object$widths), digits=0),"\n"))
  #Soil parameters related to texture
  nlayers = length(object$widths)
  dini = 0;
  dfin = 0;
  ##Water content at field capacity
  Water_WP = soil_waterWP(object, model)
  Theta_WP = soil_thetaWP(object, model)
  Water_FC = soil_waterFC(object, model)
  Theta_FC = soil_thetaFC(object, model)
  Water_SAT = soil_waterSAT(object, model)
  Theta_SAT = soil_thetaSAT(object, model)
  Water_EXTR = soil_waterExtractable(object,model)
  
  for(l in 1:nlayers) {
    dfin = dfin+object$widths[l]
    silt = 100-object$sand[l]-object$clay[l]
    if(!is.na(object$om[l])) silt = silt - object$om[l]
    usda_Type = soil_USDAType(object$clay[l],object$sand[l]);
    
    cat(paste("\nLayer ",l," [",dini," to ",dfin,"mm ]",
              "\n    clay (%):", round(object$clay[l]),"silt (%):", round(silt), "sand (%):", round(object$sand[l]), "organic matter (%):", round(object$om[l]),
              "[",   usda_Type,"]\n    Rock fragment content (%):", round(object$rfc[l]),"Macroporosity (%):", round(object$macro[l]*100),  
              "\n    Theta WP (%):", round(100*Theta_WP[l]),"Theta FC (%):", round(100*Theta_FC[l]), "Theta SAT (%):", round(100*Theta_SAT[l]), "Theta current (%)", round(100*object$W[l]*Theta_FC[l]),
              "\n    Vol. WP (mm):", round(Water_WP[l]),"Vol. FC (mm):", round(Water_FC[l]),"Vol. SAT (mm):", round(Water_SAT[l]), "Vol. current (mm):", round(object$W[l]*Water_FC[l]), 
              "\n    Temperature (Celsius):", round(object$Temp[l],1),
              "\n"))
    dini = dini+object$widths[l]
  }
  cat(paste("\nTotal soil saturated capacity (mm):", round(sum(Water_SAT), digits=0),"\n"))  
  cat(paste("Total soil water holding capacity (mm):", round(sum(Water_FC), digits=0),"\n"))  
  cat(paste("Total soil extractable water (mm):", round(sum(Water_EXTR), digits=0),"\n"))  
  cat(paste("Total soil current Volume (mm):",round(sum(object$W*Water_FC), digits=0),"\n"))
  cat(paste("Saturated water depth (mm):",round(soil_saturatedWaterDepth(object, model), digits=0),"\n\n"))
}