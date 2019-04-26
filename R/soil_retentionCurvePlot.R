soil_retentionCurvePlot<-function(soil, model="SX", layer = 1, 
                                  psi = seq(0, -6.0, by=-0.01),
                                  relative = TRUE, to = "SAT", ...) {
  model = match.arg(model, c("SX", "VG"))
  to = match.arg(to, c("SAT", "FC"))
  if(model=="SX") {
    y = unlist(lapply(as.list(psi), FUN=soil_psi2thetaSX, 
                    clay=soil$clay[layer], 
                    sand=soil$sand[layer], om = soil$om[layer]))
  } else {
   y = unlist(lapply(as.list(psi), FUN=soil_psi2thetaVG, 
                     alpha=soil$VG_alpha[layer], n = soil$VG_n[layer], 
                     theta_res = soil$VG_theta_res[layer], 
                     theta_sat = soil$VG_theta_sat[layer])) 
  }
  if(relative) {
    if(to=="SAT") {
      tfc = soil_thetaSAT(soil, model)
      y = 100*(y/tfc[layer])
      plot(-psi, y, 
           type="l",ylab="Moisture content (% saturation)", 
           xlab = "Soil water potential (-MPa)", ...)
    } else if(to=="FC") {
      tfc = soil_thetaFC(soil, model)
      y = 100*(y/tfc[layer])
      plot(-psi, y, 
           type="l",ylab="Moisture content (% field capacity)", 
           xlab = "Soil water potential (-MPa)", ...)
    }

  } else {
    plot(-psi, y, 
         type="l",ylab="Moisture content (% volume)", 
         xlab = "Soil water potential (-MPa)", ...)
  }
}