soil_retentionCurvePlot<-function(soil, model="SX", layer = 1, 
                                  psi = seq(0, -6.0, by=-0.01),
                                  relative = TRUE, to = "SAT", ...) {
  model = match.arg(model, c("SX", "VG", "both"))
  to = match.arg(to, c("SAT", "FC"))
  y_sx = unlist(lapply(as.list(psi), FUN=soil_psi2thetaSX, 
                       clay=soil$clay[layer], 
                       sand=soil$sand[layer], om = soil$om[layer]))
  y_vg = unlist(lapply(as.list(psi), FUN=soil_psi2thetaVG, 
                       alpha=soil$VG_alpha[layer], n = soil$VG_n[layer], 
                       theta_res = soil$VG_theta_res[layer], 
                       theta_sat = soil$VG_theta_sat[layer])) 
  if(relative) {
    if(to=="SAT") {
      toString = "% saturation" 
      tfc_sx = soil_thetaSAT(soil, "SX")
      tfc_vg = soil_thetaSAT(soil, "VG")
    } else if(to=="FC") {
      toString = "% field capacity" 
      tfc_sx = soil_thetaFC(soil, "SX")
      tfc_vg = soil_thetaFC(soil, "VG")
    }
    y_sx = 100*(y_sx/tfc_sx[layer])
    y_vg = 100*(y_vg/tfc_vg[layer])
    if(model=="SX") {
      plot(-psi, y_sx, 
           type="l",ylab=paste0("Moisture content (",toString,")"), 
           xlab = "Soil water potential (-MPa)", ...)
    }
    else if(model=="VG") {
      plot(-psi, y_vg, 
           type="l",ylab=paste0("Moisture content (",toString,")"), 
           xlab = "Soil water potential (-MPa)", ...)
    }
    else {
      plot(-psi, y_sx, lty=1,
           type="l",ylab=paste0("Moisture content (",toString,")"), 
           xlab = "Soil water potential (-MPa)", ...)
      lines(-psi, y_vg, lty=2)
      legend("topright", legend=c("Saxton", "Van Genuchten"), lty=c(1,2), bty="n")
    }
    
  } else {
    if(model=="SX") {
      plot(-psi, y_sx, 
           type="l",ylab="Moisture content (% volume)", 
           xlab = "Soil water potential (-MPa)", ...)
    }
    else if(model=="VG") {
      plot(-psi, y_vg, 
           type="l",ylab="Moisture content (% volume)", 
           xlab = "Soil water potential (-MPa)", ...)
    }
    else {
      plot(-psi, y_sx, 
           type="l",ylab="Moisture content (% volume)", 
           xlab = "Soil water potential (-MPa)", ...)
      lines(-psi, y_vg, lty=2)
      legend("topright", legend=c("Saxton", "Van Genuchten"), lty=c(1,2), bty="n")
    }
  }
}