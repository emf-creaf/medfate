#' @rdname soil_texture
#' 
#' @param layer Soil layer to be plotted.
#' @param relative Boolean flag to indicate that retention curve should be relative to field capacity or saturation.
#' @param to Either 'SAT' (saturation) or 'FC' (field capacity).
#'
soil_retentionCurvePlot<-function(soil, model="SX", layer = 1, 
                                  psi = seq(0, -6.0, by=-0.01),
                                  relative = TRUE, to = "SAT") {
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
    ylab=paste0("Moisture content (",toString,")")
  } else {
    ylab="Moisture content (% volume)"
  }
  df = data.frame(SX = y_sx, VG = y_vg, psi = -psi)
  xlab = "Soil water potential (-MPa)"
  g<-ggplot(df, aes(x=.data$psi))  
  if(model=="SX") {
    g<-g+ geom_path(aes(y=.data$SX), col="black")+
       ylab(ylab)+ xlab(xlab)+
       theme_bw()
    return(g)
  }
  else if(model=="VG") {
    g<-g+ geom_path(aes(y=.data$VG), col="black")+
      ylab(ylab)+ xlab(xlab)+
      theme_bw()
    return(g)  
  } 
  else {
    SX = df$SX
    VG = df$VG
    g<-g+ geom_path(aes(y=SX, linetype ="Saxton"), col="black")+
      geom_path(aes(y=VG, linetype ="Van Genuchten"), col="black")+
      scale_linetype_manual(name="", values=c("solid", "dashed"))+
      ylab(ylab)+ xlab(xlab)+
      theme_bw()
    return(g)  
  }
}