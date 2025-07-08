#' Soil water retention and conductivity plots
#' 
#' Functions to display water retention curves and conductivity curves.
#' 
#' @param layer Soil layer to be plotted.
#' @param relative Boolean flag to indicate that retention curve should be relative to field capacity or saturation.
#' @param to Either 'SAT' (saturation) or 'FC' (field capacity).
#' @param soil Initialized soil object (returned by function \code{\link{soil}}).
#' @param model model Either 'SX' or 'VG' for Saxton's or Van Genuchten's water retention models; or 'both' to plot both retention models.
#' @param psi A numeric vector specifying a sequence of water potential values.
#' 
#' @details
#' \itemize{
#'    \item{\code{soil_retentionCurvePlot()} allows plotting the water retention curve of a given soil layer.}
#'    \item{\code{soil_conductivityCurvePlot()} allows plotting the conductivity curve of a given soil layer.}
#' }
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @return An object of class ggplot.
#' 
#' @seealso  \code{\link{soil_texture}}
#' @name soil_retentionCurvePlot
soil_retentionCurvePlot<-function(soil, model="SX", layer = 1, 
                                  psi = seq(0, -6.0, by=-0.01),
                                  relative = TRUE, to = "SAT") {
  if(!inherits(soil, "soil")) stop("Parameter 'soil' should be an object of class 'soil'")
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

#' @rdname soil_retentionCurvePlot
#' @param mmol Boolean flag to indicate that saturated conductivity units should be returned in mmol/m/s/MPa. If \code{mmol = FALSE} then units are cm/day.
#' @param log Boolean to display the y-axis in logarithm units
soil_conductivityCurvePlot<-function(soil, model="SX", layer = 1, 
                                psi = seq(0, -6.0, by=-0.01),
                                relative = TRUE, to = "SAT", log = TRUE, mmol = TRUE) {
  if(!inherits(soil, "soil")) stop("Parameter 'soil' should be an object of class 'soil'")
  model = match.arg(model, c("SX", "VG", "both"))
  to = match.arg(to, c("SAT", "FC"))
  theta_sx = unlist(lapply(as.list(psi), FUN=soil_psi2thetaSX, 
                       clay=soil$clay[layer], 
                       sand=soil$sand[layer], om = soil$om[layer]))
  y_sx = unlist(lapply(as.list(theta_sx), FUN=soil_unsaturatedConductivitySX, 
                       clay=soil$clay[layer], 
                       sand=soil$sand[layer], bd = soil$bd[layer], om = soil$om[layer], mmol = TRUE))
  y_vg = unlist(lapply(as.list(psi), FUN=soil_psi2kVG, 
                       ksat= soil$Ksat[layer],
                       alpha=soil$VG_alpha[layer], n = soil$VG_n[layer], 
                       theta_res = soil$VG_theta_res[layer], 
                       theta_sat = soil$VG_theta_sat[layer])) 
  if(relative) {
    if(to=="SAT") {
      toString = "% saturation" 
      cfc_sx = soil$Ksat[layer]
      cfc_vg = soil$Ksat[layer]
    } else if(to=="FC") {
      toString = "% field capacity" 
      theta_fc_sx = soil_thetaFC(soil, "SX")
      cfc_sx = soil_unsaturatedConductivitySX(theta_fc_sx[layer],clay=soil$clay[layer], 
                                              sand=soil$sand[layer], bd = soil$bd[layer], om = soil$om[layer], mmol = TRUE)
      cfc_vg = soil_psi2kVG(ksat= soil$Ksat[layer],
                                 alpha=soil$VG_alpha[layer], n = soil$VG_n[layer], 
                                 theta_res = soil$VG_theta_res[layer], 
                                 theta_sat = soil$VG_theta_sat[layer],
                                 psi = -0.033)
    }
    y_sx = 100*(y_sx/cfc_sx)
    y_vg = 100*(y_vg/cfc_vg)
    ylab=paste0("Conductivity (",toString,")")
  } else {
    if(mmol) {
      ylab=expression(paste("Conductivity ",(mmol%.%m^{-1}%.%MPa^{-1}%.%s^{-1})))
    } else {
      cmdTOmmolm2sMPa = 655.2934
      y_sx <- y_sx/cmdTOmmolm2sMPa
      y_vg <- y_vg/cmdTOmmolm2sMPa
      ylab=expression(paste("Conductivity ",(cm%.%d^{-1})))
    }
  }

  df = data.frame(SX = y_sx, VG = y_vg, psi = -psi)
  xlab = "Soil water potential (-MPa)"
  g<-ggplot(df, aes(x=.data$psi))  
  if(model=="SX") {
    g<-g+ geom_path(aes(y=.data$SX), col="black")
  }
  else if(model=="VG") {
    g<-g+ geom_path(aes(y=.data$VG), col="black")
  } 
  else {
    SX = df$SX
    VG = df$VG
    g<-g+ geom_path(aes(y=SX, linetype ="Saxton"), col="black")+
      geom_path(aes(y=VG, linetype ="Van Genuchten"), col="black")+
      scale_linetype_manual(name="", values=c("solid", "dashed"))
  }
  if(log) g <- g+
    scale_y_log10()
  g <- g+
    ylab(ylab)+ xlab(xlab)+
    theme_bw()
  return(g)  
}