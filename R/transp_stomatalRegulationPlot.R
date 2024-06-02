#' @rdname transp_stomatalregulation
#' 
#' @param x An object of class \code{\link{spwbInput}} built using the 'Sperry' transpiration mode.
#' @param meteo A data frame with daily meteorological data series (see \code{\link{spwb}}).
#' @param day An integer to identify a day (a row) within \code{meteo}.
#' @param timestep An integer between 1 and \code{ndailysteps} specified in \code{x} (see \code{\link{defaultControl}}).
#' @param latitude Latitude (in degrees).
#' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North).
#' @param type A string with plot type, either \code{"E"} (transpiration flow), \code{"Ag"} (gross photosynthesis), \code{"An"} (net photosynthesis), \code{"Gsw"} (stomatal conductance to water vapour), \code{"T"} (temperature) or \code{"VPD"} (leaf vapour pressure deficit).
#' 
#' @keywords internal
transp_stomatalRegulationPlot<-function(x, meteo, day, timestep, latitude, elevation, slope = NA, aspect = NA,
                                        type = "E") {
  type = match.arg(type, c("E", "Ag","An" , "Gsw", "T", "VPD"))
  dctr = transp_transpirationSperry(x, meteo, day, latitude, elevation, slope, aspect,
                                    stepFunctions = timestep, 
                                    modifyInput = FALSE)
  ncoh = length(dctr$SupplyFunctions)

  l = dctr$SupplyFunctions
  cohnames = names(l)
  phsunlit = dctr$PhotoSunlitFunctions
  phshade = dctr$PhotoShadeFunctions
  pmsunlit = dctr$PMSunlitFunctions
  pmshade = dctr$PMShadeFunctions
  
  psi = numeric(0)
  E = numeric(0)
  Ag_sunlit = numeric(0)
  Ag_shade = numeric(0)
  An_sunlit = numeric(0)
  An_shade = numeric(0)
  Gsw_sunlit = numeric(0)
  Gsw_shade = numeric(0)
  Temp_sunlit = numeric(0)
  Temp_shade = numeric(0)
  VPD_sunlit = numeric(0)
  VPD_shade = numeric(0)
  cohorts = character(0)
  PM_sunlit = logical(0)
  PM_shade = logical(0)
  for(i in 1:ncoh) {
    psi = c(psi, -l[[i]]$psiLeaf)
    E = c(E, l[[i]]$E)
    Ag_sunlit = c(Ag_sunlit, phsunlit[[i]]$GrossPhotosynthesis)
    Ag_shade = c(Ag_shade, phshade[[i]]$GrossPhotosynthesis)
    An_sunlit = c(An_sunlit, phsunlit[[i]]$NetPhotosynthesis)
    An_shade = c(An_shade, phshade[[i]]$NetPhotosynthesis)
    Gsw_sunlit = c(Gsw_sunlit, phsunlit[[i]]$Gsw)
    Gsw_shade = c(Gsw_shade, phshade[[i]]$Gsw)
    Temp_sunlit = c(Temp_sunlit, phsunlit[[i]]$LeafTemperature)
    Temp_shade = c(Temp_shade, phshade[[i]]$LeafTemperature)
    VPD_sunlit = c(VPD_sunlit, phsunlit[[i]]$LeafVPD)
    VPD_shade = c(VPD_shade, phshade[[i]]$LeafVPD)
    PMsli = rep(F, length(l[[i]]$psiLeaf))
    PMshi = rep(F, length(l[[i]]$psiLeaf))
    PMsli[pmsunlit[[i]]$iMaxProfit+1] = T
    PMshi[pmshade[[i]]$iMaxProfit+1] = T
    PM_sunlit = c(PM_sunlit, PMsli)
    PM_shade = c(PM_shade, PMshi)
    cohorts = c(cohorts, rep(cohnames[i], length(l[[i]]$psiLeaf)))
  }
  df_sunlit = data.frame(psi = psi, E = E, 
                         Ag = Ag_sunlit, An = An_sunlit, Gsw = Gsw_sunlit,
                         Temp = Temp_sunlit, VPD = VPD_sunlit, 
                         PM = PM_sunlit,
                         cohort = cohorts, leaf = "sunlit", 
                         stringsAsFactors = F)
  df_shade = data.frame(psi = psi, E = E, 
                        Ag = Ag_shade, An = An_shade, Gsw = Gsw_shade,
                        Temp = Temp_shade, VPD = VPD_shade, 
                        PM = PM_shade,
                        cohort = cohorts, leaf = "shade",
                        stringsAsFactors = F)
  df = rbind(df_sunlit, df_shade)
  df$leaf = factor(df$leaf, levels = c("sunlit", "shade"))
  df_PM = df[df$PM,]
  g<-ggplot(df, aes(x=.data$psi))+
    xlab("Leaf pressure (-MPa)")+
    facet_wrap(~leaf)+
    theme_bw()
  if(type=="E") {
    g<- g + geom_path(aes(y = .data$E, col = .data$cohort))+
      geom_point(data = df_PM, aes(y=.data$E, col = .data$cohort))+
      ylab(expression(paste("Flow rate "(mmol%.%s^{-1}%.%m^{-2}))))
  } 
  else if(type=="Ag") {
    g<- g + geom_path(aes(y = .data$Ag, col = .data$cohort))+
      geom_point(data = df_PM, aes(y=.data$Ag, col = .data$cohort))+
      ylab(expression(paste("Gross photosynthesis "(mu*mol*C%.%s^{-1}%.%m^{-2}))))
  } 
  else if(type=="An") {
    g<- g + geom_path(aes(y = .data$An, col = .data$cohort))+
      geom_point(data = df_PM, aes(y=.data$An, col = .data$cohort))+
      ylab(expression(paste("Net photosynthesis "(mu*mol*C%.%s^{-1}%.%m^{-2}))))
  } 
  else if(type=="Gsw") {
    g<- g + geom_path(aes(y = .data$Gsw, col = .data$cohort))+
      geom_point(data = df_PM, aes(y=.data$Gsw, col = .data$cohort))+
      ylab(expression(paste("Stomatal conductance "(mol%.%s^{-1}%.%m^{-2}))))
  } 
  else if(type=="T") {
    g<- g + geom_path(aes(y = .data$Temp, col = .data$cohort))+
      geom_point(data = df_PM, aes(y=.data$Temp, col = .data$cohort))+
      ylab("Temperature (degrees C)")
  } 
  else if(type=="VPD") {
    g<- g + geom_path(aes(y = .data$VPD, col = .data$cohort))+
      geom_point(data = df_PM, aes(y=.data$VPD, col = .data$cohort))+
      ylab("Vapour pressure deficit (kPa)")
  } 
  g <- g + scale_color_discrete(name="")
  return(g)
}