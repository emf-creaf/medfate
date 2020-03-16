transp_stomatalRegulationPlot<-function(x, soil, meteo, day, timestep, latitude, elevation, slope = NA, aspect = NA,
                                        type = "E") {
  type = match.arg(type, c("E","An" , "Gw", "T", "VPD"))
  dctr = transp_transpirationSperry(x, soil, meteo, day, latitude, elevation, slope, aspect,
                                    stepFunctions = timestep, 
                                    modifyInputX = FALSE, modifyInputSoil = FALSE)
  ncoh = length(dctr$SupplyFunctions)

  l = dctr$SupplyFunctions
  cohnames = names(l)
  phsunlit = dctr$PhotoSunlitFunctions
  phshade = dctr$PhotoShadeFunctions
  pmsunlit = dctr$PMSunlitFunctions
  pmshade = dctr$PMShadeFunctions
  
  psi = numeric(0)
  E = numeric(0)
  An_sunlit = numeric(0)
  An_shade = numeric(0)
  Gw_sunlit = numeric(0)
  Gw_shade = numeric(0)
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
    An_sunlit = c(An_sunlit, phsunlit[[i]]$Photosynthesis)
    An_shade = c(An_shade, phshade[[i]]$Photosynthesis)
    Gw_sunlit = c(Gw_sunlit, phsunlit[[i]]$WaterVaporConductance)
    Gw_shade = c(Gw_shade, phshade[[i]]$WaterVaporConductance)
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
                         An = An_sunlit, Gw = Gw_sunlit,
                         Temp = Temp_sunlit, VPD = VPD_sunlit, 
                         PM = PM_sunlit,
                         cohort = cohorts, leaf = "sunlit", 
                         stringsAsFactors = F)
  df_shade = data.frame(psi = psi, E = E, 
                        An = An_shade, Gw = Gw_shade,
                        Temp = Temp_shade, VPD = VPD_shade, 
                        PM = PM_shade,
                        cohort = cohorts, leaf = "shade",
                        stringsAsFactors = F)
  df = rbind(df_sunlit, df_shade)
  df$leaf = factor(df$leaf, levels = c("sunlit", "shade"))
  df_PM = df[df$PM,]
  g<-ggplot(df, aes_string(x="psi"))+
    xlab("Leaf pressure (-MPa)")+
    facet_wrap(~leaf)+
    theme_bw()
  if(type=="E") {
    g<- g + geom_path(aes_string(y = "E", col = "cohort"))+
      geom_point(data = df_PM, aes_string(y="E", col = "cohort"))+
      ylab(expression(paste("Flow rate "(mmol%.%s^{-1}%.%m^{-2}))))
  } 
  else if(type=="An") {
    g<- g + geom_path(aes_string(y = "An", col = "cohort"))+
      geom_point(data = df_PM, aes_string(y="An", col = "cohort"))+
      ylab(expression(paste("Photosynthesis "(mu*mol*C%.%s^{-1}%.%m^{-2}))))
  } 
  else if(type=="Gw") {
    g<- g + geom_path(aes_string(y = "Gw", col = "cohort"))+
      geom_point(data = df_PM, aes_string(y="Gw", col = "cohort"))+
      ylab(expression(paste("Stomatal conductance "(mol%.%s^{-1}%.%m^{-2}))))
  } 
  else if(type=="T") {
    g<- g + geom_path(aes_string(y = "Temp", col = "cohort"))+
      geom_point(data = df_PM, aes_string(y="Temp", col = "cohort"))+
      ylab("Temperature (degrees C)")
  } 
  else if(type=="VPD") {
    g<- g + geom_path(aes_string(y = "VPD", col = "cohort"))+
      geom_point(data = df_PM, aes_string(y="VPD", col = "cohort"))+
      ylab("Vapour pressure deficit (kPa)")
  } 
  g <- g + scale_color_discrete(name="")
  return(g)
}