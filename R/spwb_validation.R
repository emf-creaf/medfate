spwb_validation<-function(x, measuredData, type="SWC", cohort = NULL, draw = TRUE) {
  type = match.arg(type, c("SWC", "E", "ETR", "WP"))
  evalstats<-function(obs, pred) {
    E <- pred-obs
    Bias <- mean(E, na.rm=T)
    MAE <- mean(abs(E), na.rm=T)
    R2<- cor(obs, pred, use="complete")^2
    return(c(Bias= Bias, MAE = MAE, R2 = R2))
  }
  if(type=="SWC") {
    sm = x$Soil
    d = rownames(sm)
    seld = rownames(measuredData) %in% d
    fc = soil_thetaFC(x$soilInput, model = x$spwbInput$control$soilFunctions)
    q_mod = quantile(sm$W.1, p=c(0.05,0.95), na.rm=T)
    sm_mod = (sm$W.1-q_mod[1])/(q_mod[2]-q_mod[1])
    q_obs = quantile(measuredData$SWC[seld], p=c(0.05,0.95), na.rm=T)
    sm_obs = (measuredData$SWC[seld]-q_obs[1])/(q_obs[2]-q_obs[1])
    if(draw) {
      par(mfrow=c(2,1), mar=c(4,4,1,1))
      plot(as.Date(rownames(sm)), sm$W.1*fc[1], type="l", xlab="", ylab="Soil moisture 0-30 cm (% vol)", ylim=c(0,1))
      lines(as.Date(rownames(sm)), measuredData$SWC[seld], col="red")
      legend("topright", legend =c("Measured", "Modelled"), lty=1, col=c("red", "black"), bty="n")
      plot(as.Date(rownames(sm)), sm_mod, type="l", xlab="", ylab="Soil moisture 0-30 cm (scaled)", ylim=c(0,1))
      lines(as.Date(rownames(sm)), sm_obs, col="red")
    }
    return(evalstats(sm_obs, sm_mod))
  } 
  else if(type=="ETR") {
    ET1 = x$WaterBalance$SoilEvaporation+x$WaterBalance$Transpiration
    ET2 = x$WaterBalance$Evapotranspiration
    d = rownames(x$WaterBalance)
    seld = rownames(measuredData) %in% d
    
    if(draw) {
      par(mfrow=c(1,2), mar=c(4,4,1,1))
      ETmax = ceiling(max(c(ET1, ET2, measuredData$ETR[seld]), na.rm=T))
      plot(as.Date(d), ET2, type="l", ylim=c(0,ETmax), col="gray",
           xlab = "", ylab="ETR (mm)")
      lines(as.Date(d), ET1, col="black")
      lines(as.Date(row.names(measuredData))[seld], measuredData$ETR[seld], col="red")
      legend("topright", legend = c("modelled Es+Tr", "modelled Es+Tr+In", "measured ETR"), col=c("black", "gray","red"), lty=1, bty="n")
      plot(ET2, measuredData$ETR[seld], xlab="modelled ETR (mm)", ylab="measured ETR (mm)",
           asp=1, xlim=c(0,ETmax), ylim=c(0,ETmax), pch=19, cex=0.4, col="gray")
      points(ET1, measuredData$ETR[seld], col="black", pch=19, cex=0.4)
      abline(a=0, b=1, col="black")
    }
    df = as.data.frame(rbind(evalstats(measuredData$ETR[seld], ET1),
                       evalstats(measuredData$ETR[seld], ET2)))
    row.names(df)<-c("Es+Tr", "Es+Tr+In")
    return(df)
    
  }
  else if(type=="E") {
    pt = x$PlantTranspiration
    d = rownames(pt)
    LAI = x$spwbInput$above$LAI_live
    spnames = x$spwbInput$cohorts$Name
    allcohnames = row.names(x$spwbInput$cohorts)
    seld = rownames(measuredData) %in% d
    
    if(!is.null(cohort)) {
      obscolumn = paste0("E_", cohort)
      icoh = which(allcohnames==cohort)
      
      E_mod = pt[,icoh]/LAI[icoh]
      E_obs = measuredData[[obscolumn]][seld]/LAI[icoh]
      Emax = max(c(E_mod, E_obs), na.rm=T)
      if(draw) {
        par(mfrow=c(1,2), mar=c(4,4,1,1))
        plot(as.Date(rownames(pt)), E_mod, type="l", ylab="Transpiration per leaf area (mm/day)", xlab="",
             ylim=c(0,Emax))
        lines(as.Date(row.names(measuredData))[seld], E_obs, col="red")
        legend("topright", legend =c("Measured", "Modelled"), lty=1, col=c("red", "black"), bty="n")
        plot(E_mod, E_obs, cex=0.5, pch = 19,
             xlab ="Modelled transpiration", ylab="Measured transpiration", 
             asp=1, xlim=c(0,Emax), ylim=c(0,Emax))
        abline(a=0, b=1)
      }
      return(evalstats(E_obs, E_mod))
    } else {
      obscolumns = paste0("E_", allcohnames)
      selc = obscolumns %in% names(measuredData)
      cohorts = allcohnames[selc]
      spnames = spnames[selc]
      obscolumns = obscolumns[selc]
      df = data.frame(Bias = rep(NA, length(cohorts)),
                      MAE = rep(NA, length(cohorts)),
                      R2 = rep(NA, length(cohorts)))
      row.names(df)<-cohorts
      if(draw)  par(mfrow=c(length(cohorts),2), mar=c(4,4,4,1))
      for(i in 1:length(cohorts)) {
        icoh = which(allcohnames==cohorts[i])
        obscolumn = obscolumns[i]
        E_mod = pt[,icoh]/LAI[icoh]
        E_obs = measuredData[[obscolumn]][seld]/LAI[icoh]
        Emax = max(c(E_mod, E_obs), na.rm=T)
        if(sum(!is.na(E_obs))>0) {
          if(draw) {
            plot(as.Date(rownames(pt)), E_mod, type="l", ylab="Transpiration per leaf area (mm/day)", xlab="",
                 main = paste0(cohorts[i], " (",spnames[i],")"), ylim=c(0,Emax))
            lines(as.Date(row.names(measuredData))[seld], E_obs, col="red")
            legend("topright", legend =c("Measured", "Modelled"), lty=1, col=c("red", "black"), bty="n")
            plot(E_mod, E_obs, cex=0.5, pch = 19,
                 xlab ="Modelled transpiration", ylab="Measured transpiration", 
                 asp=1, xlim=c(0,Emax), ylim=c(0,Emax))
            abline(a=0, b=1)
          }
          df[i,] = evalstats(E_obs, E_mod)
        } else {
          message(paste0("Not enough observations for ", cohorts[i]))
        }
      }
      return(df)
    }
  }
}