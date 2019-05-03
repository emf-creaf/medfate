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
      plot(wtMD[,1], frapue$measuredData$MD_T1_68[seld])
      abline(a=0,b=1)
      
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
  else if(type=="WP"){
    wtMD = x$LeafPsiMin
    wtPD = x$LeafPsiMax
    d = rownames(wtMD)
    spnames = x$spwbInput$cohorts$Name
    allcohnames = row.names(x$spwbInput$cohorts)
    seld = rownames(measuredData) %in% d
    
    if(!is.null(cohort)) {
      pdcolumn = paste0("PD_", cohort)
      mdcolumn = paste0("MD_", cohort)
      pderrcolumn = paste0("PD_", cohort, "_err")
      mderrcolumn = paste0("MD_", cohort, "_err")
      icoh = which(allcohnames==cohort)
      
      PD_mod = wtPD[,icoh]
      MD_mod = wtMD[,icoh]
      PD_obs = measuredData[[pdcolumn]][seld]
      MD_obs = measuredData[[mdcolumn]][seld]
      PD_obs_err = measuredData[[pderrcolumn]][seld]
      MD_obs_err = measuredData[[mderrcolumn]][seld]
      wpmin = min(c(PD_mod, MD_mod, PD_obs, MD_obs), na.rm=T)
      if(draw) {
        par(mfrow=c(1,2), mar=c(4,4,4,1))
        plot(as.Date(rownames(wtMD)), MD_mod, type="l", col="red", xlab="", ylab = "Leaf water potential (MPa)",
             ylim = c(wpmin,0), main = paste0(cohort, " (", spnames[icoh],")"))
        lines(as.Date(rownames(wtPD)), PD_mod, col="blue")
        arrows(x0=as.Date(row.names(measuredData)[seld]), 
               y0 = PD_obs, 
               y1 = PD_obs+1.96*PD_obs_err, col="black", length = 0.01, angle=90)
        arrows(x0=as.Date(row.names(measuredData)[seld]), 
               y0 = PD_obs, 
               y1 = PD_obs-1.96*PD_obs_err, col="black", length = 0.01, angle=90)
        points(as.Date(row.names(measuredData)[seld]), PD_obs, col="blue", pch=19)
        arrows(x0=as.Date(row.names(measuredData)[seld]), 
               y0 = MD_obs, 
               y1 = MD_obs+1.96*MD_obs_err, col="black", length = 0.01, angle=90)
        arrows(x0=as.Date(row.names(measuredData)[seld]), 
               y0 = MD_obs, 
               y1 = MD_obs-1.96*MD_obs_err, col="black", length = 0.01, angle=90)
        points(as.Date(row.names(measuredData)[seld]), MD_obs, col="red", pch=19)
        plot(PD_mod, PD_obs, col="blue", cex = 0.5, 
             xlab ="Modelled leaf water potential (MPa)", ylab="Measured leaf water potential (MPa)",
             xlim = c(wpmin,0), ylim = c(wpmin,0), asp=1, pch=19)
        points(MD_mod, MD_obs, col="red", cex=0.5, pch=19)
        legend("bottomright", col=c("blue", "red"), pch=19, legend=c("Predawn", "Midday"), bty="n")
        abline(a=0,b=1)
      }
      df = as.data.frame(rbind(evalstats(PD_obs, PD_mod),
                               evalstats(MD_obs, MD_mod)))
      row.names(df)<-c("Predawn potentials", "Midday potentials")
      return(df)
    }
  }
}