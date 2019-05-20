spwb_validation<-function(x, measuredData, type="SWC", cohort = NULL, draw = TRUE,
                          plotType = "dynamics") {
  scatterplot<-function(df, xlab="", ylab="") {
    g<-ggplot(df, aes(x=Modelled, y = Observed))+
      geom_point()+
      geom_abline(intercept=0, slope=1, col="black")+
      geom_smooth(method="lm", se = FALSE, col="gray", linetype="dashed")+
      xlab(xlab)+
      ylab(ylab)+
      theme_bw()
    return(g)
  }
  dynamicsplot<-function(df, xlab="", ylab="") {
    g<-ggplot(df, aes(x=Dates))+
      geom_path(aes(y=Observed, col="Observed"))+
      geom_path(aes(y=Modelled, col="Modelled"))+
      xlab(xlab)+
      ylab(ylab)+
      scale_color_manual(name="", values=c("Observed"="black", "Modelled"= "red"))+
      theme_bw()
    return(g)
  }
  evalstats<-function(obs, pred) {
    E <- pred-obs
    Bias <- mean(E, na.rm=T)
    MAE <- mean(abs(E), na.rm=T)
    R2<- cor(obs, pred, use="complete")^2
    return(c(n = sum(!is.na(obs) & !is.na(pred)), Bias= Bias, MAE = MAE, R2 = R2))
  }
  
  # Check arguments
  type = match.arg(type, c("SWC", "SWC_scaled","E", "ETR", "WP"))
  plotType = match.arg(plotType, c("dynamics", "scatter"))

  if(type=="SWC") {
    sm = x$Soil
    d = rownames(sm)
    fc = soil_thetaFC(x$soilInput, model = x$spwbInput$control$soilFunctions)
    df <- data.frame(Observed = NA, Modelled = sm$W.1*fc[1], Dates = as.Date(d))
    
    if(!("SWC" %in% names(measuredData))) stop(paste0("Column 'SWC' not found in measured data frame."))
    seld = rownames(measuredData) %in% d
    df$Observed[d %in% rownames(measuredData)] = measuredData$SWC[seld]
  
    if(draw) {
      if(plotType=="dynamics") {
        g<-dynamicsplot(df, ylab = "Soil moisture (% vol)")
      } else {
        g<-scatterplot(df, xlab  = "Modelled soil moisture (% vol)",
                       ylab = "Measured soil moisture (% vol)")
      }
      print(evalstats(df$Observed, df$Modelled))
      return(g)
    } else{
      print(evalstats(df$Observed, df$Modelled))
    }
  } 
  else if(type=="SWC_scaled") {
    sm = x$Soil
    d = rownames(sm)
    fc = soil_thetaFC(x$soilInput, model = x$spwbInput$control$soilFunctions)
    q_mod = quantile(sm$W.1, p=c(0.05,0.95), na.rm=T)
    df <- data.frame(Observed = NA, Modelled = (sm$W.1-q_mod[1])/(q_mod[2]-q_mod[1]), Dates = as.Date(d))
    
    if(!("SWC" %in% names(measuredData))) stop(paste0("Column 'SWC' not found in measured data frame."))
    seld = rownames(measuredData) %in% d
    q_obs = quantile(measuredData$SWC[seld], p=c(0.05,0.95), na.rm=T)
    df$Observed[d %in% rownames(measuredData)] =(measuredData$SWC[seld]-q_obs[1])/(q_obs[2]-q_obs[1])
    
    if(draw) {
      if(plotType=="dynamics") {
        g<-dynamicsplot(df, ylab = "Soil moisture (scaled)")
      } else {
        g<-scatterplot(df, xlab  = "Modelled soil moisture (scaled)",
                       ylab = "Measured soil moisture (scaled)")

      }
      print(evalstats(df$Observed, df$Modelled))
      return(g)
    } else{
      print(evalstats(df$Observed, df$Modelled))
    }
  } 
  else if(type=="E") {
    pt = x$PlantTranspiration
    d = rownames(pt)
    LAI = x$spwbInput$above$LAI_live
    spnames = x$spwbInput$cohorts$Name
    allcohnames = row.names(x$spwbInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Observed = NA, Modelled = pt[,icoh]/LAI[icoh], Dates = as.Date(d))
    ## Fill observed values
    obscolumn = paste0("E_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
    
    if(draw) {
      if(plotType=="dynamics") {
        g<-dynamicsplot(df, ylab = "Transpiration per leaf area (mm/day)")
      } else {
        g<-scatterplot(df, 
                       xlab = "Modelled transpiration per leaf area (mm/day)",
                       ylab = "Measured transpiration per leaf area (mm/day)")
      }
      print(evalstats(df$Observed, df$Modelled))
      return(g)
    } else{
      return(evalstats(df$Observed, df$Modelled))
    }
  }
  else if(type=="ETR") {
    ET1 = x$WaterBalance$SoilEvaporation+x$WaterBalance$Transpiration
    ET2 = x$WaterBalance$Evapotranspiration
    d = rownames(x$WaterBalance)
    seld = rownames(measuredData) %in% d
    ETobs = measuredData$ETR[seld]
    
    if(draw) {
      par(mfrow=c(1,2), mar=c(4,4,1,1))
      ETmax = ceiling(max(c(ET1, ET2, ETobs), na.rm=T))
      plot(as.Date(d), ET2, type="l", ylim=c(0,ETmax), col="gray",
           xlab = "", ylab="ETR (mm)")
      lines(as.Date(d), ET1, col="black")
      lines(as.Date(row.names(measuredData))[seld], ETobs, col="red")
      legend("topright", legend = c("modelled Es+Tr", "modelled Es+Tr+In", "measured ETR"), col=c("black", "gray","red"), lty=1, bty="n")
      plot(ET2, ETobs, xlab="modelled ETR (mm)", ylab="measured ETR (mm)",
           asp=1, xlim=c(0,ETmax), ylim=c(0,ETmax), pch=19, cex=0.4, col="gray")
      points(ET1, ETobs, col="black", pch=19, cex=0.4)
      abline(a=0, b=1, col="black")
      abline(lm(ETobs~ ET1), col="black", lty=2)
      abline(lm(ETobs~ ET2), col="gray", lty=2)
    }
    df = as.data.frame(rbind(evalstats(ETobs, ET1),
                       evalstats(ETobs, ET2)))
    row.names(df)<-c("Es+Tr", "Es+Tr+In")
    return(df)
    
  }
  else if(type=="E") {
    pt = x$PlantTranspiration
    d = rownames(pt)
    LAI = x$spwbInput$above$LAI_live
    spnames = x$spwbInput$cohorts$Name
    allcohnames = row.names(x$spwbInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Observed = NA, Modelled = pt[,icoh]/LAI[icoh], Dates = as.Date(d))
    ## Fill observed values
    obscolumn = paste0("E_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
    
    if(draw) {
      if(plotType=="dynamics") {
        g<-dynamicsplot(df, ylab = "Transpiration per leaf area (mm/day)")
      } else {
        g<-scatterplot(df, 
                       xlab = "Modelled transpiration per leaf area (mm/day)",
                       ylab = "Measured transpiration per leaf area (mm/day)")
      }
      print(evalstats(df$Observed, df$Modelled))
      return(g)
    } else{
      return(evalstats(df$Observed, df$Modelled))
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
        abline(lm(PD_obs~PD_mod), col="blue", lty=2)
        points(MD_mod, MD_obs, col="red", cex=0.5, pch=19)
        abline(lm(MD_obs~MD_mod), col="red", lty=2)
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