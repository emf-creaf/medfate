evaluation_table<-function(out, measuredData, type = "SWC", cohort = NULL, 
                           temporalResolution = "day", SpParams = NULL) {
  
  # Check arguments
  temporalResolution = match.arg(temporalResolution, c("day", "week", "month", "year"))
  if("spwbInput" %in% names(out)) {
    modelInput<-out[["spwbInput"]]
    type = match.arg(type, c("SWC", "REW","E", "ETR", "SE+TR", "WP", "FMC"))
  } else {
    modelInput<- out[["growthInput"]]
    type = match.arg(type, c("SWC", "REW","E", "ETR", "SE+TR", "WP", "FMC", "BAI"))
  }
  if(type=="SWC") {
    sm = out$Soil
    d = rownames(sm)
    fc = soil_thetaFC(modelInput$soil, model = modelInput$control$soilFunctions)
    mod = sm$W.1*fc[1]
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = mod)
    
    if(!("SWC" %in% names(measuredData))) stop(paste0("Column 'SWC' not found in measured data frame."))
    seld = rownames(measuredData) %in% d
    df$Observed[d %in% rownames(measuredData)] = measuredData$SWC[seld]
    
    if("SWC_err" %in% names(measuredData))  {
      df$obs_lower[d %in% rownames(measuredData)] = df$Observed[d %in% rownames(measuredData)] - 1.96*measuredData[["SWC_err"]][seld]
      df$obs_upper[d %in% rownames(measuredData)] = df$Observed[d %in% rownames(measuredData)] + 1.96*measuredData[["SWC_err"]][seld]
    }
  } 
  else if(type=="REW") {
    sm = out$Soil
    d = rownames(sm)
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = sm$W.1)
    
    if(!("SWC" %in% names(measuredData))) stop(paste0("Column 'SWC' not found in measured data frame."))
    seld = rownames(measuredData) %in% d
    q_obs = quantile(measuredData$SWC[seld], p=c(0.90), na.rm=T) # To avoid peaks over field capacity
    df$Observed[d %in% rownames(measuredData)] = measuredData$SWC[seld]/q_obs[1]
    
  } 
  else if(type=="E") {
    pt = out$Plants$Transpiration
    d = rownames(pt)
    LAI = modelInput$above$LAI_live
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = pt[,icoh]/LAI[icoh])
    ## Fill observed values
    obscolumn = paste0("E_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }
  else if(type=="ETR") {
    ET2 = out$WaterBalance$Evapotranspiration
    d = rownames(out$WaterBalance)
    df = data.frame(Dates = as.Date(d), Observed = NA, Modelled = ET2)
    
    if(!("ETR" %in% names(measuredData))) stop(paste0("Column 'ETR' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData$ETR[rownames(measuredData) %in% d]
  }
  else if(type=="SE+TR") {
    ET1 = out$WaterBalance$SoilEvaporation+out$WaterBalance$Transpiration
    d = rownames(out$WaterBalance)
    df = data.frame(Dates = as.Date(d), Observed = NA, Modelled = ET1)
    
    if(!("ETR" %in% names(measuredData))) stop(paste0("Column 'ETR' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData$ETR[rownames(measuredData) %in% d]
  }
  else if(type=="WP"){
    wtMD = out$Plants$LeafPsiMin
    wtPD = out$Plants$LeafPsiMax
    d = rownames(wtMD)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    pdcolumn = paste0("PD_", cohort)
    mdcolumn = paste0("MD_", cohort)
    pderrcolumn = paste0("PD_", cohort, "_err")
    mderrcolumn = paste0("MD_", cohort, "_err")
    icoh = which(allcohnames==cohort)
    
    df = data.frame(Dates = as.Date(d), 
                    PD_obs = NA, MD_obs = NA, PD_obs_lower = NA, PD_obs_upper = NA, MD_obs_lower = NA, MD_obs_upper = NA,
                    PD_mod = wtPD[,icoh], MD_mod = wtMD[,icoh])
    
    seld = rownames(measuredData) %in% d
    if(pdcolumn %in% names(measuredData))  df$PD_obs[d %in% rownames(measuredData)] = measuredData[[pdcolumn]][seld]
    if(mdcolumn %in% names(measuredData))  df$MD_obs[d %in% rownames(measuredData)] = measuredData[[mdcolumn]][seld]
    if(pderrcolumn %in% names(measuredData))  {
      df$PD_obs_lower[d %in% rownames(measuredData)] = df$PD_obs[d %in% rownames(measuredData)] - 1.96*measuredData[[pderrcolumn]][seld]
      df$PD_obs_upper[d %in% rownames(measuredData)] = df$PD_obs[d %in% rownames(measuredData)] + 1.96*measuredData[[pderrcolumn]][seld]
    }
    if(mderrcolumn %in% names(measuredData))  {
      df$MD_obs_lower[d %in% rownames(measuredData)] = df$MD_obs[d %in% rownames(measuredData)] - 1.96*measuredData[[mderrcolumn]][seld]
      df$MD_obs_upper[d %in% rownames(measuredData)] = df$MD_obs[d %in% rownames(measuredData)] + 1.96*measuredData[[mderrcolumn]][seld]
    }
  }
  else if(type=="FMC") {
    fmc = moisture_cohortFMC(out, SpParams)
    d = rownames(fmc)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = fmc[,icoh])
    ## Fill observed values
    obscolumn = paste0("FMC_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }
  else if(type=="BAI") {
    SAg = out$PlantGrowth$SAgrowth
    SA = out$PlantStructure$SapwoodArea
    d = rownames(SAg)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = SAg[,icoh]*SA[,icoh])
    ## Fill observed values
    obscolumn = paste0("BAI_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }  
  if(temporalResolution != "day") {
    d.cut = cut(as.Date(d), breaks=temporalResolution)
    if(type %in% c("SWC", "REW", "FMC")) {
      df = data.frame(Dates = as.Date(levels(d.cut)),
                      Observed = tapply(df$Observed, d.cut, FUN = mean, na.rm = TRUE),
                      Modelled = tapply(df$Modelled, d.cut, FUN = mean, na.rm = TRUE))
    } else if(type == "WP") {
      df = data.frame(Dates = as.Date(levels(d.cut)),
                      PD_obs = tapply(df$PD_obs, d.cut, FUN = mean, na.rm = TRUE),
                      PD_obs_lower = tapply(df$PD_obs_lower, d.cut, FUN = mean, na.rm = TRUE),
                      PD_obs_upper = tapply(df$PD_obs_upper, d.cut, FUN = mean, na.rm = TRUE),
                      MD_obs = tapply(df$MD_obs, d.cut, FUN = mean, na.rm = TRUE),
                      MD_obs_lower = tapply(df$MD_obs_lower, d.cut, FUN = mean, na.rm = TRUE),
                      MD_obs_upper = tapply(df$MD_obs_upper, d.cut, FUN = mean, na.rm = TRUE),
                      PD_mod = tapply(df$PD_mod, d.cut, FUN = mean, na.rm = TRUE),
                      MD_mod = tapply(df$MD_mod, d.cut, FUN = mean, na.rm = TRUE))
    } else { # E, ETR, BAI
      df = data.frame(Dates = as.Date(levels(d.cut)),
                      Observed = tapply(df$Observed, d.cut, FUN = sum, na.rm = TRUE),
                      Modelled = tapply(df$Modelled, d.cut, FUN = sum, na.rm = TRUE))
    }
  }
  
  return(df)
}

evaluation_stats<-function(out, measuredData, type="SWC", cohort = NULL, 
                           temporalResolution = "day", SpParams = NULL) {
  evalstats<-function(obs, pred) {
    sel_complete = !(is.na(obs) | is.na(pred))
    obs = obs[sel_complete]
    pred = pred[sel_complete]
    E <- pred-obs
    Bias <- mean(E)
    MAE <- mean(abs(E))
    r<- cor(obs, pred)
    NSE <- 1 - (sum((obs-pred)^2)/sum((obs-mean(obs))^2))
    NSEabs <- 1 - (sum(abs(obs-pred))/sum(abs(obs-mean(obs))))
    return(c(n = sum(sel_complete), Bias= Bias, MAE = MAE, r = r, NSE = NSE, NSEabs = NSEabs))
  }
  
  # Check arguments
  if("spwbInput" %in% names(out)) {
    modelInput<-out[["spwbInput"]]
    type = match.arg(type, c("SWC", "REW","E", "ETR","SE+TR", "WP", "FMC"))
  } else {
    modelInput<- out[["growthInput"]]
    type = match.arg(type, c("SWC", "REW","E", "ETR","SE+TR", "WP", "FMC", "BAI"))
  }
  
  df = evaluation_table(out = out, measuredData = measuredData, 
                        type = type, cohort = cohort, 
                        temporalResolution = temporalResolution, SpParams = SpParams)
  if(type=="SWC") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="REW") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="E") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="FMC") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="ETR") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="SE+TR") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="WP") {
    eval_res = as.data.frame(rbind(evalstats(df$PD_obs, df$PD_mod),
                                   evalstats(df$MD_obs, df$MD_mod)))
    row.names(eval_res)<-c("Predawn potentials", "Midday potentials")
  }
  else if(type=="BAI") eval_res = evalstats(df$Observed, df$Modelled)
  return(eval_res)
}

evaluation_plot<-function(out, measuredData, type="SWC", cohort = NULL, 
                          temporalResolution = "day", SpParams = NULL, 
                          plotType = "dynamics") {
  scatterplot<-function(df, xlab="", ylab="", title=NULL, err = FALSE) {
    g<-ggplot(df, aes_string(x="Modelled"))
    if(err) {
      g<-g+
        geom_pointrange(aes_string(y = "Observed", ymin = "obs_lower", ymax = "obs_upper"),cex=0.5)
    }
    g<-g + 
      geom_point(aes_string(y = "Observed"), cex=0.5)+
      geom_abline(intercept=0, slope=1, col="black")+
      geom_smooth(aes_string(y = "Observed"), method="lm", se = FALSE, col="gray", linetype="dashed")+
      xlab(xlab)+
      ylab(ylab)+
      theme_bw()
    if(!is.null(title)) g<-g+labs(title=title)
    return(g)
  }
  dynamicsplot<-function(df, xlab="", ylab="", title=NULL, err = FALSE,
                         str_obs = "Observed", str_mod = "Modelled") {
    g<-ggplot(df, aes_string(x="Dates"))
    if(err) {
      g <- g +          
        geom_ribbon(aes_(ymin=~obs_lower, ymax=~obs_upper), 
                    col="gray", alpha= 0.5)
    }
    g<-g+       
      geom_path(aes_(y=~Observed, col="Observed"))+
      geom_path(aes_(y=~Modelled, col="Modelled"))+
      xlab(xlab)+
      ylab(ylab)+
      scale_color_manual(name="", 
                         values=c("Observed"="black", "Modelled"= "red"),
                         labels =c(str_obs, str_mod))+
      theme_bw()
    if(!is.null(title)) g<-g+labs(title=title)
    return(g)
  }

  # Check arguments
  plotType = match.arg(plotType, c("dynamics", "scatter"))
  if("spwbInput" %in% names(out)) {
    modelInput<-out[["spwbInput"]]
    type = match.arg(type, c("SWC", "REW","E", "ETR","SE+TR", "WP", "FMC"))
  } else {
    modelInput<- out[["growthInput"]]
    type = match.arg(type, c("SWC", "REW","E", "ETR","SE+TR", "WP", "FMC", "BAI"))
  }
  
  df = evaluation_table(out = out, measuredData = measuredData, 
                        type = type, cohort = cohort, 
                        temporalResolution = temporalResolution, SpParams = SpParams)
  
  if(type=="SWC") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = expression(paste("Soil moisture ",(m^{3}%.%m^{-3}))),
                      err = ("SWC_err" %in% names(measuredData)))
    } else {
      g<-scatterplot(df, xlab  = expression(paste("Measured soil moisture ",(m^{3}%.%m^{-3}))),
                     ylab = expression(paste("Measured soil moisture ",(m^{3}%.%m^{-3}))), 
                     err = ("SWC_err" %in% names(measuredData)))
    }
  } 
  else if(type=="REW") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Relative extractable soil water (REW)")
    } else {
      g<-scatterplot(df, xlab  = "Modelled relative extractable soil water (REW)",
                     ylab = "Measured relative extractable soil water (REW)")
      
    }
  }
  else if(type=="E") {
    allcohnames = row.names(modelInput$cohorts)
    spnames = modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Transpiration per leaf area (l/m2/day)", 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = "Modelled transpiration per leaf area (l/m2/day)",
                     ylab = "Measured transpiration per leaf area (l/m2/day)", 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="FMC") {
    allcohnames = row.names(modelInput$cohorts)
    spnames = modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Fuel moisture content (% of dry weight)", 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = "Modelled fuel moisture content (% of dry weight)",
                     ylab = "Measured fuel moisture content (% of dry weight)", 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="BAI") {
    allcohnames = row.names(modelInput$cohorts)
    spnames = modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = paste0("Basal area increment (cm2/", temporalResolution,")"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = paste0("Modelled basal area increment (cm2/", temporalResolution,")"),
                     ylab = paste0("Measured basal area increment (cm2/", temporalResolution,")"), 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="ETR") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "ETR (mm)")
    } else {
      g<-scatterplot(df, xlab  = "Modelled ETR (mm)",
                         ylab ="Measured ETR (mm)")
    }
  }
  else if(type=="SE+TR") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "ETR or SE+TR (mm)", str_obs = "Observed (ETR)", str_mod = "Modelled (SE+TR)")
    } else {
      g<-scatterplot(df, xlab  = "Modelled SE+TR (mm)",
                     ylab ="Measured ETR (mm)")
    }
  }
  else if(type=="WP"){
    wtMD = out$Plants$LeafPsiMin
    wtPD = out$Plants$LeafPsiMax
    d = rownames(wtMD)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics"){
      g<-ggplot(df)+
        geom_path(aes_(x=~Dates, y=~PD_mod, col="Predawn", linetype="Predawn"))+
        geom_path(aes_(x=~Dates, y=~MD_mod, col="Midday", linetype="Midday"))+
        geom_pointrange(aes_(x = ~Dates, y = ~PD_obs, ymin = ~PD_obs_lower, ymax = ~PD_obs_upper, col="Predawn", linetype="Predawn"))+
        geom_pointrange(aes_(x = ~Dates, y = ~MD_obs, ymin = ~MD_obs_lower, ymax = ~MD_obs_upper, col="Midday",linetype="Midday"))+
        scale_color_manual(name="", values=c("Predawn"="blue", "Midday"= "red"))+
        scale_linetype_manual(name="", values=c("Predawn"="dashed", "Midday"= "solid"))+
        labs(title=paste0(cohort , " (",spnames[icoh],")"))+
        xlab("")+
        ylab("Leaf water potential (MPa)")+
        theme_bw()
    } else {
      g<-ggplot(df)+
        geom_abline(intercept=0, slope=1, col="black")+
        geom_pointrange(aes_(x = ~PD_mod, y = ~PD_obs, ymin = ~PD_obs_lower, ymax = ~PD_obs_upper, col="Predawn"))+
        geom_pointrange(aes_(x = ~MD_mod, y = ~MD_obs, ymin = ~MD_obs_lower, ymax = ~MD_obs_upper,col="Midday"))+
        geom_smooth(aes_(x = ~PD_mod, y = ~PD_obs, col="Predawn"), method="lm", se = FALSE, linetype="dashed")+
        geom_smooth(aes_(x = ~MD_mod, y = ~MD_obs, col="Midday"), method="lm", se = FALSE, linetype="dashed")+
        scale_color_manual(name="", values=c("Predawn"="blue", "Midday"= "red"))+
        labs(title=paste0(cohort , " (",spnames[icoh],")"))+
        xlab("Modelled leaf water potential (MPa)")+
        ylab("Measured leaf water potential (MPa)")+
        theme_bw()
    }
  }
  return(g)
}

evaluation_metric<-function(out, measuredData, type="SWC", cohort=NULL, 
                            temporalResolution = "day", SpParams = NULL,
                            metric = "loglikelihood") {
  df = evaluation_table(out = out, measuredData = measuredData, 
                        type = type, cohort = cohort, 
                        temporalResolution = temporalResolution, SpParams = SpParams)
  obs = df$Observed
  pred = df$Modelled
  sd <- sd(obs, na.rm=TRUE)
  metric<-match.arg(metric, c("loglikelihood", "NSE", "NSEabs", "MAE", "r"))
  m <- switch(metric,
         "loglikelihood" = sum(dnorm(obs, pred, sd, log=TRUE), na.rm=TRUE),
         "NSE" = 1 - (sum((obs-pred)^2, na.rm=TRUE)/sum((obs-mean(obs, na.rm=TRUE))^2, na.rm=TRUE)),
         "MAE" = mean(abs(pred-obs), na.rm=TRUE),
         "r" = cor(obs, pred),
         "NSEabs" = 1 - (sum(abs(obs-pred))/sum(abs(obs-mean(obs))))
  )
  return(m)
}