evaluation_table<-function(out, measuredData, type = "SWC", cohort = NULL) {
  
  # Check arguments
  type = match.arg(type, c("SWC", "REW","E", "ETR", "WP"))
  if(type=="SWC") {
    sm = out$Soil
    d = rownames(sm)
    fc = soil_thetaFC(out$soilInput, model = out$spwbInput$control$soilFunctions)
    df <- data.frame(Observed = NA, Modelled = sm$W.1*fc[1], Dates = as.Date(d))
    
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
    # fc = soil_thetaFC(out$soilInput, model = out$spwbInput$control$soilFunctions)
    # q_mod = quantile(sm$W.1, p=c(0.05,0.95), na.rm=T)
    # df <- data.frame(Observed = NA, Modelled = (sm$W.1-q_mod[1])/(q_mod[2]-q_mod[1]), Dates = as.Date(d))
    df <- data.frame(Observed = NA, Modelled = sm$W.1, Dates = as.Date(d))
    
    if(!("SWC" %in% names(measuredData))) stop(paste0("Column 'SWC' not found in measured data frame."))
    seld = rownames(measuredData) %in% d
    # q_obs = quantile(measuredData$SWC[seld], p=c(0.05,0.95), na.rm=T)
    q_obs = quantile(measuredData$SWC[seld], p=c(0.9), na.rm=T) # To avoid peaks over field capacity
    # df$Observed[d %in% rownames(measuredData)] =(measuredData$SWC[seld]-q_obs[1])/(q_obs[2]-q_obs[1])
    df$Observed[d %in% rownames(measuredData)] = measuredData$SWC[seld]/q_obs[1]
    
  } 
  else if(type=="E") {
    pt = out$Plants$Transpiration
    d = rownames(pt)
    LAI = out$spwbInput$above$LAI_live
    spnames = out$spwbInput$cohorts$Name
    allcohnames = row.names(out$spwbInput$cohorts)
    
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
  }
  else if(type=="ETR") {
    ET1 = out$WaterBalance$SoilEvaporation+out$WaterBalance$Transpiration
    ET2 = out$WaterBalance$Evapotranspiration
    d = rownames(out$WaterBalance)
    df = data.frame(ETobs = NA, ET1 = ET1, ET2 = ET2, Dates = as.Date(d))
    
    if(!("ETR" %in% names(measuredData))) stop(paste0("Column 'ETR' not found in measured data frame."))
    df$ETobs[d %in% rownames(measuredData)] = measuredData$ETR[rownames(measuredData) %in% d]
  }
  else if(type=="WP"){
    wtMD = out$Plants$LeafPsiMin
    wtPD = out$Plants$LeafPsiMax
    d = rownames(wtMD)
    spnames = out$spwbInput$cohorts$Name
    allcohnames = row.names(out$spwbInput$cohorts)
    
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
    
    df = data.frame(PD_obs = NA, MD_obs = NA, PD_obs_lower = NA, PD_obs_upper = NA, MD_obs_lower = NA, MD_obs_upper = NA,
                    PD_mod = wtPD[,icoh], MD_mod = wtMD[,icoh], Dates = as.Date(d))
    
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
  return(df)
}

evaluation_stats<-function(out, measuredData, type="SWC", cohort = NULL) {
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
  type = match.arg(type, c("SWC", "REW","E", "ETR", "WP"))
  
  df = evaluation_table(out, measuredData, type, cohort)
  if(type=="SWC") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="REW") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="E") eval_res = evalstats(df$Observed, df$Modelled)
  else if(type=="ETR") {
    eval_res = as.data.frame(rbind(evalstats(df$ETobs, df$ET1),
                                   evalstats(df$ETobs, df$ET2)))
    row.names(eval_res)<-c("Es+Tr", "Es+Tr+In")
  }
  else if(type=="WP") {
    eval_res = as.data.frame(rbind(evalstats(df$PD_obs, df$PD_mod),
                                   evalstats(df$MD_obs, df$MD_mod)))
    row.names(eval_res)<-c("Predawn potentials", "Midday potentials")
  }
  
  return(eval_res)
}

evaluation_plot<-function(out, measuredData, type="SWC", cohort = NULL, 
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
  dynamicsplot<-function(df, xlab="", ylab="", title=NULL, err = FALSE) {
    g<-ggplot(df, aes_string(x="Dates"))
    if(err) {
      g <- g +          
        geom_ribbon(aes_(ymin=~obs_lower, ymax=~obs_upper), col="gray", alpha= 0.5)
    }
    g<-g+       
      geom_path(aes_(y=~Observed, col="Observed"))+
      geom_path(aes_(y=~Modelled, col="Modelled"))+
      xlab(xlab)+
      ylab(ylab)+
      scale_color_manual(name="", values=c("Observed"="black", "Modelled"= "red"))+
      theme_bw()
    if(!is.null(title)) g<-g+labs(title=title)
    return(g)
  }

  # Check arguments
  type = match.arg(type, c("SWC", "REW","E", "ETR", "WP"))
  plotType = match.arg(plotType, c("dynamics", "scatter"))
  
  
  df = evaluation_table(out, measuredData, type, cohort)
  
  if(type=="SWC") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = expression(paste("Soil moisture ",(m^{3}%.%m^{-3}))), err = ("SWC_err" %in% names(measuredData)))
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
    allcohnames = row.names(out$spwbInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Transpiration per leaf area (mm/day)", 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = "Modelled transpiration per leaf area (mm/day)",
                     ylab = "Measured transpiration per leaf area (mm/day)", 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="ETR") {
    ETobs = df$ETobs
    if(plotType=="dynamics") {
      g<-ggplot(df, aes_string(x="Dates"))+
        geom_path(aes(y=ETobs, col="Measured ETR"))+
        geom_path(aes(y=ET1, col="Modelled Es+Tr"))+
        geom_path(aes(y=ET2, col="Modelled Es+Tr+In"))+
        xlab("")+
        ylab("ETR (mm)")+
        scale_color_manual(name="", values=c("Measured ETR"="red", 
                                             "Modelled Es+Tr" = "black",
                                             "Modelled Es+Tr+In"= "gray"))+
        theme_bw()
    } else {
      g<-ggplot(df, aes_string(y = "ETobs"))+
        geom_abline(intercept=0, slope=1, col="black")+
        geom_point(aes_string(x = "ET1"), col="black", cex=0.5)+
        geom_point(aes_string(x = "ET2"), col="gray", cex=0.5)+
        geom_smooth(aes_string(x = "ET1"), method="lm", se = FALSE, col="black", linetype="dashed")+
        geom_smooth(aes_string(x = "ET2"), method="lm", se = FALSE, col="gray", linetype="dashed")+
        xlab("Modelled ETR (mm)")+
        ylab("Measured ETR (mm)")+
        theme_bw()
    }
  }
  else if(type=="WP"){
    wtMD = out$Plants$LeafPsiMin
    wtPD = out$Plants$LeafPsiMax
    d = rownames(wtMD)
    spnames = out$spwbInput$cohorts$Name
    allcohnames = row.names(out$spwbInput$cohorts)
    
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

evaluation_loglikelihood<-function(out, measuredData, type="SWC", cohort=NULL) {
  df = evaluation_table(out, measuredData, type=type, cohort = cohort)
  sd = sd(df$Modelled, na.rm=T)
  ll = sum(dnorm(df$Observed, df$Modelled, sd, log=T), na.rm=T)
  return(ll)
}

evaluation_loglikelihood_function<-function(x, soil,
                                            cohNames, parNames, 
                                            measuredData, type = "SWC", cohort = NULL, 
                                            ...) {
  l = list(...)
  
  if(inherits(x, "spwbInput")) model = "spwb"
  else model = "growth"
  
  custom = data.frame(Cohort = unique(cohNames), as.data.frame(matrix(NA, nrow = length(unique(cohNames)), ncol = length(parNames))))
  names(custom) = c("Cohort", parNames)
  x$control$verbose = FALSE
  
  llf<-function(v) {
    resetInputs(x, soil)
    for(i in 1:length(parNames)) custom[custom$Cohort==cohNames[i], parNames[i]] = v[i]
    x = modifyCohortParams(x, custom, soil)
    if(model=="spwb") {
      S = spwb(x, soil, 
               meteo = l[["meteo"]], 
               latitude = l[["latitude"]], elevation = l[["elevation"]],
               slope  = l[["slope"]],aspect = l[["aspect"]])
    } 
    else {
      S = growth(x, soil, 
                 meteo = l[["meteo"]], 
                 latitude = l[["latitude"]], elevation = l[["elevation"]],
                 slope  = l[["slope"]], aspect = l[["aspect"]])
    }
    ll = evaluation_loglikelihood(S, measuredData, type, cohort)
    return(ll)
  }
  return(llf)
}