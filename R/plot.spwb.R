plot.spwb<-function(x, type="PET_Precipitation", bySpecies = FALSE,
                   yearAxis=FALSE, xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                   add=FALSE, ...) {
  dates = as.Date(rownames(x$WaterBalance))
  transpMode = x$Input$control$transpirationMode
  WaterBalance = x$WaterBalance
  Soil = x$Soil
  nlayers = x$NumSoilLayers
  TYPES = c("PET_Precipitation","PET_NetRain","Snow","Evapotranspiration","SoilPsi","SoilTheta","SoilVol", "Export", "LAI", "WTD",
            "PlantLAI",
            "PlantStress", "PlantPsi","PlantPhotosynthesis", "PlantTranspiration",
            "PlantPhotosynthesisLeaf","PlantTranspirationLeaf")
  if(transpMode=="Complex") {
    TYPES = c("PET_Precipitation","PET_NetRain","Snow","Evapotranspiration","SoilPsi","SoilTheta","SoilVol", "Export", "LAI", "WTD",
              "PlantLAI",
              "SoilPlantConductance",
              "PlantStress", "PlantPhotosynthesis", "PlantTranspiration",
              "PlantPhotosynthesisLeaf","PlantTranspirationLeaf", 
              "LeafPsi","StemPsi","RootPsi","StemRWC", "LeafRWC", "PlantWaterBalance",
              "PlantAbsorbedSWR", "PlantAbsorbedSWRLeaf",
              "PlantAbsorbedLWR", "PlantAbsorbedLWRLeaf",
              "AirTemperature","SoilTemperature", "CanopyTemperature",
              "CanopyEnergyBalance", "SoilEnergyBalance")
  } 
  type = match.arg(type,TYPES)  
  numDays = length(dates)
  numYears = round(numDays/365)
  firstYear=as.numeric(format(dates[1],"%Y"))
  cohortnames = colnames(x$PlantTranspiration)
  plotAxes<-function(){
    if(!yearAxis) axis.Date(1, dates)
    else {
      axis(1, at = (0:numYears)*365, labels=FALSE)
      axis(1, at = -182+365*(1:numYears), tick = FALSE, line=FALSE, labels=firstYear:(firstYear+numYears-1))
    }
    axis(2)    
  }
  mnp = max(WaterBalance$Precipitation)
  if(is.null(xlab)) xlab = ifelse(yearAxis,"Year", "Date")  
  if(type=="PET_Precipitation") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})
    if(!is.null(xlim)) span = xlim[1]:xlim[2]
    else span = 1:numDays
    if(is.null(ylim)) ylim = c(0,mnp)
    barplot(WaterBalance$Precipitation[span], ylim=ylim, col="black",space=0, ylab=ylab, 
            xlab=xlab, axes=FALSE)
    barplot(WaterBalance$Snow[span], col="red", border = "red", add=TRUE)
    plotAxes()
    lines(1:length(span), WaterBalance$PET[span], col="gray")    
    legend("topleft", bty="n", col=c("black","red", "gray"),lty=1, lwd=2,
           legend=c("Precipitation","Snow", "PET"))
    
  } 
  else if(type=="PET_NetRain") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})    
    if(!is.null(xlim)) span = xlim[1]:xlim[2]
    else span = 1:numDays
    if(is.null(ylim)) ylim = c(0,mnp)
    barplot(WaterBalance$NetRain[span], ylim=ylim, col="black",space=0, ylab=ylab, 
            xlab=xlab, axes=FALSE)
    plotAxes()
    lines(1:length(span), WaterBalance$PET[span], col="gray")    
    legend("topleft", bty="n", col=c("black","gray"),lty=c(1,1), lwd=2,
           legend=c("NetRain","PET"))        
  } 
  else if(type=="Snow") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})    
    if(!is.null(xlim)) span = xlim[1]:xlim[2]
    else span = 1:numDays
    mnp = max(c(WaterBalance$Snow[span], Soil$SWE[span]))
    if(is.null(ylim)) ylim = c(0,mnp)
    barplot(WaterBalance$Snow[span], ylim=ylim, col="black",space=0, ylab=ylab, 
            xlab=xlab, axes=FALSE)
    plotAxes()
    lines(1:length(span), Soil$SWE[span], col="red", lwd=1.5)    
    legend("topleft", bty="n", col=c("black","red"),lty=c(1,1), lwd=2,
           legend=c("Snow","Snowpack (SWE)"))        
  } 
  else if(type=="Evapotranspiration") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})
    if(is.null(ylim)) ylim = c(0,max(WaterBalance$Evapotranspiration))
    plot(dates, WaterBalance$Evapotranspiration, ylim=ylim, type="l", ylab=ylab, 
         xlab=xlab, xlim=xlim,frame=FALSE, col="black", axes=FALSE, lwd=2)
    plotAxes()
    lines(dates, WaterBalance$Transpiration, col="gray", lty=2, lwd=1.5)
    lines(dates, WaterBalance$SoilEvaporation, col="black", lty=3, lwd=1.5)
    legend("topleft", bty="n", col=c("black","gray","black"),lty=c(1,2,3), lwd=c(2,1.5,1.5),
           legend=c("Total evapotranspiration","Transpiration","Bare soil evaporation"))
  } 
  else if(type=="PlantWaterBalance") {
    pwb = WaterBalance$PlantExtraction - WaterBalance$Transpiration
    if(is.null(ylab)) ylab = expression(paste("Extraction - transpiration (",L%.%m^{-2},")"))
    if(is.null(ylim)) ylim = c(min(pwb,na.rm=T),max(pwb,na.rm=T))
    plot(dates, pwb, ylim=ylim, type="l", ylab=ylab, 
         xlab=xlab, xlim=xlim,frame=FALSE, col="black", axes=FALSE, lwd=2)
    plotAxes()
  } 
  else if(type=="LAI") {
    if(is.null(ylab)) ylab = expression(paste("Leaf Area Index   ",(m^{2}%.%m^{-2})))
    if(is.null(ylim)) ylim = c(0,max(WaterBalance$LAIcell))
    plot(dates, WaterBalance$LAIcell, ylim=ylim, type="l", ylab=ylab, 
         xlab=xlab, xlim=xlim,frame=FALSE, col="black", axes=FALSE, lwd=1)
    plotAxes()
    lines(dates, WaterBalance$LAIcelldead, lty=2)
  } 
  else if(type=="WTD") {
    if(is.null(ylab)) ylab = expression(paste("Water table depth  (mm)"))
    if(is.null(ylim)) ylim = c(max(Soil$WTD),0)
    plot(dates, Soil$WTD, ylim=ylim, type="l", ylab=ylab, 
         xlab=xlab, xlim=xlim,frame=FALSE, col="black", axes=FALSE, lwd=1)
    plotAxes()
    lines(dates, Soil$WTD, lty=2)
  } 
  else if(type=="SoilPsi") {
    PsiM = Soil[,paste("psi",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Soil water potential (MPa)"    
    if(is.null(ylim)) ylim =c(min(PsiM),0)
    matplot(dates, PsiM, lwd=1.5,
         ylim=ylim, type="l", ylab=ylab, xlab=xlab, xlim=xlim,
         frame=FALSE, lty=c(1,2,3,4,5), col="black",axes=FALSE)
    plotAxes()
    legend("bottomleft", bty="n", lty=c(1,2,3,4,5), col="black", lwd=1.5,
           legend=paste("Layer", 1:nlayers))
    
  } 
  else if(type=="SoilTheta") {
    WM = Soil[,paste("W",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Soil moisture (% field capacity)"
    if(is.null(ylim)) ylim = c(0,100*max(WM,na.rm = T))
    matplot(dates, WM*100, lwd=1.5,
            ylim=ylim, type="l", ylab=ylab, xlab=xlab, xlim=xlim,
            frame=FALSE, lty=c(1,2,3,4,5), col="black",axes=FALSE)
    plotAxes()
    legend("bottomleft", bty="n", col="black",lty=1:nlayers, lwd=1.5,
           legend=paste("Layer", 1:nlayers))
  } 
  else if(type=="SoilVol") {
    MLM= Soil[,paste("ML",1:nlayers,sep=".")]
    MLTot = Soil$MLTot
    if(is.null(ylim)) ylim =c(0,max(MLTot)*1.3)
    if(is.null(ylab)) ylab = "Soil water content (mm)"
    plot(dates, MLTot, ylim=ylim, lwd=2, type="l",xlim=xlim,
         ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    matlines(dates, MLM, lty = 1:nlayers, col="black", lwd=1.5)
    legend("topleft", bty="n", col="black",lty=c(1,1:nlayers), lwd=c(2,rep(1.5,5)),
           legend=c("Total",paste("Layer", 1:nlayers)))
    
  } 
  else if(type=="Export") {
    if(is.null(ylab)) ylab =  expression(L%.%m^{-2})    
    mnp = max(WaterBalance$DeepDrainage+WaterBalance$Runoff)    
    if(is.null(ylim)) ylim = c(0,mnp)
    plot(dates, WaterBalance$DeepDrainage+WaterBalance$Runoff, ylim=ylim, col="black", type="l", 
         ylab=ylab, xlab=xlab, xlim=xlim,
         frame=FALSE, axes=FALSE)
    lines(dates, WaterBalance$DeepDrainage, col="blue")
    lines(dates, WaterBalance$Runoff, col="red")
    plotAxes()
    legend("topright", bty="n", col=c("black","blue","red"),lty=c(1,1,1), lwd=c(1.5,1,1),
           legend=c("DD+R","Deep drainage (DD)","Runoff (R)"))        
  } 
  else if(type=="PlantLAI") {
    OM = x$PlantLAI
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = expression(paste("Leaf Area Index   ",(m^{2}%.%m^{-2})))
    if(is.null(ylim)) ylim = c(0,max(OM, na.rm=TRUE))
    matplot(dates, OM, ylim = ylim, lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="SoilPlantConductance") {
    OM = x$dEdP
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, x$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, x$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = expression(paste("Average soil-plant conductance ",(mmol%.%m^{-2}%.%s^{-1})))
    if(is.null(ylim)) ylim = c(min(OM, na.rm=T),max(OM, na.rm=T))
    matplot(dates, OM, lty=1:length(cohortnames), col = 1:length(cohortnames),
            ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantStress") {
    OM = x$PlantStress
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, x$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, x$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = "Drought stress [0-1]"
    if(is.null(ylim)) ylim = c(0,1)
    matplot(dates, OM, lty=1:length(cohortnames), col = 1:length(cohortnames),
            ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="StemRWC") {
    OM = x$PlantRWCstem*100
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, x$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, x$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = "Relative water content in stem [%]"
    if(is.null(ylim)) ylim = c(min(OM, na.rm=T),max(OM, na.rm=T))
    matplot(dates, OM, lty=1:length(cohortnames), col = 1:length(cohortnames),
            ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="LeafRWC") {
    OM = x$PlantRWCleaf*100
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, x$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, x$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = "Relative water content in leaf [%]"
    if(is.null(ylim)) ylim = c(min(OM, na.rm=T),max(OM, na.rm=T))
    matplot(dates, OM, lty=1:length(cohortnames), col = 1:length(cohortnames),
            ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantPsi") {
    OM = x$PlantPsi
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, x$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = "Plant water potential (MPa)"
    if(is.null(ylim)) ylim = c(min(OM, na.rm = TRUE),0)
    matplot(dates, OM, ylim = ylim, lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="LeafPsi") {
    OM = x$LeafPsi
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, x$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = "Midday leaf water potential (MPa)"
    if(is.null(ylim)) ylim = c(min(OM, na.rm = TRUE),0)
    matplot(dates, OM, ylim = ylim, lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="StemPsi") {
    OM = x$StemPsi
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, x$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = "Midday (upper) stem water potential (MPa)"
    if(is.null(ylim)) ylim = c(min(OM, na.rm = TRUE),0)
    matplot(dates, OM, ylim = ylim, lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="RootPsi") {
    OM = x$RootPsi
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, x$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = "Midday root crown water potential (MPa)"
    if(is.null(ylim)) ylim = c(min(OM, na.rm = TRUE),0)
    matplot(dates, OM, ylim = ylim, lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantTranspiration") {
    OM = x$PlantTranspiration
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      cohortnames = colnames(OM)
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant transpiration   ",(L%.%m^{-2})))
    if(is.null(ylim)) ylim = c(0,max(OM, na.rm=T))
    matplot(dates, OM, ylim = ylim,
            lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()      
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantTranspirationLeaf") {
    df = x$PlantTranspiration
    if(bySpecies) {
      m1 = apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$PlantLAI,1,tapply, x$cohorts$Name, sum, na.rm=T)
      df = t(m1/lai1)
      cohortnames = colnames(df)
    } else {
      df = df/x$PlantLAI
      df[x$PlantLAI==0] = NA
    }
    if(is.null(ylab)) ylab = expression(paste("Plant transpiration per leaf area  ",(L%.%m^{-2})))
    if(is.null(ylim)) ylim = c(0,max(df, na.rm=T))
    matplot(dates, df, ylim = ylim, 
            lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()      
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantPhotosynthesis") {
    df = x$PlantPhotosynthesis
    if(bySpecies) {
      df = t(apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T))
      cohortnames = colnames(df)
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant photosynthesis   ",(g*C%.%m^{-2})))
    if(is.null(ylim)) ylim = c(min(df, na.rm=T),max(df, na.rm=T))
    matplot(dates, df, ylim = ylim, 
            lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()     
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantPhotosynthesisLeaf") {
    df = x$PlantPhotosynthesis
    if(bySpecies) {
      m1 = apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$PlantLAI,1,tapply, x$cohorts$Name, sum, na.rm=T)
      df = t(m1/lai1)
      cohortnames = colnames(df)
    } else {
      df = df/x$PlantLAI
      df[x$PlantLAI==0] = NA
    }
    if(is.null(ylab)) ylab = expression(paste("Plant photosynthesis per leaf area   ",(g*C%.%m^{-2})))
    if(is.null(ylim)) ylim = c(min(df, na.rm=T),max(df, na.rm=T))
    matplot(dates, df, ylim = ylim, 
            lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()     
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantAbsorbedSWR") {
    df = x$PlantAbsorbedSWR
    if(bySpecies) {
      df = t(apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T))
      cohortnames = colnames(df)
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant absorbed SWR  ",(MJ%.%m^{-2})))
    if(is.null(ylim)) ylim = c(0,max(df))
    matplot(dates, df, 
            lty=1:length(cohortnames), col = 1:length(cohortnames), 
            ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()     
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantAbsorbedSWRLeaf") {
    df = x$PlantAbsorbedSWR
    if(bySpecies) {
      m1 = apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$PlantLAI,1,tapply, x$cohorts$Name, sum, na.rm=T)
      df = t(m1/lai1)
      cohortnames = colnames(df)
    } else {
      df = df/x$PlantLAI
      df[x$PlantLAI==0] = NA
    }
    if(is.null(ylab)) ylab = expression(paste("Plant absorbed SWR per leaf area  ",(MJ%.%m^{-2})))
    if(is.null(ylim)) ylim = c(min(df, na.rm=T),max(df, na.rm=T))
    matplot(dates, df, ylim = ylim, lwd=1, type="l", 
            lty=1:length(cohortnames), col = 1:length(cohortnames),xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()     
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantAbsorbedLWR") {
    df = x$PlantAbsorbedLWR
    if(bySpecies) {
      df = t(apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T))
      cohortnames = colnames(df)
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant absorbed LWR  ",(MJ%.%m^{-2})))
    if(is.null(ylim)) ylim = c(0,max(df))
    matplot(dates, df, ylim = ylim, 
            lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()     
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="PlantAbsorbedLWRLeaf") {
    df = x$PlantAbsorbedLWR
    if(bySpecies) {
      m1 = apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$PlantLAI,1,tapply, x$cohorts$Name, sum, na.rm=T)
      df = t(m1/lai1)
      cohortnames = colnames(df)
    } else {
      df = df/x$PlantLAI
      df[x$PlantLAI==0] = NA
    }
    if(is.null(ylab)) ylab = expression(paste("Plant absorbed LWR per leaf area  ",(MJ%.%m^{-2})))
    if(is.null(ylim)) ylim = c(min(df, na.rm=T),max(df, na.rm=T))
    matplot(dates, df, ylim = ylim, 
            lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()     
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  } 
  else if(type=="AirTemperature") {
    if(is.null(ylab)) ylab = "Above-canopy temperature (Celsius)"
    if(is.null(ylim)) ylim = c(min(x$Temperature$Tatm_min),max(x$Temperature$Tatm_max))
    if(!add) {
      plot(dates, x$Temperature$Tatm_mean, ylim = ylim, type="l", xlim=xlim,
         ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE, ...)
      plotAxes()   
    } else {
      lines(dates, x$Temperature$Tatm_min, col="black", ...)
    }
    lines(dates, x$Temperature$Tatm_min, col="blue", ...)
    lines(dates, x$Temperature$Tatm_max, col="red", ...)
  } 
  else if(type=="CanopyTemperature") {
    if(is.null(ylab)) ylab = "Canopy temperature (Celsius)"
    if(is.null(ylim)) ylim = c(min(x$Temperature$Tcan_min),max(x$Temperature$Tcan_max))
    if(!add) {
      plot(dates, x$Temperature$Tcan_mean, ylim = ylim, type="l", xlim=xlim,
         ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE, ...)
      plotAxes()   
    } else {
      lines(dates, x$Temperature$Tcan_mean, col="black", ...)
    }
    lines(dates, x$Temperature$Tcan_min, col="blue", ...)
    lines(dates, x$Temperature$Tcan_max, col="red", ...)
  } 
  else if(type=="SoilTemperature") {
    if(is.null(ylab)) ylab = "Soil temperature (Celsius)"
    if(is.null(ylim)) ylim = c(min(x$Temperature$Tsoil_min),max(x$Temperature$Tsoil_max))
    if(!add) {
      plot(dates, x$Temperature$Tsoil_mean, ylim = ylim, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE, ...)
      plotAxes()   
    } else {
      lines(dates, x$Temperature$Tsoil_mean, col="black", ...)
    }
    lines(dates, x$Temperature$Tsoil_min, col="blue", ...)
    lines(dates, x$Temperature$Tsoil_max, col="red", ...)
  } 
  else if(type=="CanopyEnergyBalance") {
    if(is.null(ylab)) ylab = expression(MJ%.%m^{-2})    
    mxin = max(c(x$EnergyBalance$SWRcanin,x$EnergyBalance$LWRcanin,x$EnergyBalance$LWRsoilcan,-x$EnergyBalance$LEcan,-x$EnergyBalance$Hcan))    
    mxout = max(c(x$EnergyBalance$LWRcanout,x$EnergyBalance$LEcan,x$EnergyBalance$Hcan))    
    if(is.null(ylim)) ylim = c(-mxout,mxin)
    plot(dates, x$EnergyBalance$Ebalcan, ylim=ylim, type="n", 
         ylab=ylab, xlab=xlab, xlim=xlim,
         frame=FALSE, axes=FALSE,...)
    plotAxes()
    abline(h=0, col="black", lwd=1.5)
    lines(dates, x$EnergyBalance$Ebalcan, col="black",...)
    lines(dates, x$EnergyBalance$SWRcanin, col="red",...)
    lines(dates, x$EnergyBalance$LWRcanin, col="brown",...)
    lines(dates, -x$EnergyBalance$LWRcanout, col="blue",...)
    lines(dates, x$EnergyBalance$LWRsoilcan, col="orange",...)
    lines(dates, -x$EnergyBalance$LEcan, col="green",...)
    lines(dates, -x$EnergyBalance$Hcan, col="gray",...)
    lines(dates, -x$EnergyBalance$Hcansoil, col="dark gray",...)
    legend("topright", bty="n", col=c("red","brown","orange", "blue","green", "gray", "dark gray", "black"), lty=1,
           legend=c("SWR abs. from atm.","LWR abs. from atm.","LWR abs. from soil","LWR emmited", "Latent heat (L)",
                    "Convection can./atm.","Convection soil/can.", "Balance"),...)        
  } 
  else if(type=="SoilEnergyBalance") {
    if(is.null(ylab)) ylab = expression(MJ%.%m^{-2})    
    mxin = max(c(x$EnergyBalance$SWRsoilin, x$EnergyBalance$LWRsoilin,x$EnergyBalance$LWRcanout,x$EnergyBalance$Hcansoil))    
    mxout = max(c(x$EnergyBalance$LWRsoilout,-x$EnergyBalance$Hcansoil))    
    if(is.null(ylim)) ylim = c(-mxout,mxin)
    plot(dates, x$EnergyBalance$Ebalsoil, ylim=ylim, type="n", 
         ylab=ylab, xlab=xlab, xlim=xlim,
         frame=FALSE, axes=FALSE,...)
    plotAxes()
    abline(h=0, col="black", lwd=1.5)
    lines(dates, x$EnergyBalance$Ebalsoil, col="black",...)
    lines(dates, x$EnergyBalance$SWRsoilin, col="red",...)
    lines(dates, x$EnergyBalance$LWRsoilin, col="brown",...)
    lines(dates, x$EnergyBalance$LWRcanout, col="orange",...)
    lines(dates, -x$EnergyBalance$LWRsoilout, col="blue",...)
    lines(dates, x$EnergyBalance$Hcansoil, col="gray",...)
    legend("topright", bty="n", col=c("red","brown","orange", "blue", "gray", "black"), lty=1,
           legend=c("SWR abs. from atm.","LWR abs. from atm.", "LWR abs. from canopy","LWR emmited", "Convection soil/can.", "Balance"),...)        
  }
}
