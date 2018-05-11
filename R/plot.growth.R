plot.growth<-function(x, type="PET_Precipitation", bySpecies = FALSE, 
                      yearAxis=FALSE, xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,...) {
  TYPES_SWB = c("PET_Precipitation","PET_NetPrec","ET","Psi","Theta","Vol", "LAI", "PlantLAI",
                "PlantStress", "PlantPsi","PlantPhotosynthesis","PlantTranspiration",
                "PlantPhotosynthesisLeaf","PlantTranspirationLeaf",
                "Export")
  TYPES = c(TYPES_SWB, "PlantRespiration", "PlantRespirationLeaf", "PlantCBalance", "PlantCBalanceLeaf",
            "PlantCstorageFast", "PlantCstorageSlow", "PlantSAgrowth", "PlantRelativeSAgrowth", "PlantSA",
            "PlantLAIlive","PlantLAIdead")

  type = match.arg(type,TYPES)  

  if(type %in% TYPES_SWB) {
    x$PlantLAI = x$PlantLAIlive
    plot.swb(x,type, bySpecies, yearAxis, xlim, ylim, xlab, ylab, ...)
  } else {
    dates = as.Date(rownames(x$DailyBalance))
    transpMode = x$control$transpirationMode
    DailyBalance = x$DailyBalance
    SoilWaterBalance = x$SoilWaterBalance
    nlayers = x$NumSoilLayers
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
    if(is.null(xlab)) xlab = ifelse(yearAxis,"Year", "Date")  
    if(type=="PlantRespiration") {
      df = x$PlantRespiration
      if(bySpecies) {
        df = t(apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T))
        cohortnames = colnames(df)
      } 
      if(is.null(ylab)) ylab = expression(paste("Plant respiration   ",(g*C%.%m^{-2})))
      if(is.null(ylim)) ylim = c(min(df),max(df))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    } 
    else if(type=="PlantRespirationLeaf") {
      df = x$PlantRespiration
      if(bySpecies) {
        m1 = apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T)
        lai1 = apply(x$PlantLAIlive,1,tapply, x$cohorts$Name, sum, na.rm=T)
        df = t(m1/lai1)
        cohortnames = colnames(df)
      } else {
        df = df/x$PlantLAIlive
        df[x$PlantLAIlive==0] = NA
      }
      if(is.null(ylab)) ylab = expression(paste("Plant respiration per leaf area   ",(g*C%.%m^{-2})))
      if(is.null(ylim)) ylim = c(min(df, na.rm=T),max(df, na.rm=T))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    } 
    else if(type=="PlantCBalance") {
      df = x$PlantPhotosynthesis-x$PlantRespiration
      if(bySpecies) {
        df = t(apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T))
        cohortnames = colnames(df)
      } 
      if(is.null(ylab)) ylab = expression(paste("Plant C balance   ",(g*C%.%m^{-2})))
      if(is.null(ylim)) ylim = c(min(df),max(df))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
      abline(h=0, col="gray")
    } 
    else if(type=="PlantCBalanceLeaf") {
      df = (x$PlantPhotosynthesis-x$PlantRespiration)
      if(bySpecies) {
        m1 = apply(df,1, tapply, x$cohorts$Name, sum, na.rm=T)
        lai1 = apply(x$PlantLAIlive,1,tapply, x$cohorts$Name, sum, na.rm=T)
        df = t(m1/lai1)
        cohortnames = colnames(df)
      } else {
        df = df/x$PlantLAIlive
        df[x$PlantLAIlive==0] = NA
      }
      if(is.null(ylab)) ylab = expression(paste("Plant C balance per leaf area  ",(g*C%.%m^{-2})))
      if(is.null(ylim)) ylim = c(min(df),max(df))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
      abline(h=0, col="gray")
      
    }
    else if(type=="PlantCstorageFast") {
      df = x$PlantCstorageFast
      if(bySpecies) {
        df = t(apply(df,1, tapply, x$cohorts$Name, mean, na.rm=T))
        cohortnames = colnames(df)
      } 
      if(is.null(ylab)) ylab = expression(paste("Fast C storage pool  [0-1]"))
      if(is.null(ylim)) ylim = c(min(df),max(df))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    } 
    else if(type=="PlantCstorageSlow") {
      df = x$PlantCstorageSlow
      if(bySpecies) {
        df = t(apply(df,1, tapply, x$cohorts$Name, mean, na.rm=T))
        cohortnames = colnames(df)
      } 
      if(is.null(ylab)) ylab = expression(paste("Slow C storage pool  [0-1]"))
      if(is.null(ylim)) ylim = c(min(df),max(df))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    } 
    else if(type=="PlantSAgrowth") {
      df = x$PlantSAgrowth
      if(bySpecies) {
        df = t(apply(df,1, tapply, x$cohorts$Name, mean, na.rm=T))
        cohortnames = colnames(df)
      } 
      if(is.null(ylab)) ylab = expression(paste("Sapwood area growth  ",(cm^2)))
      if(is.null(ylim)) ylim = c(min(df),max(df))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    } 
    else if(type=="PlantRelativeSAgrowth") {
      df = x$PlantSAgrowth/x$PlantSA
      if(bySpecies) {
        df = t(apply(df,1, tapply, x$cohorts$Name, mean, na.rm=T))
        cohortnames = colnames(df)
      } 
      if(is.null(ylab)) ylab = expression(paste("Relative sapwood area growth  ",(cm^2%.%cm^{-2})))
      if(is.null(ylim)) ylim = c(min(df),max(df))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    } 
    else if(type=="PlantSA") {
      df = x$PlantSA
      if(bySpecies) {
        df = t(apply(df,1, tapply, x$cohorts$Name, mean, na.rm=T))
        cohortnames = colnames(df)
      } 
      if(is.null(ylab)) ylab = expression(paste("Sapwood area  ",(cm^2)))
      if(is.null(ylim)) ylim = c(min(df),max(df))
      matplot(dates, df, ylim = ylim, 
              lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()     
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    }
    else if(type=="PlantLAIlive") {
      OM = x$PlantLAIlive
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
        cohortnames = colnames(OM)
      } 
      if(is.null(ylab)) ylab = expression(paste("Live Leaf Area Index   ",(m^{2}%.%m^{-2})))
      if(is.null(ylim)) ylim = c(0,max(OM, na.rm=TRUE))
      matplot(dates, OM, ylim = ylim, lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    } 
    else if(type=="PlantLAIdead") {
      OM = x$PlantLAIdead
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
        cohortnames = colnames(OM)
      } 
      if(is.null(ylab)) ylab = expression(paste("Dead Leaf Area Index   ",(m^{2}%.%m^{-2})))
      if(is.null(ylim)) ylim = c(0,max(OM, na.rm=TRUE))
      matplot(dates, OM, ylim = ylim, lty=1:length(cohortnames), col = 1:length(cohortnames),
              lwd=1, type="l", xlim=xlim,
              ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
      plotAxes()
      legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
             col = 1:length(cohortnames), bty="n")
    } 
    
  }
}
