plot.swb<-function(x, type="PET_Precipitation", yearAxis=FALSE, xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,...) {
  dates = as.Date(rownames(x$DailyBalance))
  DailyBalance = x$DailyBalance
  SoilWaterBalance = x$SoilWaterBalance
  nlayers = x$NumSoilLayers
  TYPES = c("PET_Precipitation","PET_NetPrec","ET","Psi","Theta","Vol", "LAI",
            "PlantStress", "PlantPsi","PlantPhotosynthesis","PlantTranspiration",
            "AirTemperature","SoilTemperature", "CanopyTemperature","Export")
  type = match.arg(type,TYPES)  
  numDays = length(dates)
  numYears = round(numDays/365)
  firstYear=as.numeric(format(dates[1],"%Y"))
  
  plotAxes<-function(){
    if(!yearAxis) axis.Date(1, dates)
    else {
      axis(1, at = (0:numYears)*365, labels=FALSE)
      axis(1, at = -182+365*(1:numYears), tick = FALSE, line=FALSE, labels=firstYear:(firstYear+numYears-1))
    }
    axis(2)    
  }
  mnp = max(DailyBalance$Precipitation)
  if(is.null(xlab)) xlab = ifelse(yearAxis,"Year", "Date")  
  if(type=="PET_Precipitation") {
    if(is.null(ylab)) ylab = "mm water"
    if(!is.null(xlim)) span = xlim[1]:xlim[2]
    else span = 1:numDays
    if(is.null(ylim)) ylim = c(0,mnp)
    barplot(DailyBalance$Precipitation[span], ylim=ylim, col="black",space=0, ylab=ylab, 
            xlab=xlab, axes=FALSE)
    plotAxes()
    lines(1:length(span), DailyBalance$PET[span], col="gray")    
    legend("topleft", bty="n", col=c("black","gray"),lty=c(1,1), lwd=2,
           legend=c("Precipitation","PET"))
    
  } else if(type=="PET_NetPrec") {
    if(is.null(ylab)) ylab = "mm water"    
    if(!is.null(xlim)) span = xlim[1]:xlim[2]
    else span = 1:numDays
    if(is.null(ylim)) ylim = c(0,mnp)
    barplot(DailyBalance$NetPrec[span], ylim=ylim, col="black",space=0, ylab=ylab, 
            xlab=xlab, axes=FALSE)
    plotAxes()
    lines(1:length(span), DailyBalance$PET[span], col="gray")    
    legend("topleft", bty="n", col=c("black","gray"),lty=c(1,1), lwd=2,
           legend=c("NetPrec","PET"))        
  } else if(type=="ET") {
    if(is.null(ylab)) ylab = "mm water"
    if(is.null(ylim)) ylim = c(0,max(DailyBalance$Etot))
    plot(dates, DailyBalance$Etot, ylim=ylim, type="l", ylab=ylab, 
         xlab=xlab, xlim=xlim,frame=FALSE, col="black", axes=FALSE, lwd=2)
    plotAxes()
    lines(dates, DailyBalance$Eplanttot, col="gray", lty=2, lwd=1.5)
    lines(dates, DailyBalance$Esoil, col="black", lty=3, lwd=1.5)
    legend("topleft", bty="n", col=c("black","gray","black"),lty=c(1,2,3), lwd=c(2,1.5,1.5),
           legend=c("Total evapotranspiration","Plant transpiration","Bare soil evaporation"))
  } else if(type=="LAI") {
    if(is.null(ylab)) ylab = "LAI (m2/m2)"
    if(is.null(ylim)) ylim = c(0,max(DailyBalance$LAIcell))
    plot(dates, DailyBalance$LAIcell, ylim=ylim, type="l", ylab=ylab, 
         xlab=xlab, xlim=xlim,frame=FALSE, col="black", axes=FALSE, lwd=1)
    plotAxes()
  } else if(type=="Psi") {
    PsiM = SoilWaterBalance[,paste("psi",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Soil water potential (MPa)"    
    if(is.null(ylim)) ylim =c(min(PsiM),0)
    matplot(dates, PsiM, lwd=1.5,
         ylim=ylim, type="l", ylab=ylab, xlab=xlab, xlim=xlim,
         frame=FALSE, lty=c(1,2,3,4,5), col="black",axes=FALSE)
    plotAxes()
    legend("bottomleft", bty="n", lty=c(1,2,3,4,5), col="black", lwd=1.5,
           legend=paste("Layer", 1:nlayers))
    
  } else if(type=="Theta") {
    WM = SoilWaterBalance[,paste("W",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "% field capacity"
    if(is.null(ylim)) ylim = c(0,100)
    matplot(dates, WM*100, lwd=1.5,
            ylim=ylim, type="l", ylab=ylab, xlab=xlab, xlim=xlim,
            frame=FALSE, lty=c(1,2,3,4,5), col="black",axes=FALSE)
    plotAxes()
    legend("bottomleft", bty="n", col="black",lty=1:nlayers, lwd=1.5,
           legend=paste("Layer", 1:nlayers))
  } else if(type=="Vol") {
    MLM= SoilWaterBalance[,paste("ML",1:nlayers,sep=".")]
    MLTot = SoilWaterBalance$MLTot
    if(is.null(ylim)) ylim =c(0,max(MLTot)*1.3)
    if(is.null(ylab)) ylab = "mm soil water"
    plot(dates, MLTot, ylim=ylim, lwd=2, type="l",xlim=xlim,
         ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    matlines(dates, MLM, lty = 1:nlayers, col="black", lwd=1.5)
    legend("topleft", bty="n", col="black",lty=c(1,1:nlayers), lwd=c(2,rep(1.5,5)),
           legend=c("Total",paste("Layer", 1:nlayers)))
    
  } else if(type=="PlantStress") {
    if(is.null(ylab)) ylab = "Drought stress [0-1]"
    if(is.null(ylim)) ylim = c(0,1)
    matplot(dates, x$PlantStress, ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()        
  } else if(type=="PlantPsi") {
    if(is.null(ylab)) ylab = "Plant water potential (MPa)"
    if(is.null(ylim)) ylim = c(min(x$PlantPsi),0)
    matplot(dates, x$PlantPsi, ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()      
  } else if(type=="PlantTranspiration") {
    if(is.null(ylab)) ylab = "Plant transpiration (mm)"
    if(is.null(ylim)) ylim = c(0,max(x$PlantTranspiration))
    matplot(dates, x$PlantTranspiration, ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()      
  } else if(type=="PlantPhotosynthesis") {
    if(is.null(ylab)) ylab = "Plant photosynthesis (gC/m2)"
    if(is.null(ylim)) ylim = c(min(x$PlantPhotosynthesis),max(x$PlantPhotosynthesis))
    matplot(dates, x$PlantPhotosynthesis, ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()     
  } else if(type=="AirTemperature") {
    if(is.null(ylab)) ylab = "Above-canopy temperature (Celsius)"
    if(is.null(ylim)) ylim = c(min(x$EnergyBalance$Tatm_min),max(x$EnergyBalance$Tatm_max))
    plot(dates, x$EnergyBalance$Tatm_mean, ylim = ylim, lwd=1, type="l", xlim=xlim,
         ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    lines(dates, x$EnergyBalance$Tatm_min, col="blue")
    lines(dates, x$EnergyBalance$Tatm_max, col="red")
    plotAxes()   
  } else if(type=="CanopyTemperature") {
    if(is.null(ylab)) ylab = "Canopy temperature (Celsius)"
    if(is.null(ylim)) ylim = c(min(x$EnergyBalance$Tcan_min),max(x$EnergyBalance$Tcan_max))
    plot(dates, x$EnergyBalance$Tcan_mean, ylim = ylim, lwd=1, type="l", xlim=xlim,
         ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    lines(dates, x$EnergyBalance$Tcan_min, col="blue")
    lines(dates, x$EnergyBalance$Tcan_max, col="red")
    plotAxes()   
  } else if(type=="SoilTemperature") {
    if(is.null(ylab)) ylab = "Soil temperature (Celsius)"
    if(is.null(ylim)) ylim = c(min(x$EnergyBalance$Tsoil_min),max(x$EnergyBalance$Tsoil_max))
    plot(dates, x$EnergyBalance$Tsoil_mean, ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    lines(dates, x$EnergyBalance$Tsoil_min, col="blue")
    lines(dates, x$EnergyBalance$Tsoil_max, col="red")
    plotAxes()   
  } else if(type=="Export") {
    if(is.null(ylab)) ylab = "mm water"    
    mnp = max(DailyBalance$DeepDrainage+DailyBalance$Runoff)    
    if(is.null(ylim)) ylim = c(0,mnp)
    plot(dates, DailyBalance$DeepDrainage+DailyBalance$Runoff, ylim=ylim, col="black", type="l", 
         ylab=ylab, xlab=xlab, xlim=xlim,
         frame=FALSE, axes=FALSE)
    lines(dates, DailyBalance$DeepDrainage, col="blue")
    lines(dates, DailyBalance$Runoff, col="red")
    plotAxes()
    legend("topright", bty="n", col=c("black","blue","red"),lty=c(1,1,1), lwd=c(1.5,1,1),
           legend=c("DD+R","Deep drainage (DD)","Runoff (R)"))        
  }
}
