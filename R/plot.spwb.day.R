plot.spwb.day<-function(x, type="PlantPsi", bySpecies = FALSE, xlab = NULL, ylab = NULL, ...) {
  if(!("EnergyBalance" %in% names(x))) stop("Plotting function available for transpirationMode = 'Complex' only.")
  EB = x$EnergyBalance
  Plants = x$Plants
  PlantsInst = x$PlantsInst
  TYPES = c("LeafPsi","RootPsi",
            "StemPLC","StemRWC", "LeafRWC",
            "PlantTranspiration","PlantPhotosynthesis","PlantAbsorbedSWR",
            "LeafTranspiration","LeafPhotosynthesis", "LeafAbsorbedSWR",
            "LeafVPD","LeafStomatalConductance", "LeafTemperature",
            "Temperature","CanopyEnergyBalance", "SoilEnergyBalance")
  type = match.arg(type,TYPES)  
  cohortnames = row.names(x$cohorts)
  timesteps = as.numeric(colnames(x$PlantsInst$PsiLeaf))
  if(is.null(xlab)) xlab = "Time step"
  if(type=="LeafPsi") {
    OM = PlantsInst$PsiLeaf
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM)
    } 
    if(is.null(ylab)) ylab = "Leaf water potential (MPa)"
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="RootPsi") {
    OM = PlantsInst$PsiRoot
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM)
    } 
    if(is.null(ylab)) ylab = "Root crown water potential (MPa)"
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="StemPLC") {
    OM = PlantsInst$PLCstem*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM)
    } 
    if(is.null(ylab)) ylab = "Stem percent of conductance loss (%)"
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="StemRWC") {
    OM = PlantsInst$RWCstem*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM)
    } 
    if(is.null(ylab)) ylab = "Stem symplasm relative water content (%)"
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="LeafRWC") {
    OM = PlantsInst$RWCleaf*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM)
    } 
    if(is.null(ylab)) ylab = "Leaf symplasm relative water content (%)"
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="PlantTranspiration") {
    OM = PlantsInst$E
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      cohortnames = rownames(OM)
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant transpiration   ",(L%.%m^{-2})))
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="LeafTranspiration") {
    OM = PlantsInst$E
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM)
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = expression(paste("Transpiration per leaf area   ",(L%.%m^{-2})))
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="PlantPhotosynthesis") {
    OM = PlantsInst$An
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      cohortnames = rownames(OM)
    } 
    if(is.null(ylab)) ylab = expression(paste("Net photosynthesis  ",(g*C%.%m^{-2})))
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="LeafPhotosynthesis") {
    OM = PlantsInst$An
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM)
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = expression(paste("Net photosynthesis per leaf area  ",(g*C%.%m^{-2})))
    matplot(timesteps, t(OM), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  else if(type=="PlantAbsorbedSWR") {
    old = par(mfrow=c(1,2), mar=c(5,5,3,1))
    OM_SL = PlantsInst$AbsRad$SWR_SL
    OM_SH = PlantsInst$AbsRad$SWR_SH
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM_SL)
    } 
    matplot(timesteps, t(OM_SL), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=expression(paste("Absorbed SWR sunlit    ",(W%.%m^{-2}))), 
            xlab=xlab, frame=FALSE, ...)
    matplot(timesteps, t(OM_SH), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=expression(paste("Absorbed SWR shade    ",(W%.%m^{-2}))), 
            xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
    par(old)
  }
  else if(type=="LeafAbsorbedSWR") {
    old = par(mfrow=c(1,2), mar=c(5,5,3,1))
    OM_SL = PlantsInst$AbsRad$SWR_SL
    OM_SH = PlantsInst$AbsRad$SWR_SH
    if(bySpecies) {
      m1 = apply(OM_SL,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(Plants$LAIsunlit,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = m1/lai1
      m1 = apply(OM_SH,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(Plants$LAIshade,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = m1/lai1
      cohortnames = rownames(OM_SL)
    } else {
      OM_SL = OM_SL/PlantsInst$LAIsunlit
      OM_SH = OM_SH/PlantsInst$LAIshade
    }
    matplot(timesteps, t(OM_SL), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=expression(paste("Absorbed SWR sunlit per leaf   ",(W%.%m^{-2}))), 
            xlab=xlab, frame=FALSE, ...)
    matplot(timesteps, t(OM_SH), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=expression(paste("Absorbed SWR shade per leaf   ",(W%.%m^{-2}))), 
            xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
    par(old)
  }
  else if(type=="LeafTemperature") {
    old = par(mfrow=c(1,2), mar=c(5,5,3,1))
    OM_SL = PlantsInst$Tempsunlit
    OM_SH = PlantsInst$Tempshade
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM_SL)
    } 
    matplot(timesteps, t(OM_SL), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab="Leaf temperature sunlit (degrees C)", xlab=xlab, frame=FALSE, ...)
    matplot(timesteps, t(OM_SH), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab="Leaf temperature shade (degrees C)", xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
    par(old)
  }
  else if(type=="LeafVPD") {
    old = par(mfrow=c(1,2), mar=c(5,5,3,1))
    OM_SL = PlantsInst$VPDsunlit
    OM_SH = PlantsInst$VPDshade
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM_SL)
    } 
    matplot(timesteps, t(OM_SL), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab="Vapour pressure deficit sunlit (kPa)", xlab=xlab, frame=FALSE, ...)
    matplot(timesteps, t(OM_SH), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab="Vapour pressure deficit shade (kPa)", xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
    par(old)
  }
  else if(type=="LeafStomatalConductance") {
    old = par(mfrow=c(1,2), mar=c(5,5,3,1))
    OM_SL = PlantsInst$GWsunlit
    OM_SH = PlantsInst$GWshade
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
      cohortnames = rownames(OM_SL)
    } 
    
    matplot(timesteps, t(OM_SL), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=expression(paste("Stomatal conductance sunlit   ", (mmol%.%m^{-2}%.%s^{-1}))), 
            xlab=xlab, frame=FALSE, ...)
    matplot(timesteps, t(OM_SH), lty=1:length(cohortnames), col = 1:length(cohortnames),
            lwd=1, type="l", ylab=expression(paste("Stomatal conductance shade   ", (mmol%.%m^{-2}%.%s^{-1}))), 
            xlab=xlab, frame=FALSE, ...)
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
    par(old)
  }
  else if(type=="Temperature") {
    
    if(is.null(ylab)) ylab = "Temperature (degrees C)"
    plot(timesteps, EB$Temperature$Tatm, 
         ylim = c(min(c(EB$Temperature$Tatm,EB$Temperature$Tcan,EB$Temperature$Tsoil.1)),
                  max(c(EB$Temperature$Tatm,EB$Temperature$Tcan,EB$Temperature$Tsoil.1))),
         lty=1, col = "black",
          type="l", ylab=ylab, xlab=xlab, frame=FALSE, ...)
    lines(timesteps, EB$Temperature$Tcan, col="red", lty=2)
    lines(timesteps, EB$Temperature$Tsoil.1, col="gray", lty=3)
    
    legend("topleft", legend = c("Above-canopy", "Inside canopy",
                                  "Soil"), 
           lty=1:3, 
           col = c("black", "red", "gray"), bty="n")
  } 
  else if(type=="CanopyEnergyBalance") {
    if(is.null(ylab)) ylab = expression(W%.%m^{-2})     
    mxmn = max(abs(EB$CanopyEnergyBalance))    
    plot(timesteps, EB$CanopyEnergyBalance$Ebalcan, type="n", 
         ylab=ylab, xlab=xlab, ylim = c(-mxmn,mxmn),
         frame=FALSE,...)
    abline(h=0, col="black", lwd=1.5)
    lines(timesteps, EB$CanopyEnergyBalance$Ebalcan, col="black",...)
    lines(timesteps, EB$CanopyEnergyBalance$SWRcanin, col="red",...)
    lines(timesteps, EB$CanopyEnergyBalance$LWRcanin, col="brown",...)
    lines(timesteps, -EB$CanopyEnergyBalance$LWRcanout, col="blue",...)
    lines(timesteps, EB$CanopyEnergyBalance$LWRsoilcan, col="orange",...)
    lines(timesteps, -EB$CanopyEnergyBalance$LEcan, col="green",...)
    lines(timesteps, -EB$CanopyEnergyBalance$Hcan, col="gray",...)
    lines(timesteps, -EB$SoilEnergyBalance$Hcansoil, col="dark gray",...)
    legend("topright", bty="n", col=c("red","brown","orange", "blue","green", "gray", "dark gray", "black"), lty=1,
           legend=c("SWR abs. from atm.","LWR abs. from atm.","LWR abs. from soil","LWR emmited", "Latent heat (L)",
                    "Convection can./atm.","Convection soil/can.", "Balance"),...)        
  }
  else if(type=="SoilEnergyBalance") {
    if(is.null(ylab)) ylab = expression(W%.%m^{-2})    
    mxmn = max(abs(EB$SoilEnergyBalance))    
    plot(timesteps, EB$SoilEnergyBalance$Ebalsoil, type="n", 
         ylab=ylab, xlab=xlab, ylim = c(-mxmn,mxmn),
         frame=FALSE,...)
    abline(h=0, col="black", lwd=1.5)
    lines(timesteps, EB$SoilEnergyBalance$Ebalsoil, col="black",...)
    lines(timesteps, EB$SoilEnergyBalance$SWRsoilin, col="red",...)
    lines(timesteps, EB$SoilEnergyBalance$LWRsoilin, col="brown",...)
    lines(timesteps, EB$CanopyEnergyBalance$LWRcanout, col="orange",...)
    lines(timesteps, -EB$SoilEnergyBalance$LWRsoilout, col="blue",...)
    lines(timesteps, EB$SoilEnergyBalance$Hcansoil, col="gray",...)
    legend("topright", bty="n", col=c("red","brown","orange", "blue", "gray", "black"), lty=1,
           legend=c("SWR abs. from atm.","LWR abs. from atm.", "LWR abs. from canopy","LWR emmited", "Convection soil/can.", "Balance"),...)        
  }
}
