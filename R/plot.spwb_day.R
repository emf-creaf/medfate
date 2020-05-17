plot.growth_day<-function(x, type="PlantTranspiration", bySpecies = FALSE, 
                          xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...) {
  GROWTH_TYPES = .getSubdailyGROWTHPlotTypes()
  PWB_TYPES = .getSubdailySPWBPlotTypes()
  type = match.arg(type,GROWTH_TYPES)
  if(type %in% PWB_TYPES) {
    plot.pwb_day(x, type, bySpecies, xlim, ylim, xlab, ylab,...)
  }
  CBInst = x$PlantCBInst
  if(type=="SugarLeaf") {
    OM = CBInst$SugarLeaf
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
}

plot.spwb_day<-function(x, type="PlantTranspiration", bySpecies = FALSE, 
                        xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...) {
  plot.pwb_day(x, type, bySpecies, xlim, ylim, xlab, ylab,...)
}


plot.pwb_day<-function(x, type="PlantTranspiration", bySpecies = FALSE, 
                       xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...) {
  if(!("EnergyBalance" %in% names(x))) stop("Plotting function available for transpirationMode = 'Sperry' only.")
  EB = x$EnergyBalance
  Plants = x$Plants
  PlantsInst = x$PlantsInst
  SunlitLeavesInst = x$SunlitLeavesInst
  ShadeLeavesInst = x$ShadeLeavesInst
  type = match.arg(type,.getSubdailySPWBPlotTypes())  
  cohortnames = row.names(x$cohorts)
  timesteps = as.numeric(colnames(x$PlantsInst$LeafPsi))
  if(is.null(xlab)) xlab = "Time step"
  if(type=="LeafPsiAverage") {
    OM = PlantsInst$LeafPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Leaf water potential (MPa)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemPsi") {
    OM = PlantsInst$StemPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Stem water potential (MPa)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafSympPsi") {
    OM = PlantsInst$LeafSympPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Leaf symplastic water potential (MPa)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemSympPsi") {
    OM = PlantsInst$StemSympPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Stem symplastic water potential (MPa)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="RootPsi") {
    OM = PlantsInst$RootPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Root crown water potential (MPa)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemPLC") {
    OM = PlantsInst$StemPLC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Stem percent of conductance loss (%)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemRWC") {
    OM = PlantsInst$StemRWC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Stem relative water content (%)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemSympRWC") {
    OM = PlantsInst$StemSympRWC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Stem symplasm relative water content (%)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafRWC") {
    OM = PlantsInst$LeafRWC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Leaf relative water content (%)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafSympRWC") {
    OM = PlantsInst$LeafSympRWC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = "Leaf symplasm relative water content (%)"
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="SoilPlantConductance") {
    OM = PlantsInst$dEdP
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
    } 
    if(is.null(ylab)) ylab = expression(paste("Soil-plant conductance ",(mmol%.%m^{-2}%.%MPa^{-1}%.%s^{-1})))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantWaterBalance") {
    OM = PlantsInst$PWB
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab =  expression(paste("Extraction - transpiration (",L%.%m^{-2},")"))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="WaterBalancePerLeaf") {
    OM = PlantsInst$PWB
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab =  expression(paste("Extraction - transpiration per leaf area (",L%.%m^{-2},")"))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantExtraction") {
    OM = x$ExtractionInst
    nlayers = nrow(OM)
    if(is.null(ylab)) ylab = expression(paste("Extraction from soil layers   ",(L%.%m^{-2})))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantTranspiration") {
    OM = PlantsInst$E
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant transpiration   ",(L%.%m^{-2})))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="TranspirationPerLeaf") {
    OM = PlantsInst$E
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = expression(paste("Transpiration per leaf area   ",(L%.%m^{-2})))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantGrossPhotosynthesis") {
    OM = PlantsInst$Ag
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant gross photosynthesis  ",(gC%.%m^{-2})))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantNetPhotosynthesis") {
    OM = PlantsInst$An
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant net photosynthesis  ",(gC%.%m^{-2})))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="GrossPhotosynthesisPerLeaf") {
    OM = PlantsInst$Ag
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = expression(paste("Gross photosynthesis per leaf area ",(gC%.%m^{-2})))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="NetPhotosynthesisPerLeaf") {
    OM = PlantsInst$An
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = expression(paste("Net photosynthesis per leaf area ",(gC%.%m^{-2})))
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafTranspiration") {
    OM_SL = SunlitLeavesInst$VPD*SunlitLeavesInst$GW
    OM_SH = ShadeLeavesInst$VPD*ShadeLeavesInst$GW
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=expression(paste("Leaf transpiration ",(mmol%.%m^{-2}%.%s^{-1})))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafNetPhotosynthesis") {
    OM_SL = SunlitLeavesInst$An
    OM_SH = ShadeLeavesInst$An
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=expression(paste("Leaf net photosynthesis ",(mu%.%mol%.%m^{-2}%.%s^{-1})))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafGrossPhotosynthesis") {
    OM_SL = SunlitLeavesInst$Ag
    OM_SH = ShadeLeavesInst$Ag
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=expression(paste("Leaf gross photosynthesis ",(mu%.%mol%.%m^{-2}%.%s^{-1})))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantAbsorbedSWR") {
    OM_SL = SunlitLeavesInst$Abs_SWR
    OM_SH = ShadeLeavesInst$Abs_SWR
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=expression(paste("Plant absorbed SWR ",(W%.%m^{-2})))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafAbsorbedSWR") {
    OM_SL = SunlitLeavesInst$Abs_SWR
    OM_SH = ShadeLeavesInst$Abs_SWR
    if(bySpecies) {
      m1 = apply(OM_SL,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$SunlitLeaves$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = m1/lai1
      m1 = apply(OM_SH,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$ShadeLeaves$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = m1/lai1
    } else {
      OM_SL = OM_SL/x$SunlitLeaves$LAI
      OM_SH = OM_SH/x$ShadeLeaves$LAI
    }
    if(is.null(ylab)) ylab=expression(paste("Absorbed SWR per leaf area ",(W%.%m^{-2})))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafAbsorbedLWR") {
    OM_SL = SunlitLeavesInst$Abs_LWR
    OM_SH = ShadeLeavesInst$Abs_LWR
    if(bySpecies) {
      m1 = apply(OM_SL,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$SunlitLeaves$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = m1/lai1
      m1 = apply(OM_SH,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$ShadeLeaves$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = m1/lai1
    } else {
      OM_SL = OM_SL/x$SunlitLeaves$LAI
      OM_SH = OM_SH/x$ShadeLeaves$LAI
    }
    if(is.null(ylab)) ylab=expression(paste("Absorbed LWR per leaf area ",(W%.%m^{-2})))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafPsi") {
    OM_SL = SunlitLeavesInst$Psi
    OM_SH = ShadeLeavesInst$Psi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab="Leaf water potential (-MPa)"
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafTemperature") {
    OM_SL = SunlitLeavesInst$Temp
    OM_SH = ShadeLeavesInst$Temp
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab="Leaf temperature (degrees C)"
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafVPD") {
    OM_SL = SunlitLeavesInst$VPD
    OM_SH = ShadeLeavesInst$VPD
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab="Leaf vapour pressure deficit (kPa)"
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafStomatalConductance") {
    OM_SL = SunlitLeavesInst$GW
    OM_SH = ShadeLeavesInst$GW
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=expression(paste("Stomatal conductance ", (mmol%.%m^{-2}%.%s^{-1})))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafIntrinsicWUE") {
    OM_SL = SunlitLeavesInst$An/SunlitLeavesInst$GW
    OM_SH = ShadeLeavesInst$An/ShadeLeavesInst$GW
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=expression(paste("iWUE  ", (mu%.%mol%.%mol^{-1})))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafCi") {
    OM_SL = SunlitLeavesInst$Ci
    OM_SH = ShadeLeavesInst$Ci
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=expression(paste("Intercellular CO2 concentration ", (ppm)))
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="Temperature") {
    
    if(is.null(ylab)) ylab = "Temperature (degrees C)"
    df = data.frame(row.names=timesteps)
    df[["Above-canopy"]] = EB$Temperature$Tatm
    df[["Inside-canopy"]] = EB$Temperature$Tcan
    df[["Soil"]] = EB$Temperature$Tsoil.1
    g<-.multiple_subday_dynamics(as.matrix(df), ylab=ylab, ylim = ylim)
    return(g)
  } 
  else if(type=="CanopyEnergyBalance") {
    if(is.null(ylab)) ylab = expression(W%.%m^{-2})
    df = data.frame(row.names=timesteps)
    df[["Balance"]] = EB$CanopyEnergyBalance$Ebalcan
    df[["SWR abs. from atm."]] = EB$CanopyEnergyBalance$SWRcanin 
    df[["LWR abs. from atm."]] = EB$CanopyEnergyBalance$LWRcanin
    df[["LWR abs. from soil"]] = EB$CanopyEnergyBalance$LWRsoilcan
    df[["LWR emmited"]] = -EB$CanopyEnergyBalance$LWRcanout
    df[["Latent heat"]] = -EB$CanopyEnergyBalance$LEcan
    df[["Convection can./atm."]] = -EB$CanopyEnergyBalance$Hcan
    df[["Convection soil/can."]] = -EB$SoilEnergyBalance$Hcansoil
    g<-.multiple_subday_dynamics(as.matrix(df), ylab=ylab, ylim = ylim)
    # cols = c("red","brown","orange", "blue","green", "gray", "dark gray", "black")
    # names(cols) = c("SWR abs. from atm.","LWR abs. from atm.","LWR abs. from soil","LWR emmited", "Latent heat (L)",
    #                 "Convection can./atm.","Convection soil/can.", "Balance")
    return(g)
  }
  else if(type=="SoilEnergyBalance") {
    if(is.null(ylab)) ylab = expression(W%.%m^{-2})    
    df = data.frame(row.names=timesteps)
    df[["Balance"]] = EB$SoilEnergyBalance$Ebalsoil
    df[["SWR abs. from atm."]] = EB$SoilEnergyBalance$SWRsoilin
    df[["LWR abs. from atm."]] = EB$SoilEnergyBalance$LWRsoilin
    df[["LWR abs. from canopy"]] = EB$CanopyEnergyBalance$LWRcanout
    df[["LWR emmited"]] = -EB$SoilEnergyBalance$LWRsoilout
    df[["Convection soil/can."]] = EB$SoilEnergyBalance$Hcansoil
    df[["Latent heat"]] = -EB$SoilEnergyBalance$LEsoil
    g<-.multiple_subday_dynamics(as.matrix(df), ylab=ylab, ylim = ylim)
    
    # legend("topright", bty="n", col=c("red","brown","orange", "blue", "gray", "black"), lty=1,
    #        legend=c("SWR abs. from atm.","LWR abs. from atm.", "LWR abs. from canopy","LWR emmited", 
    #                   "Convection soil/can.", "Balance"),...)        
    return(g)
  }
}
