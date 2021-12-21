.getDailyPWBPlotTypes<-function(transpirationMode = "Granier") {
  if(transpirationMode=="Granier") {
    TYPES = c("SoilPsi","SoilTheta", "SoilRWC", "LAI",
              "PlantExtraction","PlantLAI", 
              "PlantStress", "PlantPsi","StemPLC",
              "PlantGrossPhotosynthesis", "PlantTranspiration",
              "GrossPhotosynthesisPerLeaf","TranspirationPerLeaf")
  } else {
    TYPES = c("SoilPsi","SoilTheta", "SoilRWC", "LAI",
              "PlantExtraction","HydraulicRedistribution",
              "PlantLAI",
              "SoilPlantConductance","PlantStress", 
              "PlantNetPhotosynthesis", "PlantGrossPhotosynthesis", "PlantTranspiration","PlantWUE",
              "NetPhotosynthesisPerLeaf","GrossPhotosynthesisPerLeaf","TranspirationPerLeaf",
              "TempMin_SL", "TempMin_SH", "TempMax_SL","TempMax_SH",
              "GSWMin_SL", "GSWMin_SH", "GSWMax_SL", "GSWMax_SH", "LeafPsiRange",
              "LeafPsiMin", "LeafPsiMax", "LeafPsiMin_SL", "LeafPsiMax_SL", "LeafPsiMin_SH", "LeafPsiMax_SH",
              "StemPsi","RootPsi","StemPLC", "StemRWC", "LeafRWC", "StemSympRWC", "LeafSympRWC", 
              "PlantWaterBalance",
              "PlantAbsorbedSWR", "AbsorbedSWRPerLeaf",
              "PlantNetLWR", "NetLWRPerLeaf",
              "Temperature","TemperatureRange", "AirTemperature","SoilTemperature", "CanopyTemperature",
              "CanopyEnergyBalance", "SoilEnergyBalance")
  }
  return(TYPES)
}
.getDailySPWBPlotTypes<-function(transpirationMode = "Granier") {
  TYPES = c(.getDailyPWBPlotTypes(transpirationMode),
            "PET_Precipitation","PET_NetRain","Snow","Evapotranspiration",
            "SoilVol", 
            "Export", "WTD")
  return(TYPES)
}

.getSubdailySPWBPlotTypes<-function(){
  TYPES = c("LeafPsi","LeafPsiAverage","RootPsi", "StemPsi", 
            "LeafSympPsi", "StemSympPsi",
            "StemPLC","StemRWC", "LeafRWC","StemSympRWC", "LeafSympRWC",
            "SoilPlantConductance",
            "PlantExtraction","PlantTranspiration", "TranspirationPerLeaf",
            "PlantGrossPhotosynthesis","GrossPhotosynthesisPerLeaf","PlantNetPhotosynthesis","NetPhotosynthesisPerLeaf", 
            "PlantAbsorbedSWR",
            "LeafTranspiration","LeafNetPhotosynthesis", "LeafGrossPhotosynthesis", 
            "LeafAbsorbedSWR","LeafNetLWR",
            "LeafCi", "LeafIntrinsicWUE",
            "LeafVPD","LeafStomatalConductance", "LeafTemperature",
            "Temperature","CanopyEnergyBalance", "SoilEnergyBalance", 
            "PlantWaterBalance", "WaterBalancePerLeaf")
  return(TYPES)
}
.getUniqueDailyGROWTHPlotTypes<-function(transpirationMode = "Granier"){
  TYPES = c("GrossPhotosynthesis","MaintenanceRespiration","GrowthCosts", "CarbonBalance",
            "SugarTransport", "LeafPI0", "StemPI0",
            "SugarLeaf", "SugarSapwood", "StarchLeaf", "StarchSapwood","SugarTransport", "RootExudation",
            "SapwoodArea", "LeafArea", 
            "SapwoodBiomass", "LeafBiomass", "FineRootBiomass",
            "LabileBiomass", "TotalLivingBiomass",
            "SAgrowth", "LAgrowth", 
            "HuberValue")
  if(transpirationMode=="Sperry") {
    TYPES = c(TYPES, "FRAgrowth", "FineRootArea","RootAreaLeafArea")
  }
  return(TYPES)
}
.getDailyGROWTHPlotTypes<-function(transpirationMode = "Granier"){
  TYPES = c(.getUniqueDailyGROWTHPlotTypes(transpirationMode),
            .getDailySPWBPlotTypes(transpirationMode))
  return(TYPES)
}
.getSubdailyGROWTHPlotTypes<-function(){
  TYPES = c("GrossPhotosynthesis","MaintenanceRespiration","GrowthCosts", "CarbonBalance",
            "SugarTransport",
            "SugarLeaf", "SugarSapwood", "StarchLeaf", "StarchSapwood","SugarTransport",
            .getSubdailySPWBPlotTypes())
  return(TYPES)
}
.getYLab<-function(type) {
  ylab="Unknown"
  if(type=="PlantTranspiration") ylab = expression(paste("Plant transpiration ",(L%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantPhotosynthesis") ylab = expression(paste("Plant photosynthesis ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantGrossPhotosynthesis") ylab = expression(paste("Plant gross photosynthesis ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantNetPhotosynthesis") ylab = expression(paste("Plant net photosynthesis ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantWUE") ylab = expression(paste("Plant water use efficiency ",(g*C%.%L^{-1})))
  else if(type=="PlantAbsorbedSWR") ylab = expression(paste("Plant absorbed SWR ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantNetLWR") ylab = expression(paste("Plant net LWR ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantWaterBalance") ylab = expression(paste("Plant water balance   ",(L%.%m^{-2})))
  else if(type=="TranspirationPerLeaf") ylab = expression(paste("Transpiration per leaf area ",(L%.%m^{-2}%.%d^{-1})))
  else if(type=="PhotosynthesisPerLeaf") ylab = expression(paste("Photosynthesis per leaf area ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="GrossPhotosynthesisPerLeaf") ylab = expression(paste("Gross photosynthesis per leaf area ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="NetPhotosynthesisPerLeaf") ylab = expression(paste("Net photosynthesis per leaf area ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="AbsorbedSWRPerLeaf") ylab = expression(paste("Absorbed SWR per leaf area ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="NetLWRPerLeaf") ylab = expression(paste("Net LWR per leaf area ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="GrossPhotosynthesis") ylab=expression(paste("Gross photosynthesis ", (gGluc%.%gdry^{-1}%.%d^{-1})))
  else if(type=="MaintenanceRespiration") ylab=expression(paste("Maintenance respiration ", (gGluc%.%gdry^{-1}%.%d^{-1})))
  else if(type=="GrowthCosts") ylab=expression(paste("Growth costs ", (gGluc%.%gdry^{-1}%.%d^{-1})))
  else if(type=="CarbonBalance") ylab=expression(paste("Carbon balance ", (gGluc%.%gdry^{-1}%.%d^{-1})))
  else if(type=="SugarLeaf") ylab=expression(paste("Leaf sugar concentration  ", (mol%.%L^{-1})))
  else if(type=="StarchLeaf") ylab=expression(paste("Leaf starch concentration  ", (mol%.%L^{-1})))
  else if(type=="SugarSapwood") ylab=expression(paste("Sapwood sugar concentration  ", (mol%.%L^{-1})))
  else if(type=="StarchSapwood") ylab=expression(paste("Sapwood starch concentration  ", (mol%.%L^{-1})))
  else if(type=="SugarTransport") ylab=expression(paste("Floem sugar transport rate ", (mmol%.%s^{-1})))
  else if(type=="RootExudation") ylab=expression(paste("Root exudation ", (gGluc%.%gdry^{-1})))
  else if(type=="SapwoodArea")  ylab = expression(paste("Sapwood area  ",(cm^2)))
  else if(type=="LeafArea")  ylab = expression(paste("Leaf area  ",(m^2)))
  else if(type=="FineRootArea")  ylab = expression(paste("Fine root area  ",(m^2)))
  else if(type=="SapwoodBiomass")  ylab = expression(paste("Sapwood biomass  ", (gdry%.%ind^{-1})))
  else if(type=="LeafBiomass")  ylab = expression(paste("Leaf biomass  ", (gdry%.%ind^{-1})))
  else if(type=="FineRootBiomass")  ylab = expression(paste("Fine root biomass  ", (gdry%.%ind^{-1})))
  else if(type=="LabileBiomass")  ylab = expression(paste("Labile C biomass  ", (gdry%.%ind^{-1})))
  else if(type=="TotalLivingBiomass")  ylab = expression(paste("Total living biomass  ", (gdry%.%ind^{-1})))
  else if(type=="HuberValue")  ylab = expression(paste("Huber value  ",(cm^2 %.% m^{-2})))
  else if(type=="RootAreaLeafArea")  ylab = expression(paste("Root area / Leaf area  ",(m^2 %.% m^{-2})))
  else if(type=="SAgrowth") ylab = expression(paste("Sapwood area growth rate ",(cm^2 %.% cm^{-2} %.% d^{-1})))
  else if(type=="LAgrowth") ylab = expression(paste("Leaf area growth rate ",(m^2 %.% cm^{-2} %.% d^{-1})))
  else if(type=="FRAgrowth") ylab = expression(paste("Fine root area growth rate ",(m^2 %.% cm^{-2} %.% d^{-1})))
  else if(type=="LeafPI0")  ylab = expression(paste("Leaf osmotic potential at full turgor  ",(MPa)))
  else if(type=="StemPI0")  ylab = expression(paste("Stem osmotic potential at full turgor  ",(MPa)))
  else if(type=="StemPLC") ylab = "Percent loss conductance in stem [%]"
  else if(type=="StemRWC") ylab = "Relative water content in stem [%]"
  else if(type=="StemSympRWC") ylab = "Relative water content in stem symplasm [%]"
  else if(type=="LeafRWC") ylab = "Relative water content in leaf [%]"
  else if(type=="LeafSympRWC") ylab = "Relative water content in leaf symplasm [%]"
  else if(type=="PlantPsi") ylab = "Plant water potential (MPa)"
  else if(type=="PlantStress") ylab = "Drought stress [0-1]"
  else if(type=="StemPsi") ylab = "Midday stem water potential (MPa)"
  else if(type=="RootPsi") ylab = "Midday root crown water potential (MPa)"
  else if(type=="SoilPlantConductance") ylab = expression(paste("Average soil-plant conductance ",(mmol%.%m^{-2}%.%s^{-1})))
  else if(type=="LeafPsiMin") ylab = "Minimum (midday) leaf water potential (MPa)"
  else if(type=="LeafPsiMax") ylab = "Maximum (predawn) leaf water potential (MPa)"
  else if(type=="LeafPsiMin_SL") ylab = "Minimum (midday) sunlit leaf water potential (MPa)"
  else if(type=="LeafPsiMax_SL") ylab = "Maximum (predawn) sunlit leaf water potential (MPa)"
  else if(type=="LeafPsiMin_SH") ylab = "Minimum (midday) shade leaf water potential (MPa)"
  else if(type=="LeafPsiMax_SH") ylab = "Maximum (predawn) shade leaf water potential (MPa)"
  else if(type=="GSWMin_SH") ylab = expression(paste("Minimum shade leaf stomatal conductance ",(mol%.%m^{-2}%.%s^{-1})))
  else if(type=="GSWMax_SH") ylab = expression(paste("Maximum shade leaf stomatal conductance ",(mol%.%m^{-2}%.%s^{-1})))
  else if(type=="GSWMin_SL") ylab = expression(paste("Minimum sunlit leaf stomatal conductance ",(mol%.%m^{-2}%.%s^{-1})))
  else if(type=="GSWMax_SL") ylab = expression(paste("Maximum sunlit leaf stomatal conductance ",(mol%.%m^{-2}%.%s^{-1})))
  else if(type=="TempMin_SH") ylab = expression(paste("Minimum shade leaf temperature (Celsius)"))
  else if(type=="TempMax_SH") ylab = expression(paste("Maximum shade leaf temperature (Celsius)"))
  else if(type=="TempMin_SL") ylab = expression(paste("Minimum sunlit leaf temperature (Celsius)"))
  else if(type=="TempMax_SL") ylab = expression(paste("Maximum sunlit leaf temperature (Celsius)"))
  else if(type=="StandBasalArea") ylab = expression(paste("Stand basal area ",(m^{-2}%.%ha^{-1})))
  else if(type=="StandLAI") ylab = expression(paste("Stand leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="StandDensity") ylab = expression(paste("Stand tree density ",(ind%.%ha^{-1})))
  else if(type=="SpeciesBasalArea") ylab = expression(paste("Species basal area ",(m^{-2}%.%ha^{-1})))
  else if(type=="SpeciesLAI") ylab = expression(paste("Species leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="SpeciesDensity") ylab = expression(paste("Species tree density ",(ind%.%ha^{-1})))
  else if(type=="CohortBasalArea") ylab = expression(paste("Cohort basal area ",(m^{-2}%.%ha^{-1})))
  else if(type=="CohortLAI") ylab = expression(paste("Cohort leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="CohortDensity") ylab = expression(paste("Cohort tree density ",(ind%.%ha^{-1})))
  return(ylab)
}

.averageBySpecies<-function(OM, spnames) {
  if(ncol(OM)>1) OM = t(apply(OM,1, tapply, spnames, mean, na.rm=T))
  else colnames(OM) = spnames[1]
  return(OM)
}
.sumBySpecies<-function(OM, spnames) {
  if(ncol(OM)>1) OM = t(apply(OM,1, tapply, spnames, sum, na.rm=T))
  else colnames(OM) = spnames[1]
  return(OM)
}
.averageByLAISpecies<-function(OM, PlantsLAI, spnames) {
  if(ncol(OM)>1) {
    lai1 = t(apply(PlantsLAI,1, tapply, spnames, sum))
    m1 = t(apply(PlantsLAI * OM,1, tapply, spnames, sum))
    OM = m1/lai1
    OM[lai1==0] = NA
  }
  else colnames(OM) = spnames[1]
  return(OM)
}
.temporalSummary<-function(OM, summary.freq, FUN = mean, ...) {
  varnames = colnames(OM)
  date.factor = cut(as.Date(rownames(OM)), breaks=summary.freq)
  df = data.frame(row.names = as.Date(as.character(levels(date.factor))))
  for(i in 1:length(varnames)) {
    df[[varnames[i]]] = tapply(OM[,i],INDEX=date.factor, FUN=FUN, ...)
  }
  return(as.matrix(df))
}


.plot_wb<-function(WaterBalance, Soil, input_soil, type,  
                   dates = NULL, 
                   xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                   summary.freq = NULL, ...) {
  WaterBalance = as.data.frame(WaterBalance)
  Soil = as.data.frame(Soil)
  nlayers = length(input_soil$W)
  if(type=="PET_Precipitation") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2}) 
    df = data.frame(row.names=row.names(WaterBalance))
    df[["PET"]] = WaterBalance$PET
    df[["Precipitation"]] = WaterBalance$Precipitation
    df[["Snow"]] = WaterBalance$Snow
    df[["Date"]] = as.Date(row.names(WaterBalance))
    if(!is.null(dates)) df = df[df$Date %in% dates,]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(df$Date), breaks=summary.freq)
      df = data.frame(Date = as.Date(as.character(levels(date.factor))),
                      Precipitation = tapply(df$Precipitation,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      Snow = tapply(df$Snow,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      PET = tapply(df$PET,INDEX=date.factor, FUN=sum, na.rm=TRUE))
    }
    g<-ggplot(df)+
      geom_area(aes_string(x="Date", y="Precipitation", fill='"Precipitation"'))+
      geom_area(aes_string(x="Date", y="Snow", fill='"Snow"'))+
      geom_path(aes_string(x="Date", y="PET", col='"PET"'))+
      scale_fill_manual(name="", values=c("Precipitation"="black", "Snow"="red"))+
      scale_color_manual(name="", values=c("PET"="gray"))+
      ylab(ylab)+ xlab(xlab)+
      theme_bw()
    return(g)
  } 
  else if(type=="PET_NetRain") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2}) 
    df = data.frame(row.names=row.names(WaterBalance))
    df[["PET"]] = WaterBalance$PET
    df[["NetRain"]] = WaterBalance$NetRain
    df[["Date"]] = as.Date(row.names(WaterBalance))
    if(!is.null(dates)) df = df[df$Date %in% dates,]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(df$Date), breaks=summary.freq)
      df = data.frame(Date = as.Date(as.character(levels(date.factor))),
                      NetRain = tapply(df$NetRain,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      PET = tapply(df$PET,INDEX=date.factor, FUN=sum, na.rm=TRUE))
    }
    g<-ggplot(df)+
      geom_area(aes_string(x="Date", y="NetRain", fill='"NetRain"'))+
      geom_path(aes_string(x="Date", y="PET", col='"PET"'))+
      scale_fill_manual(name="", values=c("NetRain"="black"))+
      scale_color_manual(name="", values=c("PET"="gray"))+
      ylab(ylab)+xlab(xlab)+
      theme_bw()
    return(g)
  } 
  else if(type=="Evapotranspiration") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})
    df = data.frame(row.names=row.names(WaterBalance))
    df[["Total evapotranspiration"]] = WaterBalance$Evapotranspiration
    df[["Interception evaporation"]] = WaterBalance$Interception
    df[["Plant transpiration"]] = WaterBalance$Transpiration
    df[["Bare soil evaporation"]] = WaterBalance$SoilEvaporation
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(row.names(df)), breaks=summary.freq)
      df = data.frame(row.names = as.Date(as.character(levels(date.factor))),
                      "Total evapotranspiration" = tapply(df[["Total evapotranspiration"]],INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      "Interception evaporation" = tapply(df[["Interception evaporation"]],INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      "Plant transpiration" = tapply(df[["Plant transpiration"]],INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      "Bare soil evaporation" = tapply(df[["Bare soil evaporation"]],INDEX=date.factor, FUN=sum, na.rm=TRUE))
    }
    return(.multiple_dynamics(as.matrix(df), ylab=ylab, xlab=xlab, ylim = ylim))
  } 
  else if(type=="Snow") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})  
    df = data.frame(row.names=row.names(WaterBalance))
    df[["Snow"]] = WaterBalance$Snow
    df[["Snowpack"]] = Soil$SWE
    df[["Date"]] = as.Date(row.names(WaterBalance))
    if(!is.null(dates)) df = df[df$Date %in% dates,]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(df$Date), breaks=summary.freq)
      df = data.frame(Date = as.Date(as.character(levels(date.factor))),
                      Snow = tapply(df$Snow,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      Snowpack = tapply(df$Snowpack,INDEX=date.factor, FUN=mean, na.rm=TRUE))
    }
    g<-ggplot(df)+
      geom_area(aes_string(x="Date", y="Snow", fill='"Snow"'))+
      geom_path(aes_string(x="Date", y="Snowpack", col='"Snowpack"'))+
      scale_fill_manual(name="", values=c("Snow"="black"))+
      scale_color_manual(name="", values=c("Snowpack"="gray"))+
      ylab(ylab)+xlab(xlab)+
      theme_bw()
    return(g)
  } 
  else if(type=="WTD") {
    if(is.null(ylab)) ylab = expression(paste("Water table depth  (mm)"))
    xv = Soil$WTD
    names(xv) = row.names(Soil)
    if(!is.null(dates)) xv = xv[names(xv) %in% as.character(dates)]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(names(xv)), breaks=summary.freq)
      xv = tapply(xv,INDEX=date.factor, FUN=mean, na.rm=TRUE)
      names(xv) = as.character(levels(date.factor))
    }      
    return(.single_dynamics(xv, ylab = ylab, ylim = ylim))
  } 
  else if(type=="Export") {
    if(is.null(ylab)) ylab =  expression(L%.%m^{-2})    
    df = data.frame(row.names=row.names(WaterBalance))
    df[["Export"]] = WaterBalance$DeepDrainage + WaterBalance$Runoff
    df[["DeepDrainage"]] = WaterBalance$DeepDrainage
    df[["Runoff"]] = WaterBalance$Runoff 
    df[["Date"]]= as.Date(row.names(WaterBalance))
    if(!is.null(dates)) df = df[df$Date %in% dates,]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(df$Date), breaks=summary.freq)
      df = data.frame(Date = as.Date(as.character(levels(date.factor))),
                      Export = tapply(df$Export,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      DeepDrainage = tapply(df$DeepDrainage,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      Runoff = tapply(df$Runoff,INDEX=date.factor, FUN=sum, na.rm=TRUE))
    }
    g<-ggplot(df)+
      geom_line(aes_string(x="Date", y="Export", col='"Export"'))+
      geom_line(aes_string(x="Date", y="DeepDrainage", col='"Deep drainage"'))+
      geom_line(aes_string(x="Date", y="Runoff", col='"Runoff"'))+
      scale_color_manual(name="", values=c("Export"="black", "Deep drainage" = "blue", "Runoff" = "red"))+
      ylab(ylab)+ xlab(xlab)+
      theme_bw()
    return(g)
  } 
  else if(type=="SoilVol") {
    if(is.null(ylab)) ylab = "Soil water content (mm)"
    MLM = data.frame("Total" = Soil$MLTot, 
                     Soil[,paste("ML",1:nlayers,sep=".")])
    if(!is.null(dates)) MLM = MLM[row.names(MLM) %in% as.character(dates),]
    if(!is.null(summary.freq)) MLM = .temporalSummary(MLM, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(MLM), ylab = ylab, ylim = ylim,
                              xlab=xlab, labels = c("Total", paste("Layer", 1:nlayers))))
  } 
}
.plot_soil<-function(Soil, input_soil, input_control, type,  
                     dates = NULL, 
                     xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                     summary.freq = NULL, ...) {
  Soil = as.data.frame(Soil)
  nlayers = length(input_soil$W)
  if(type=="SoilPsi") {
    PsiM = Soil[,paste("psi",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Soil water potential (MPa)"    
    if(!is.null(dates)) PsiM = PsiM[row.names(PsiM) %in% as.character(dates),]
    if(!is.null(summary.freq)) PsiM = .temporalSummary(PsiM, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(PsiM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="SoilTheta") {
    WM = Soil[,paste("W",1:nlayers,sep=".")]
    if(!is.null(dates)) WM = WM[row.names(WM) %in% as.character(dates),]
    theta_FC = soil_thetaFC(input_soil, model = input_control$soilFunctions)
    WM = 100*sweep(WM, 2,theta_FC, "*")
    if(!is.null(summary.freq)) WM = .temporalSummary(WM, summary.freq, mean, na.rm=TRUE)
    if(is.null(ylab)) ylab = "Soil moisture (% volume)"
    return(.multiple_dynamics(as.matrix(WM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="SoilRWC") {
    WM = Soil[,paste("W",1:nlayers,sep=".")]
    if(!is.null(dates)) WM = WM[row.names(WM) %in% as.character(dates),]
    if(!is.null(summary.freq)) WM = .temporalSummary(WM, summary.freq, mean, na.rm=TRUE)
    if(is.null(ylab)) ylab = "Soil moisture (% field capacity)"
    return(.multiple_dynamics(as.matrix(WM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="PlantExtraction") {
    extrBal = Soil[,paste("PlantExt",1:nlayers,sep=".")]
    if(!is.null(dates)) extrBal = extrBal[row.names(extrBal) %in% as.character(dates),]
    if(is.null(ylab)) ylab = "Extraction from soil layer (mm)"    
    if(!is.null(summary.freq)) extrBal = .temporalSummary(extrBal, summary.freq, sum, na.rm=TRUE)
    g<-.multiple_dynamics(as.matrix(extrBal),  xlab = xlab, ylab = ylab, ylim = ylim,
                          labels = paste("Layer", 1:nlayers))
    g<-g+geom_abline(slope=0, intercept=0, col="gray")
    return(g)
  } 
  else if(type=="HydraulicRedistribution") {
    hydrIn = Soil[,paste("HydraulicInput",1:nlayers,sep=".")]
    if(!is.null(dates)) hydrIn = hydrIn[row.names(hydrIn) %in% as.character(dates),]
    if(is.null(ylab)) ylab = "Hydraulic input (mm)"    
    if(!is.null(summary.freq)) hydrIn = .temporalSummary(hydrIn, summary.freq, sum, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(hydrIn),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
}
.plot_stand<-function(Stand, type,  
                      dates = NULL, 
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                      summary.freq = NULL, ...) {
  Stand = as.data.frame(Stand)
  if(type=="LAI") {
    if(is.null(ylab)) ylab = expression(paste("Leaf Area Index   ",(m^{2}%.%m^{-2})))
    df = Stand[,c("LAI", "LAIexpanded", "LAIdead")]
    names(df)<-c("Total (live+dead)", "Live unfolded","Dead standing")
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df), ylab = ylab, ylim = ylim))
  } 
}
.plot_plant_om<-function(OM, PlantsLAI, spnames,
                         type,  bySpecies = FALSE,
                         dates = NULL, 
                         xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                         summary.freq = NULL, ...) {
  OM = as.data.frame(OM)
  if(bySpecies) OM = .averageByLAISpecies(OM, PlantsLAI, spnames)
  if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
  if(!is.null(summary.freq)) OM = .temporalSummary(OM, summary.freq, mean, na.rm=TRUE)
  if(is.null(ylab)) ylab = .getYLab(type)
  return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
}
.plot_plant_om_sum<-function(OM, spnames,
                             type,  bySpecies = FALSE,
                             dates = NULL, 
                             xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                             summary.freq = NULL, ...) {
  OM = as.data.frame(OM)
  if(bySpecies) OM = .sumBySpecies(OM, spnames)
  if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
  if(!is.null(summary.freq)) OM = .temporalSummary(OM, summary.freq, mean, na.rm=TRUE)
  if(is.null(ylab)) ylab = .getYLab(type)
  return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
}
.plot_temperature<-function(Temperature, type,  
                      dates = NULL, 
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                      summary.freq = NULL, ...) {
  Temperature = as.data.frame(Temperature)
  if(type=="Temperature") {
    if(is.null(ylab)) ylab = "Temperature (Celsius)"
    df = data.frame(row.names=row.names(Temperature))
    df[["Above-canopy"]] = Temperature$Tatm_mean
    df[["Inside-canopy"]] = Temperature$Tcan_mean
    df[["Soil"]] = Temperature$Tsoil_mean
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  }
  else if(type=="TemperatureRange") {
    if(is.null(ylab)) ylab = "Temperature (Celsius)"
    df1 = data.frame(row.names=row.names(Temperature))
    df1[["Above-canopy"]] = Temperature$Tatm_min
    df1[["Inside-canopy"]] = Temperature$Tcan_min
    df1[["Soil"]] = Temperature$Tsoil_min
    df2 = data.frame(row.names=row.names(Temperature))
    df2[["Above-canopy"]] = Temperature$Tatm_max
    df2[["Inside-canopy"]] = Temperature$Tcan_max
    df2[["Soil"]] = Temperature$Tsoil_max
    if(!is.null(dates)) {
      df1 = df1[row.names(df1) %in% as.character(dates),]
      df2 = df2[row.names(df2) %in% as.character(dates),]
    }
    if(!is.null(summary.freq)) {
      df1 = .temporalSummary(df1, summary.freq, mean, na.rm=TRUE)
      df2 = .temporalSummary(df2, summary.freq, mean, na.rm=TRUE)
    }
    return(.multiple_dynamics_range(as.matrix(df1),  as.matrix(df2), xlab = xlab, ylab=ylab, ylim = ylim))
  }
  else if(type=="AirTemperature") {
    if(is.null(ylab)) ylab = "Above-canopy temperature (Celsius)"
    df = data.frame(row.names=row.names(Temperature))
    df[["Mean"]] = Temperature$Tatm_mean
    df[["Minimum"]] = Temperature$Tatm_min
    df[["Maximum"]] = Temperature$Tatm_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="CanopyTemperature") {
    if(is.null(ylab)) ylab = "Canopy temperature (Celsius)"
    df = data.frame(row.names=row.names(Temperature))
    df[["Mean"]] = Temperature$Tcan_mean
    df[["Minimum"]] = Temperature$Tcan_min
    df[["Maximum"]] = Temperature$Tcan_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="SoilTemperature") {
    if(is.null(ylab)) ylab = "Soil temperature (Celsius)"
    df = data.frame(row.names=row.names(Temperature))
    df[["Mean"]] = Temperature$Tsoil_mean
    df[["Minimum"]] = Temperature$Tsoil_min
    df[["Maximum"]] = Temperature$Tsoil_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
}
.plot_energybalance<-function(EnergyBalance, type,  
                            dates = NULL, 
                            xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                            summary.freq = NULL, ...) {
  EnergyBalance = as.data.frame(EnergyBalance)
  if(type=="CanopyEnergyBalance") {
    if(is.null(ylab)) ylab = expression(MJ%.%m^{-2})    
    df = data.frame(row.names=row.names(EnergyBalance))
    df[["Balance"]] = EnergyBalance$Ebalcan
    df[["SWR abs."]] = EnergyBalance$SWRcan 
    df[["Net LWR"]] = EnergyBalance$LWRcan
    df[["Latent heat"]] = -EnergyBalance$LEcan
    df[["Convection can./atm."]] = -EnergyBalance$Hcan
    df[["Convection soil/can."]] = -EnergyBalance$Hcansoil
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    g<-.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim)
    return(g)
  } 
  else if(type=="SoilEnergyBalance") {
    if(is.null(ylab)) ylab = expression(MJ%.%m^{-2})    
    df = data.frame(row.names=row.names(EnergyBalance))
    df[["Balance"]] = EnergyBalance$Ebalsoil
    df[["SWR abs."]] = EnergyBalance$SWRsoil
    df[["Net LWR"]] = EnergyBalance$LWRsoil
    df[["Convection soil/can."]] = EnergyBalance$Hcansoil
    df[["Latent heat"]] = -EnergyBalance$LEsoil
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    g<-.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim)
    return(g)
  }
}

.sumSubdailyBySpecies<-function(OM, spnames) {
  if(ncol(OM)>2) OM = cbind(OM[,1],t(apply(OM[,-1],1, tapply, spnames, sum, na.rm=T)))
  else colnames(OM)[2] = spnames[1]
  return(OM)
}
.averageSubdailyBySpecies<-function(OM, spnames) {
  if(ncol(OM)>2) OM = cbind(OM[,1],t(apply(OM[,-1],1, tapply, spnames, mean, na.rm=T)))
  else colnames(OM)[2] = spnames[1]
  return(OM)
}
.plotsubdaily<-function(x, type="PlantTranspiration", cohorts = NULL, bySpecies = FALSE,
                        dates = NULL, xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL) {
  if(("spwb" %in% class(x)) || ("pwb" %in% class(x))) {
    input = x$spwbInput
    type = match.arg(type, .getSubdailySPWBPlotTypes())  
  } else {
    input = x$growthInput
    type = match.arg(type, .getSubdailyGROWTHPlotTypes())  
  }
  
  if(is.null(cohorts)) {cohorts = row.names(input$cohorts)}
  else cohorts = {row.names(input$cohorts)[cohorts]}
  spnames = as.character(input$cohorts[cohorts,"Name"])
  
  if(type=="PlantTranspiration") {
    m = extractSubdaily(x, "E", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = expression(paste("Plant transpiration ",(L%.%m^{-2})))
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantGrossPhotosynthesis") {
    m = extractSubdaily(x, "Ag", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = expression(paste("Gross photosynthesis ",(gC%.%m^{-2})))
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantNetPhotosynthesis") {
    m = extractSubdaily(x, "An", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = expression(paste("Net photosynthesis ",(gC%.%m^{-2})))
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsiAverage") {
    m = extractSubdaily(x, "LeafPsi", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = expression(paste("Leaf water potential ",(MPa)))
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("StemPsi", "RootPsi", "LeafRWC", "StemRWC","LeafSympRWC", "StemSympRWC")) {
    m = extractSubdaily(x, type, dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantWaterBalance") {
    m = extractSubdaily(x, "PWB", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = expression(paste("Plant water balance ",(L%.%m^{-2})))
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="SoilPlantConductance") {
    m = extractSubdaily(x, "dEdP", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = expression(paste("Soil-plant conductance ",(mmol%.%m^{-2}%.%s^{-1})))
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="Temperature") {
    m = extractSubdaily(x, "Temperature", dates)
    m = m[,c("datetime","Tatm", "Tcan", "Tsoil.1")]
    names(m) = c("datetime", "Above-canopy","Inside-canopy", "Soil")
    if(is.null(ylab)) ylab = "Temperature (Celsius)"
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantExtraction") {
    m = extractSubdaily(x, "ExtractionInst", dates)
    if(is.null(ylab)) ylab = expression(paste("Extraction from soil layers ",(L%.%m^{-2})))
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsi") {
    mSu = extractSubdaily(x, "SunlitLeaves$Psi", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$Psi", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Leaf water potential ", (MPa)))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafAbsorbedSWR") {
    mSu = extractSubdaily(x, "SunlitLeaves$Abs_SWR", dates)
    mSh = extractSubdaily(x, "ShadeLeaves$Abs_SWR", dates)
    if(is.null(ylab)) ylab=expression(paste("Absorbed SWR per leaf area ",(W%.%m^{-2})))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafNetLWR") {
    mSu = extractSubdaily(x, "SunlitLeaves$Net_LWR", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$Net_LWR", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Net LWR per leaf area ",(W%.%m^{-2})))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafTranspiration") {
    mSu = extractSubdaily(x, "SunlitLeaves$E", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$E", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Leaf transpiration ",(mmol%.%m^{-2}%.%s^{-1})))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafGrossPhotosynthesis") {
    mSu = extractSubdaily(x, "SunlitLeaves$Ag", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$Ag", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Leaf gross photosynthesis ",(mu%.%mol%.%m^{-2}%.%s^{-1})))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafNetPhotosynthesis") {
    mSu = extractSubdaily(x, "SunlitLeaves$An", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$An", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Leaf net photosynthesis ",(mu%.%mol%.%m^{-2}%.%s^{-1})))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafStomatalConductance") {
    mSu = extractSubdaily(x, "SunlitLeaves$Gsw", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$Gsw", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Stomatal conductance ", (mol%.%m^{-2}%.%s^{-1})))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafTemperature") {
    mSu = extractSubdaily(x, "SunlitLeaves$Temp", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$Temp", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab="Leaf temperature (degrees C)"
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafVPD") {
    mSu = extractSubdaily(x, "SunlitLeaves$VPD", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$VPD", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab="Leaf vapour pressure deficit (kPa)"
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafCi") {
    mSu = extractSubdaily(x, "SunlitLeaves$Ci", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$Ci", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Intercellular CO2 concentration  ", (ppm)))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafIntrinsicWUE") {
    mSu = extractSubdaily(x, "SunlitLeaves$iWUE", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$iWUE", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("iWUE  ", (mu%.%mol%.%mol^{-1})))
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("GrossPhotosynthesis", "MaintenanceRespiration", "GrowthCosts", "CarbonBalance",
                      "SugarLeaf","StarchLeaf","SugarSapwood","StarchSapwood", "SugarTransport")) {
    m = extractSubdaily(x, type, dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
}