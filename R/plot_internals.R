.getWaterBalancePlotTypes<-function() {
  TYPES = c("PET & Precipitation" = "PET_Precipitation",
            "PET and Net rain" = "PET_NetRain",
            "Snow" = "Snow",
            "Water exported" = "Export",
            "Evapotranspiration" = "Evapotranspiration")
  return(TYPES)
}
.getSoilPlotTypes<-function(model = "pwb", transpirationMode = "Granier") {
  TYPES = c(
    "Soil water potential" = "SoilPsi",
    "Soil relative water content" = "SoilRWC",
    "Soil relative extractable water" = "SoilREW",
    "Soil moisture (m3/m3) content" = "SoilTheta",
    "Soil volume (mm) content" = "SoilVol")
  TYPES = c(TYPES, 
            "Plant extraction from soil"= "PlantExtraction",
            "Hydraulic redistribution" = "HydraulicRedistribution")
  return(TYPES)
}
.getStandPlotTypes<-function(model = "pwb") {
  TYPES = c("Stand LAI"="LAI",
            "Ground-level irradiance" = "GroundIrradiance")
  if(model %in% c("growth", "fordyn")) {
    TYPES = c(TYPES, 
              "Carbon balance" = "CarbonBalance",
              "Biomass balance" = "BiomassBalance")
  }
  return(TYPES)
}
.getPlantPlotTypes<-function(transpirationMode = "Granier") {
  TYPES = c("Plant LAI" = "PlantLAI",
            "Plant LAI (live)" = "PlantLAIlive",
            "Transpiration" = "PlantTranspiration",
            "Transpiration per leaf" = "TranspirationPerLeaf",
            "Plant water balance" = "PlantWaterBalance",
            "Gross photosynthesis" = "PlantGrossPhotosynthesis",
            "Gross photosynthesis per leaf" = "GrossPhotosynthesisPerLeaf")
  if(transpirationMode %in% c("Sperry","Sureau")) {
    TYPES <-c(TYPES,
              "Net photosynthesis" = "PlantNetPhotosynthesis",
              "Net photosynthesis per leaf" = "NetPhotosynthesisPerLeaf",
              "Water use efficiency" = "PlantWUE",
              "Absorbed short-wave radiation" = "PlantAbsorbedSWR",
              "Absorbed short-wave radiation per leaf" = "AbsorbedSWRPerLeaf",
              "Net long-wave radiation" = "PlantNetLWR",
              "Net long-wave radiation per leaf" = "NetLWRPerLeaf")
  } else {
    TYPES <-c(TYPES,
              "Fraction of PAR" = "FPAR",
              "Absorbed SWR fraction" = "AbsorbedSWRFraction")
  }           
  
  if(transpirationMode %in% c("Sperry","Sureau")) {
    TYPES <-c(TYPES,
              "Minimum leaf water potential" = "LeafPsiMin",
              "Maximum leaf water potential" = "LeafPsiMax",
              "Leaf water potential range" = "LeafPsiRange",
              "Midday upper stem water potential" = "StemPsi",
              "Midday root crown water potential" = "RootPsi",
              "Stem relative water content" = "StemRWC",
              "Leaf relative water content" = "LeafRWC",
              "Live fuel moisture content" = "LFMC")
  } else {
    TYPES <-c(TYPES,
              "Plant water potential" = "PlantPsi",
              "Stem relative water content" = "StemRWC",
              "Leaf relative water content" = "LeafRWC",
              "Live fuel moisture content" = "LFMC")
  }
  if(transpirationMode %in% c("Sperry","Sureau")) {
    TYPES <-c(TYPES,
              "Soil-plant conductance" = "SoilPlantConductance")
  }
  TYPES <-c(TYPES,
            "Plant stress" = "PlantStress",
            "Leaf percent conductance loss" = "LeafPLC",
            "Stem percent conductance loss" = "StemPLC")
  return(TYPES)
}
.getSunlitShadePlotTypes<-function(transpirationMode = "Granier"){
  TYPES = character(0)
  if(transpirationMode %in% c("Sperry","Sureau")) {
    TYPES = c(TYPES,
              "Minimum leaf temperature (sunlit)" ="TempMin_SL", 
              "Minimum leaf temperature (shade)" = "TempMin_SH", 
              "Maximum leaf temperature (sunlit)" ="TempMax_SL",
              "Maximum leaf temperature (shade)" ="TempMax_SH",
              "Minimum stomatal conductance (sunlit)" ="GSWMin_SL", 
              "Minimum stomatal conductance (shade)" ="GSWMin_SH", 
              "Maximum stomatal conductance (sunlit)" ="GSWMax_SL", 
              "Maximum stomatal conductance (shade)" ="GSWMax_SH",
              "Minimum leaf water potential (sunlit)" ="LeafPsiMin_SL", 
              "Minimum leaf water potential (shade)" ="LeafPsiMin_SH", 
              "Maximum leaf water potential (sunlit)" ="LeafPsiMax_SL", 
              "Maximum leaf water potential (shade)" ="LeafPsiMax_SH")
  }
  return(TYPES)
}
.getEnergyPlotTypes<-function(transpirationMode = "Granier") {
  TYPES = character(0)
  if(transpirationMode %in% c("Sperry","Sureau")) {
    TYPES = c(TYPES,
              "Temperature" = "Temperature",
              "Temperature range" = "TemperatureRange",
              "Above-canopy air temperature" = "AirTemperature",
              "Within-canopy air temperature" = "CanopyTemperature",
              "Soil surface temperature" = "SoilTemperature",
              "Canopy energy balance components" = "CanopyEnergyBalance",
              "Soil energy balance components" = "SoilEnergyBalance")
  }
  return(TYPES)
}
.getStructuralGROWTHPlotTypes<-function(transpirationMode = "Granier"){
  TYPES = c("Leaf biomass per individual" = "LeafBiomass",
            "Sapwood biomass per individual" = "SapwoodBiomass",
            "Fine root biomass per individual" = "FineRootBiomass",
            "Sapwood area" ="SapwoodArea",
            "Leaf area" = "LeafArea", 
            "Fine root area" ="FineRootArea")
  TYPES = c(TYPES,             
            "Diameter at breast height" = "DBH",
            "Plant height" =  "Height", 
            "Sapwood area / Leaf area" ="HuberValue",
            "Fine root area / Leaf area" ="RootAreaLeafArea")
  return(TYPES)
}
.getGrowthMortalityGROWTHPlotTypes<-function(transpirationMode = "Granier"){
  TYPES = c("Sapwood area growth rate" ="SAgrowth",
            "Leaf area growth rate" = "LAgrowth",
            "Leaf area growth rate" = "LAgrowth", 
            "Fine root area growth rate" ="FRAgrowth")
  TYPES = c(TYPES,
            "Starvation mortality rate" = "StarvationRate", 
            "Dessication mortality rate" = "DessicationRate", 
            "Mortality rate" = "MortalityRate")
  return(TYPES)
}
.getCohortBiomassBalanceGROWTHPlotTypes<-function(transpirationMode = "Granier"){
  return(c("Structural biomass balance (g/m2)" = "StructuralBiomassBalance",
           "Labile biomass balance (g/m2)" = "LabileBiomassBalance",
           "Plant biomass balance (g/m2)" = "PlantBiomassBalance",
           "Mortality biomass loss (g/m2)" = "MortalityBiomassLoss",
           "Cohort biomass balance (g/m2)" = "CohortBiomassBalance"))
}
.getLabileGROWTHPlotTypes<-function(transpirationMode = "Granier") {
  TYPES= c("Gross photosynthesis per dry" = "GrossPhotosynthesis",
           "Maintenance respiration per dry" = "MaintenanceRespiration",
           "Photosynthesis-maintenance ratio" = "PhotosynthesisMaintenanceRatio",
           "Growth costs per dry" = "GrowthCosts",
           "Labile carbon balance per dry" = "LabileCarbonBalance",
           "Leaf sugar concentration" = "SugarLeaf",
           "Leaf starch concentration" = "StarchLeaf",
           "Sapwood sugar concentration" = "SugarSapwood",
           "Sapwood starch concentration" = "StarchSapwood",
           "Sugar transport" = "SugarTransport",
           "Root exudation" = "RootExudation")
  return(TYPES)
}

.getDailyPWBPlotTypes<-function(transpirationMode = "Granier") {
  TYPES = c(.getSoilPlotTypes("pwb", transpirationMode),
            .getStandPlotTypes("pwb"),
            .getPlantPlotTypes(transpirationMode),
            .getSunlitShadePlotTypes(transpirationMode),
            .getEnergyPlotTypes(transpirationMode))
  return(TYPES)
}

.getDailySPWBPlotTypes<-function(transpirationMode = "Granier") {
  TYPES = c(.getWaterBalancePlotTypes(),
            .getSoilPlotTypes("spwb", transpirationMode),
            .getStandPlotTypes("spwb"),
            .getPlantPlotTypes(transpirationMode),
            .getSunlitShadePlotTypes(transpirationMode),
            .getEnergyPlotTypes(transpirationMode))
  return(TYPES)
}

.getUniqueDailyGROWTHPlotTypes<-function(transpirationMode = "Granier"){
  TYPES = c("Carbon balance" = "CarbonBalance",
            "Biomass balance" = "BiomassBalance",
            .getLabileGROWTHPlotTypes(transpirationMode),
            .getCohortBiomassBalanceGROWTHPlotTypes(transpirationMode),
            .getStructuralGROWTHPlotTypes(transpirationMode),
            .getGrowthMortalityGROWTHPlotTypes(transpirationMode))
  return(TYPES)
}
.getDailyGROWTHPlotTypes<-function(transpirationMode = "Granier"){
  TYPES = c(.getWaterBalancePlotTypes(),
            .getSoilPlotTypes("growth", transpirationMode),
            .getStandPlotTypes("growth"),
            .getPlantPlotTypes(transpirationMode),
            .getSunlitShadePlotTypes(transpirationMode),
            .getEnergyPlotTypes(transpirationMode),
            .getUniqueDailyGROWTHPlotTypes(transpirationMode))
  return(TYPES)
}

.getUniqueFORDYNPlotTypes<-function(transpirationMode = "Granier"){
  return(c("Stand basal area" = "StandBasalArea", 
           "Dominant tree height" = "DominantTreeHeight",
           "Dominant tree diameter" = "DominantTreeDiameter",
           "Quadratic mean tree diameter" = "QuadraticMeanTreeDiameter",
           "Hart-Becking index" = "HartBeckingIndex",
           "Basal area of trees by species" = "SpeciesBasalArea", 
           "Basal area of trees by cohort" = "CohortBasalArea", 
           "Stand density of trees" = "StandDensity",
           "Density of trees by species" = "SpeciesDensity",
           "Density of trees by cohort" = "CohortDensity",
           "Stand shrub cover" = "StandShrubCover",
           "Shrub cover by species" = "SpeciesShrubCover",
           "Shrub cover by cohort" = "CohortShrubCover",
           "Number of tree species" = "NumTreeSpecies",
           "Number of shrub species" = "NumShrubSpecies",
           "Number of tree cohorts" = "NumTreeCohorts",
           "Number of shrub cohorts" = "NumShrubCohorts",
           "Number of cohorts by species" = "NumCohortsSpecies"))
}

.getSubdailySoilPlotTypes<-function(){
  return(c("Plant extraction from soil" = "PlantExtraction"))
}
.getSubdailyPlantPlotTypes<-function(){
  TYPES = c("Average leaf water potential" = "LeafPsiAverage",
            "Root crown water potential" = "RootPsi", 
            "Stem water potential" = "StemPsi", 
            "Leaf symplastic water potential" = "LeafSympPsi", 
            "Stem symplastic water potential" = "StemSympPsi",
            "Leaf percent conductance loss" = "LeafPLC",
            "Stem percent conductance loss" = "StemPLC",
            "Stem relative water content" = "StemRWC", 
            "Leaf relative water content" = "LeafRWC",
            "Stem symplastic relative water content" = "StemSympRWC", 
            "Leaf symplastic relative water content" = "LeafSympRWC",
            "Soil-plant conductance" = "SoilPlantConductance",
            "Plant transpiration" = "PlantTranspiration", 
            "Transpiration per leaf area" = "TranspirationPerLeaf",
            "Plant gross photosynthesis" = "PlantGrossPhotosynthesis",
            "Gross photosynthesis per leaf area" = "GrossPhotosynthesisPerLeaf",
            "Plant net photosynthesis" = "PlantNetPhotosynthesis",
            "Net photosynthesis per leaf area" = "NetPhotosynthesisPerLeaf", 
            "Plant water balance" = "PlantWaterBalance", 
            "Water balance per leaf area"= "WaterBalancePerLeaf")
  return(TYPES)
}
.getSubdailySunlitShadePlotTypes<-function(){
  TYPES = c("Leaf area" = "LeafLAI",
            "Leaf water potential" = "LeafPsi",
            "Leaf transpiration" = "LeafTranspiration",
            "Leaf gross photosynthesis" = "LeafGrossPhotosynthesis", 
            "Leaf net photosynthesis" = "LeafNetPhotosynthesis",
            "Absorbed SWR per leaf area" = "LeafAbsorbedSWR",
            "Absorbed PAR per leaf area" = "LeafAbsorbedPAR",
            "Net LWR per leaf area" = "LeafNetLWR",
            "Leaf internal CO2" = "LeafCi", 
            "Leaf intrinsic WUE" = "LeafIntrinsicWUE",
            "Leaf vapour pressure deficit" = "LeafVPD",
            "Leaf stomatal conductance" = "LeafStomatalConductance", 
            "Leaf temperature" = "LeafTemperature")
  return(TYPES)
}
.getSubdailyEnergyBalancePlotTypes<-function(){
  TYPES = c("Air/canopy/soil temperature" ="Temperature",
            "CanopyEnergyBalance",
            "SoilEnergyBalance")
  return(TYPES)
}
.getSubdailyLabilePlotTypes<-function(){
  TYPES = c("Gross photosynthesis per dry" = "GrossPhotosynthesis",
            "Maintenance respiration per dry" = "MaintenanceRespiration",
            "Growth costs per dry" = "GrowthCosts",
            "Root exudation per dry" = "RootExudation",
            "Labile carbon balance per dry" = "LabileCarbonBalance",
            "Leaf sugar concentration" = "SugarLeaf",
            "Leaf starch concentration" = "StarchLeaf",
            "Sapwood sugar concentration" = "SugarSapwood",
            "Sapwood starch concentration" = "StarchSapwood",
            "Sugar transport" = "SugarTransport")
  return(TYPES)
}

.getSubdailySPWBPlotTypes<-function(){
  TYPES = c(.getSubdailySoilPlotTypes(),
            .getSubdailyPlantPlotTypes(),
            .getSubdailySunlitShadePlotTypes(),
            .getSubdailyEnergyBalancePlotTypes())
  return(TYPES)
}
.getSubdailyGROWTHPlotTypes<-function(){
  TYPES = c(.getSubdailySoilPlotTypes(),
            .getSubdailyPlantPlotTypes(),
            .getSubdailySunlitShadePlotTypes(),
            .getSubdailyEnergyBalancePlotTypes(),
            .getSubdailyLabilePlotTypes())
  return(TYPES)
}


.getYLab<-function(type) {
  ylab="Unknown"
  if(type=="PlantTranspiration") ylab = expression(paste("Plant transpiration ",(L%.%m^{-2}%.%d^{-1})))
  else if(type=="TranspirationPerLeaf") ylab = expression(paste("Transpiration per leaf area ",(L%.%m^{-2}%.%d^{-1})))
  else if(type=="LeafTranspiration") ylab = expression(paste("Leaf transpiration ",(mmol%.%m^{-2}%.%s^{-1})))
  else if(type=="PlantPhotosynthesis") ylab = expression(paste("Plant photosynthesis ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantGrossPhotosynthesis") ylab = expression(paste("Plant gross photosynthesis ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="LeafGrossPhotosynthesis") ylab = expression(paste("Leaf gross photosynthesis ",(mu%.%mol%.%m^{-2}%.%s^{-1})))
  else if(type=="PlantNetPhotosynthesis") ylab = expression(paste("Plant net photosynthesis ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="LeafNetPhotosynthesis") ylab = expression(paste("Leaf net photosynthesis ",(mu%.%mol%.%m^{-2}%.%s^{-1})))
  else if(type=="PhotosynthesisPerLeaf") ylab = expression(paste("Photosynthesis per leaf area ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="GrossPhotosynthesisPerLeaf") ylab = expression(paste("Gross photosynthesis per leaf area ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="NetPhotosynthesisPerLeaf") ylab = expression(paste("Net photosynthesis per leaf area ",(g*C%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantWUE") ylab = expression(paste("Plant water use efficiency ",(g*C%.%L^{-1})))
  else if(type=="FPAR") ylab = "Average fraction of PAR (%)"
  else if(type=="AbsorbedSWRFraction") ylab = "Fraction of absorbed SWR (%)"
  else if(type=="PlantAbsorbedSWR") ylab = expression(paste("Plant absorbed SWR ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantNetLWR") ylab = expression(paste("Plant net LWR ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="PlantExtraction") ylab = expression(paste("Extraction from soil layers   ",(L%.%m^{-2})))
  else if(type=="PlantWaterBalance") ylab = expression(paste("Plant water balance   ",(L%.%m^{-2})))
  else if(type=="WaterBalancePerLeaf") ylab = expression(paste("Water balance per leaf area ",(L%.%m^{-2})))
  else if(type=="PlantLAI") ylab = expression(paste("Leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="PlantLAIlive") ylab = expression(paste("(Live) leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="AbsorbedSWRPerLeaf") ylab = expression(paste("Absorbed SWR per leaf area ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="AbsorbedPARPerLeaf") ylab = expression(paste("Absorbed PAR per leaf area ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="LeafLAI") ylab = expression(paste("Leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="LeafAbsorbedSWR") ylab = expression(paste("Absorbed SWR per leaf area ",(W%.%m^{-2})))
  else if(type=="LeafAbsorbedPAR") ylab = expression(paste("Absorbed PAR per leaf area ",(W%.%m^{-2})))
  else if(type=="NetLWRPerLeaf") ylab = expression(paste("Net LWR per leaf area ",(MJ%.%m^{-2}%.%d^{-1})))
  else if(type=="LeafNetLWR") ylab = expression(paste("Net LWR per leaf area ",(W%.%m^{-2})))
  else if(type=="GrossPhotosynthesis") ylab=expression(paste("Gross photosynthesis ", (gGluc%.%gdry^{-1}%.%d^{-1})))
  else if(type=="MaintenanceRespiration") ylab=expression(paste("Maintenance respiration ", (gGluc%.%gdry^{-1}%.%d^{-1})))
  else if(type=="PhotosynthesisMaintenanceRatio") ylab="Photosynthesis-maintenance ratio "
  else if(type=="GrowthCosts") ylab=expression(paste("Growth costs ", (gGluc%.%gdry^{-1}%.%d^{-1})))
  else if(type=="LabileCarbonBalance") ylab=expression(paste("Labile carbon balance ", (gGluc%.%gdry^{-1}%.%d^{-1})))
  else if(type=="SugarLeaf") ylab=expression(paste("Leaf sugar concentration  ", (mol%.%L^{-1})))
  else if(type=="StarchLeaf") ylab=expression(paste("Leaf starch concentration  ", (mol%.%L^{-1})))
  else if(type=="SugarSapwood") ylab=expression(paste("Sapwood sugar concentration  ", (mol%.%L^{-1})))
  else if(type=="StarchSapwood") ylab=expression(paste("Sapwood starch concentration  ", (mol%.%L^{-1})))
  else if(type=="SugarTransport") ylab=expression(paste("Floem sugar transport rate ", (mmol%.%s^{-1})))
  else if(type=="RootExudation") ylab=expression(paste("Root exudation ", (gGluc%.%gdry^{-1})))
  else if(type=="StructuralBiomassBalance") ylab=expression(paste("Structural biomass change ", (g%.%m^{-2})))
  else if(type=="LabileBiomassBalance") ylab=expression(paste("Labile biomass change ", (g%.%m^{-2})))
  else if(type=="PlantBiomassBalance") ylab=expression(paste("Plant biomass balance ", (g%.%m^{-2})))
  else if(type=="MortalityBiomassLoss") ylab=expression(paste("Mortality biomass loss ", (g%.%m^{-2})))
  else if(type=="CohortBiomassBalance") ylab=expression(paste("Cohort biomass balance ", (g%.%m^{-2})))
  else if(type=="SapwoodArea")  ylab = expression(paste("Sapwood area  ",(cm^2)))
  else if(type=="LeafArea")  ylab = expression(paste("Leaf area  ",(m^2)))
  else if(type=="FineRootArea")  ylab = expression(paste("Fine root area  ",(m^2)))
  else if(type=="LeafBiomass")  ylab = expression(paste("Leaf structural biomass  ", (gdry%.%ind^{-1})))
  else if(type=="SapwoodBiomass")  ylab = expression(paste("Sapwood structural biomass  ", (gdry%.%ind^{-1})))
  else if(type=="FineRootBiomass")  ylab = expression(paste("Fine root biomass  ", (gdry%.%ind^{-1})))
  else if(type=="DBH")  ylab = expression(paste("Diameter at breast height  ", (cm)))
  else if(type=="Height")  ylab = expression(paste("Plant height  ", (cm)))
  else if(type=="HuberValue")  ylab = expression(paste("Huber value  ",(cm^2 %.% m^{-2})))
  else if(type=="RootAreaLeafArea")  ylab = expression(paste("Root area / Leaf area  ",(m^2 %.% m^{-2})))
  else if(type=="SAgrowth") ylab = expression(paste("Sapwood area growth rate ",(cm^2  %.% d^{-1})))
  else if(type=="LAgrowth") ylab = expression(paste("Leaf area growth rate ",(m^2 %.% d^{-1})))
  else if(type=="FRAgrowth") ylab = expression(paste("Fine root area growth rate ",(m^2 %.% d^{-1})))
  else if(type=="DessicationRate") ylab = expression(paste("Dessication mortality rate ",(ind %.% d^{-1})))
  else if(type=="StarvationRate") ylab = expression(paste("Starvation mortality rate ",(ind %.% d^{-1})))
  else if(type=="MortalityRate") ylab = expression(paste("Mortality rate ",(ind %.% d^{-1})))
  else if(type=="LeafPI0")  ylab = expression(paste("Leaf osmotic potential at full turgor  ",(MPa)))
  else if(type=="StemPI0")  ylab = expression(paste("Stem osmotic potential at full turgor  ",(MPa)))
  else if(type=="LeafPLC") ylab = "Percent loss conductance in leaves [%]"
  else if(type=="StemPLC") ylab = "Percent loss conductance in stem [%]"
  else if(type=="StemRWC") ylab = "Relative water content in stem [%]"
  else if(type=="StemSympRWC") ylab = "Relative water content in stem symplasm [%]"
  else if(type=="StemSympPsi") ylab = "Stem symplastic water potential (MPa)"
  else if(type=="LeafRWC") ylab = "Relative water content in leaf [%]"
  else if(type=="LFMC") ylab = "Live fuel moisture content [%]"
  else if(type=="LeafSympRWC") ylab = "Relative water content in leaf symplasm [%]"
  else if(type=="LeafSympPsi") ylab = "Leaf symplastic water potential (MPa)"
  else if(type=="PlantPsi") ylab = "Plant water potential (MPa)"
  else if(type=="PlantStress") ylab = "Drought stress [%]"
  else if(type=="StemPsi") ylab = "Stem water potential (MPa)"
  else if(type=="RootPsi") ylab = "Root crown water potential (MPa)"
  else if(type=="LeafPsiAverage") ylab = "Average leaf water potential (MPa)"
  else if(type=="LeafPsi") ylab = "Leaf water potential (MPa)"
  else if(type=="LeafPsiRange") ylab = "Leaf water potential range (MPa)"
  else if(type=="SoilPlantConductance") ylab = expression(paste("Soil-plant conductance ",(mmol%.%m^{-2}%.%s^{-1})))
  else if(type=="LeafPsiMin") ylab = "Minimum (midday) leaf water potential (MPa)"
  else if(type=="LeafPsiMax") ylab = "Maximum (predawn) leaf water potential (MPa)"
  else if(type=="LeafPsiMin_SL") ylab = "Minimum (midday) sunlit leaf water potential (MPa)"
  else if(type=="LeafPsiMax_SL") ylab = "Maximum (predawn) sunlit leaf water potential (MPa)"
  else if(type=="LeafPsiMin_SH") ylab = "Minimum (midday) shade leaf water potential (MPa)"
  else if(type=="LeafPsiMax_SH") ylab = "Maximum (predawn) shade leaf water potential (MPa)"
  else if(type=="LeafStomatalConductance") ylab = expression(paste("Stomatal conductance ", (mol%.%m^{-2}%.%s^{-1})))
  else if(type=="LeafCi") ylab = expression(paste("Intercellular CO2 concentration  ", (ppm)))
  else if(type=="LeafVPD") ylab = "Leaf vapour pressure deficit (kPa)"
  else if(type=="LeafTemperature") ylab ="Leaf temperature (degrees C)"
  else if(type=="LeafIntrinsicWUE") ylab = expression(paste("iWUE  ", (mu%.%mol%.%mol^{-1})))
  else if(type=="GSWMin_SH") ylab = expression(paste("Minimum shade leaf stomatal conductance ",(mol%.%m^{-2}%.%s^{-1})))
  else if(type=="GSWMax_SH") ylab = expression(paste("Maximum shade leaf stomatal conductance ",(mol%.%m^{-2}%.%s^{-1})))
  else if(type=="GSWMin_SL") ylab = expression(paste("Minimum sunlit leaf stomatal conductance ",(mol%.%m^{-2}%.%s^{-1})))
  else if(type=="GSWMax_SL") ylab = expression(paste("Maximum sunlit leaf stomatal conductance ",(mol%.%m^{-2}%.%s^{-1})))
  else if(type=="TempMin_SH") ylab = expression(paste("Minimum shade leaf temperature (Celsius)"))
  else if(type=="TempMax_SH") ylab = expression(paste("Maximum shade leaf temperature (Celsius)"))
  else if(type=="TempMin_SL") ylab = expression(paste("Minimum sunlit leaf temperature (Celsius)"))
  else if(type=="TempMax_SL") ylab = expression(paste("Maximum sunlit leaf temperature (Celsius)"))
  else if(type=="NumTreeSpecies") ylab = expression(paste("Number of tree species"))
  else if(type=="NumShrubSpecies") ylab = expression(paste("Number of shrub species"))
  else if(type=="NumTreeCohorts") ylab = expression(paste("Number of tree cohorts"))
  else if(type=="NumShrubCohorts") ylab = expression(paste("Number of shrub cohorts"))
  else if(type=="NumCohortsSpecies") ylab = expression(paste("Number of cohorts by species"))
  else if(type=="StandBasalArea") ylab = expression(paste("Stand basal area ",(m^{-2}%.%ha^{-1})))
  else if(type=="DominantTreeHeight") ylab = expression(paste("Dominant tree height ",(cm)))
  else if(type=="DominantTreeDiameter") ylab = expression(paste("Dominant tree diameter ",(cm)))
  else if(type=="QuadraticMeanTreeDiameter") ylab = expression(paste("Quadratic mean tree diameter ",(cm)))
  else if(type=="HartBeckingIndex") ylab = expression(paste("Hart-becking index "))
  else if(type=="StandLAI") ylab = expression(paste("Stand leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="StandDensity") ylab = expression(paste("Stand tree density ",(ind%.%ha^{-1})))
  else if(type=="StandShrubCover") ylab = expression(paste("Stand shrub cover (%)"))
  else if(type=="SpeciesBasalArea") ylab = expression(paste("Species basal area ",(m^{-2}%.%ha^{-1})))
  else if(type=="SpeciesLAI") ylab = expression(paste("Species leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="SpeciesDensity") ylab = expression(paste("Species tree density ",(ind%.%ha^{-1})))
  else if(type=="SpeciesShrubCover") ylab = expression(paste("Shrub cover by species (%)"))
  else if(type=="CohortBasalArea") ylab = expression(paste("Cohort basal area ",(m^{-2}%.%ha^{-1})))
  else if(type=="CohortLAI") ylab = expression(paste("Cohort leaf area index ",(m^{-2}%.%m^{-2})))
  else if(type=="CohortDensity") ylab = expression(paste("Cohort tree density ",(ind%.%ha^{-1})))
  else if(type=="CohortShrubCover") ylab = expression(paste("Shrub cover by cohort (%)"))
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
.averageByLAISpecies<-function(OM, PlantsLAIlive, spnames) {
  if(ncol(OM)>1) {
    lai1 = t(apply(PlantsLAIlive,1, tapply, spnames, sum))
    m1 = t(apply(PlantsLAIlive * OM,1, tapply, spnames, sum))
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


.plot_wb<-function(WaterBalance, Snow, type,  
                   dates = NULL, 
                   xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                   summary.freq = NULL, ...) {
  WaterBalance = as.data.frame(WaterBalance)
  Snow = as.data.frame(Snow)
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
      geom_bar(aes(x=.data$Date, y=.data$Precipitation, fill="Precipitation"), stat = "identity")+
      geom_bar(aes(x=.data$Date, y=.data$Snow, fill="Snow"), stat = "identity")+
      geom_path(aes(x=.data$Date, y=.data$PET, col="PET"))+
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
      geom_area(aes(x=.data$Date, y=.data$NetRain, fill="NetRain"))+
      geom_path(aes(x=.data$Date, y=.data$PET, col="PET"))+
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
    df[["Woody transpiration"]] = WaterBalance$Transpiration
    df[["Herbaceous transpiration"]] = WaterBalance$HerbTranspiration
    df[["Bare soil evaporation"]] = WaterBalance$SoilEvaporation
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(row.names(df)), breaks=summary.freq)
      df = data.frame(row.names = as.Date(as.character(levels(date.factor))),
                      "Total evapotranspiration" = tapply(df[["Total evapotranspiration"]],INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      "Interception evaporation" = tapply(df[["Interception evaporation"]],INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      "Woody transpiration" = tapply(df[["Woody transpiration"]],INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      "Herbaceous transpiration" = tapply(df[["Herbaceous transpiration"]],INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      "Bare soil evaporation" = tapply(df[["Bare soil evaporation"]],INDEX=date.factor, FUN=sum, na.rm=TRUE))
    }
    return(.multiple_dynamics(as.matrix(df), ylab=ylab, xlab=xlab, ylim = ylim))
  } 
  else if(type=="Snow") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})  
    df = data.frame(row.names=row.names(WaterBalance))
    df[["Snow"]] = WaterBalance$Snow
    df[["Snowpack"]] = Snow$SWE
    df[["Date"]] = as.Date(row.names(WaterBalance))
    if(!is.null(dates)) df = df[df$Date %in% dates,]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(df$Date), breaks=summary.freq)
      df = data.frame(Date = as.Date(as.character(levels(date.factor))),
                      Snow = tapply(df$Snow,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      Snowpack = tapply(df$Snowpack,INDEX=date.factor, FUN=mean, na.rm=TRUE))
    }
    g<-ggplot(df)+
      geom_area(aes(x=.data$Date, y=.data$Snow, fill="Snow"))+
      geom_path(aes(x=.data$Date, y=.data$Snowpack, col="Snowpack"))+
      scale_fill_manual(name="", values=c("Snow"="black"))+
      scale_color_manual(name="", values=c("Snowpack"="gray"))+
      ylab(ylab)+xlab(xlab)+
      theme_bw()
    return(g)
  } 
  else if(type=="Export") {
    if(is.null(ylab)) ylab =  expression(L%.%m^{-2})    
    df = data.frame(row.names=row.names(WaterBalance))
    df[["InfiltrationExcess"]] = WaterBalance$InfiltrationExcess
    df[["CapillarityRise"]] = -WaterBalance$CapillarityRise
    df[["SaturationExcess"]] = WaterBalance$SaturationExcess
    df[["Export"]] = WaterBalance$DeepDrainage + WaterBalance$Runoff
    df[["DeepDrainage"]] = WaterBalance$DeepDrainage
    df[["Runoff"]] = WaterBalance$Runoff 
    df[["Date"]]= as.Date(row.names(WaterBalance))
    if(!is.null(dates)) df = df[df$Date %in% dates,]
    if(!is.null(summary.freq)) {
      date.factor = cut(as.Date(df$Date), breaks=summary.freq)
      df = data.frame(Date = as.Date(as.character(levels(date.factor))),
                      InfiltrationExcess = tapply(df$InfiltrationExcess,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      SaturationExcess = tapply(df$SaturationExcess,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      CapillarityRise = tapply(df$CapillarityRise,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      Export = tapply(df$Export,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      DeepDrainage = tapply(df$DeepDrainage,INDEX=date.factor, FUN=sum, na.rm=TRUE),
                      Runoff = tapply(df$Runoff,INDEX=date.factor, FUN=sum, na.rm=TRUE))
    }
    g<-ggplot(df)+
      geom_line(aes(x=.data$Date, y=.data$Export, col="Export"))+
      geom_line(aes(x=.data$Date, y=.data$DeepDrainage, col="Deep drainage"))+
      geom_line(aes(x=.data$Date, y=.data$InfiltrationExcess, col="Infiltration excess"))+
      geom_line(aes(x=.data$Date, y=.data$SaturationExcess, col="Saturation excess"))+
      geom_line(aes(x=.data$Date, y=.data$CapillarityRise, col="Capillarity rise"))+
      geom_line(aes(x=.data$Date, y=.data$Runoff, col="Runoff"))+
      scale_color_manual(name="", values=c("Export"="black", "Deep drainage" = "blue", 
                                           "Infiltration excess" = "yellow", "Saturation excess" = "lightblue",
                                           "Runoff" = "red", 
                                           "Capillarity rise" = "darkgreen"))+
      ylab(ylab)+ xlab(xlab)+
      theme_bw()
    return(g)
  } 
}
.plot_soil<-function(Soil, input_soil, input_control, type,  
                     dates = NULL, 
                     xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                     summary.freq = NULL, ...) {
  nlayers = length(input_soil$W)
  if(type=="SoilPsi") {
    PsiM = as.data.frame(Soil$Psi)
    if(is.null(ylab)) ylab = "Soil water potential (MPa)"    
    if(!is.null(dates)) PsiM = PsiM[row.names(PsiM) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) PsiM = .temporalSummary(PsiM, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(PsiM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = c(paste("Layer", 1:nlayers), "Overall")))
  } 
  else if(type=="SoilTheta") {
    SWCM = as.data.frame(Soil$SWC)
    if(!is.null(dates)) SWCM = SWCM[row.names(SWCM) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) SWCM = .temporalSummary(SWCM, summary.freq, mean, na.rm=TRUE)
    if(is.null(ylab)) ylab = "Soil water content (% volume)"
    return(.multiple_dynamics(as.matrix(SWCM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = c(paste("Layer", 1:nlayers), "Overall")))
  } 
  else if(type=="SoilREW") {
    REWM = as.data.frame(100*Soil$REW)
    if(!is.null(dates)) REWM = REWM[row.names(REWM) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) REWM = .temporalSummary(REWM, summary.freq, mean, na.rm=TRUE)
    if(is.null(ylab)) ylab = "Relative extractable water (%)"
    return(.multiple_dynamics(as.matrix(REWM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = c(paste("Layer", 1:nlayers), "Overall")))
  } 
  else if(type=="SoilRWC") {
    RWCM = as.data.frame(100*Soil$RWC)
    if(!is.null(dates)) RWCM = RWCM[row.names(RWCM) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) RWCM = .temporalSummary(RWCM, summary.freq, mean, na.rm=TRUE)
    if(is.null(ylab)) ylab = "Relative water content (% field capacity)"
    return(.multiple_dynamics(as.matrix(RWCM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = c(paste("Layer", 1:nlayers), "Overall")))
  } 
  else if(type=="SoilVol") {
    if(is.null(ylab)) ylab = "Soil water content (mm)"
    MLM = as.data.frame(Soil$ML)
    if(!is.null(dates)) MLM = MLM[row.names(MLM) %in% as.character(dates),]
    if(!is.null(summary.freq)) MLM = .temporalSummary(MLM, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(MLM), ylab = ylab, ylim = ylim,
                              xlab=xlab, labels = c(paste("Layer", 1:nlayers), "Overall")))
  } 
  else if(type=="PlantExtraction") {
    extrBal = as.data.frame(Soil$PlantExt)
    if(!is.null(dates)) extrBal = extrBal[row.names(extrBal) %in% as.character(dates),,drop = FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    if(!is.null(summary.freq)) extrBal = .temporalSummary(extrBal, summary.freq, sum, na.rm=TRUE)
    g<-.multiple_dynamics(as.matrix(extrBal),  xlab = xlab, ylab = ylab, ylim = ylim,
                          labels = c(paste("Layer", 1:nlayers), "Overall"))
    g<-g+geom_abline(slope=0, intercept=0, col="gray")
    return(g)
  } 
  else if(type=="HydraulicRedistribution") {
    hydrIn = as.data.frame(Soil$HydraulicInput)
    if(!is.null(dates)) hydrIn = hydrIn[row.names(hydrIn) %in% as.character(dates),,drop = FALSE]
    if(is.null(ylab)) ylab = "Hydraulic input (mm)"    
    if(!is.null(summary.freq)) hydrIn = .temporalSummary(hydrIn, summary.freq, sum, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(hydrIn),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = c(paste("Layer", 1:nlayers), "Overall")))
  } 
}
.plot_stand<-function(Stand, type,  
                      dates = NULL, 
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                      summary.freq = NULL, ...) {
  Stand = as.data.frame(Stand)
  if(type=="LAI") {
    if(is.null(ylab)) ylab = expression(paste("Leaf Area Index   ",(m^{2}%.%m^{-2})))
    df = Stand[,c("LAI", "LAIherb", "LAIexpanded", "LAIdead")]
    names(df)<-c("Total (herb+unfolded+dead)", "Herbaceous", "Woody plants unfolded","Woody plants dead")
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df), ylab = ylab, ylim = ylim))
  } 
  else if(type=="GroundIrradiance") {
    if(is.null(ylab)) ylab = expression(paste("Ground-level irradiance (%) "))
    df = Stand[,c("LgroundPAR", "LgroundSWR")]
    names(df)<-c("Photosynthetically-active radiation","Short-wave radiation")
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df), ylab = ylab, ylim = ylim))
  } 
}
.plot_plant_om<-function(OM, PlantsLAIlive, spnames,
                         type,  bySpecies = FALSE,
                         dates = NULL, 
                         xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                         summary.freq = NULL, ...) {
  OM = as.data.frame(OM)
  if(bySpecies) OM = .averageByLAISpecies(OM, PlantsLAIlive, spnames)
  if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),,drop =FALSE]
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
  if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),, drop = FALSE]
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
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
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
      df1 = df1[row.names(df1) %in% as.character(dates),,drop = FALSE]
      df2 = df2[row.names(df2) %in% as.character(dates),,drop = FALSE]
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
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="CanopyTemperature") {
    if(is.null(ylab)) ylab = "Canopy temperature (Celsius)"
    df = data.frame(row.names=row.names(Temperature))
    df[["Mean"]] = Temperature$Tcan_mean
    df[["Minimum"]] = Temperature$Tcan_min
    df[["Maximum"]] = Temperature$Tcan_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="SoilTemperature") {
    if(is.null(ylab)) ylab = "Soil temperature (Celsius)"
    df = data.frame(row.names=row.names(Temperature))
    df[["Mean"]] = Temperature$Tsoil_mean
    df[["Minimum"]] = Temperature$Tsoil_min
    df[["Maximum"]] = Temperature$Tsoil_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
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
    df[["Latent heat vaporisation"]] = -EnergyBalance$LEVcan
    df[["Latent heat fusion"]] = -EnergyBalance$LEFsnow
    df[["Convection can./atm."]] = -EnergyBalance$Hcan
    df[["Convection soil/can."]] = -EnergyBalance$Hcansoil
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
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
    df[["Latent heat vaporisation"]] = -EnergyBalance$LEVsoil
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    g<-.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim)
    return(g)
  }
}
.plot_biomass<-function(BiomassBalance, type,  
                      dates = NULL, 
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                      summary.freq = NULL, ...) {
  df<-as.data.frame(BiomassBalance)
  names(df)<-c("Structural balance", "Labile balance","Plant individual balance", "Mortality loss",
               "Cohort balance")
  if(is.null(ylab))  ylab = expression(g%.%m^{-2})    
  if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
  if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
  return(.multiple_dynamics(as.matrix(df), ylab = ylab, ylim = ylim))
}
.plot_carbon<-function(CarbonBalance, type,  
                        dates = NULL, 
                        xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                        summary.freq = NULL, ...) {
  df<-as.data.frame(CarbonBalance)
  df[,2] = - df[,2]
  df[,3] = - df[,3]
  names(df)<-c("Gross primary production", "Maintenance respiration","Synthesis respiration", "Net primary production")
  if(is.null(ylab))  ylab = expression(gC%.%m^{-2})    
  if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),,drop = FALSE]
  if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
  return(.multiple_dynamics(as.matrix(df), ylab = ylab, ylim = ylim))
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
  
  if(is.null(cohorts)) {
    cohorts = row.names(input$cohorts)
  } else if(is.numeric(cohorts)){
    cohorts = row.names(input$cohorts)[cohorts]
  }
  spnames = as.character(input$cohorts[cohorts,"Name"])
  
  if(type=="PlantTranspiration") {
    m = .extractSubdaily(x, "E", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="TranspirationPerLeaf") {
    m = .extractSubdaily(x, "E", dates)[,c("datetime", cohorts), drop=FALSE]
    lai = .extractSubdaily(x, "PlantLAI", dates)[,c("datetime", cohorts), drop=FALSE]
    m[,-1] = m[,-1, drop = FALSE]/lai[,-1,drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantGrossPhotosynthesis") {
    m = .extractSubdaily(x, "Ag", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="GrossPhotosynthesisPerLeaf") {
    m = .extractSubdaily(x, "Ag", dates)[,c("datetime", cohorts), drop=FALSE]
    lai = .extractSubdaily(x, "PlantLAI", dates)[,c("datetime", cohorts), drop=FALSE]
    m[,-1] = m[,-1, drop = FALSE]/lai[,-1,drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantNetPhotosynthesis") {
    m = .extractSubdaily(x, "An", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="NetPhotosynthesisPerLeaf") {
    m = .extractSubdaily(x, "An", dates)[,c("datetime", cohorts), drop=FALSE]
    lai = .extractSubdaily(x, "PlantLAI", dates)[,c("datetime", cohorts), drop=FALSE]
    m[,-1] = m[,-1, drop = FALSE]/lai[,-1,drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsiAverage") {
    m = .extractSubdaily(x, "LeafPsi", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("LeafRWC", "StemRWC", "LeafPLC", "StemPLC",
                      "LeafSympRWC", "StemSympRWC")) {
    m = .extractSubdaily(x, type, dates)[,c("datetime", cohorts), drop=FALSE]
    m[, -1] <- 100*m[, -1]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("StemPsi", "RootPsi", "LeafSympPsi", "StemSympPsi")) {
    m = .extractSubdaily(x, type, dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantWaterBalance") {
    m = .extractSubdaily(x, "PWB", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="WaterBalancePerLeaf") {
    m = .extractSubdaily(x, "PWB", dates)[,c("datetime", cohorts), drop=FALSE]
    m[,-1] = m[,-1, drop = FALSE]/x$Plants$LAI
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="SoilPlantConductance") {
    m = .extractSubdaily(x, "dEdP", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="Temperature") {
    m = .extractSubdaily(x, "Temperature", dates)
    if("Tsoil.1" %in% colnames(m)) {
      m = m[,c("datetime","Tatm", "Tcan", "Tsoil.1")]
    } else {
      sm = .extractSubdaily(x, "SoilTemperature", dates)
      m = cbind(m, sm[,"1"])
    }
    names(m) = c("datetime", "Above-canopy","Inside-canopy", "Soil")
    if(is.null(ylab)) ylab = "Temperature (Celsius)"
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="CanopyEnergyBalance") {
    if(is.null(ylab)) ylab = expression(W%.%m^{-2})
    ceb = .extractSubdaily(x, "CanopyEnergyBalance", dates)
    seb = .extractSubdaily(x, "SoilEnergyBalance", dates)
    df = data.frame(datetime=ceb$datetime)
    df[["Balance"]] = ceb$Ebalcan
    df[["SWR abs."]] = ceb$SWRcan 
    df[["LWR net"]] = ceb$LWRcan
    df[["Latent heat vaporisation"]] = -ceb$LEVcan
    df[["Latent heat fusion"]] = -seb$LEFsnow
    df[["Convection can./atm."]] = -ceb$Hcan
    df[["Convection soil/can."]] = -seb$Hcansoil
    return(.multiple_dynamics_subdaily(df,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="SoilEnergyBalance") {
    if(is.null(ylab)) ylab = expression(W%.%m^{-2})    
    seb = .extractSubdaily(x, "SoilEnergyBalance", dates)
    df = data.frame(datetime=seb$datetime)
    df[["Balance"]] = seb$Ebalsoil
    df[["SWR abs."]] = seb$SWRsoil
    df[["LWR net"]] = seb$LWRsoil
    df[["Convection soil/can."]] = seb$Hcansoil
    df[["Latent heat vaporisation"]] = -seb$LEVsoil
    return(.multiple_dynamics_subdaily(df,  xlab = xlab, ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantExtraction") {
    m = .extractSubdaily(x, "ExtractionInst", dates)
    if(is.null(ylab)) ylab =.getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafLAI") {
    mSu = .extractSubdaily(x, "SunlitLeaves$LAI", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$LAI", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsi") {
    mSu = .extractSubdaily(x, "SunlitLeaves$Psi", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$Psi", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafAbsorbedSWR") {
    mSu = .extractSubdaily(x, "SunlitLeaves$Abs_SWR", dates)
    mSh = .extractSubdaily(x, "ShadeLeaves$Abs_SWR", dates)
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafAbsorbedPAR") {
    mSu = .extractSubdaily(x, "SunlitLeaves$Abs_PAR", dates)
    mSh = .extractSubdaily(x, "ShadeLeaves$Abs_PAR", dates)
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafNetLWR") {
    mSu = .extractSubdaily(x, "SunlitLeaves$Net_LWR", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$Net_LWR", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafTranspiration") {
    mSu = .extractSubdaily(x, "SunlitLeaves$E", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$E", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafGrossPhotosynthesis") {
    mSu = .extractSubdaily(x, "SunlitLeaves$Ag", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$Ag", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafNetPhotosynthesis") {
    mSu = .extractSubdaily(x, "SunlitLeaves$An", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$An", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafStomatalConductance") {
    mSu = .extractSubdaily(x, "SunlitLeaves$Gsw", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$Gsw", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafTemperature") {
    mSu = .extractSubdaily(x, "SunlitLeaves$Temp", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$Temp", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafVPD") {
    mSu = .extractSubdaily(x, "SunlitLeaves$VPD", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$VPD", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafCi") {
    mSu = .extractSubdaily(x, "SunlitLeaves$Ci", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$Ci", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafIntrinsicWUE") {
    mSu = .extractSubdaily(x, "SunlitLeaves$iWUE", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = .extractSubdaily(x, "ShadeLeaves$iWUE", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily_sunlit_shade(mSu, mSh, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("GrossPhotosynthesis", "MaintenanceRespiration", "GrowthCosts", "RootExudation" , "LabileCarbonBalance",
                      "SugarLeaf","StarchLeaf","SugarSapwood","StarchSapwood", "SugarTransport")) {
    m = .extractSubdaily(x, type, dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
}
