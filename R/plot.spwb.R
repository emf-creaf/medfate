plot.spwb<-function(x, type="PET_Precipitation", cohorts = NULL, bySpecies = FALSE,
                    dates = NULL, subdaily = FALSE, 
                    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                    summary.freq = NULL, ...) {
  
  if(subdaily) return(.plotsubdaily(x,type, cohorts, bySpecies, dates, 
                                    xlim, ylim, xlab, ylab))
  
  if(("spwb" %in% class(x)) || ("pwb" %in% class(x))) {
    input = x$spwbInput
  } else {
    input = x$growthInput
  }
  TYPES = .getDailySPWBPlotTypes(input$control$transpirationMode)  

  type = match.arg(type,TYPES)  
  if(is.null(xlab)) xlab = ""
  if(type %in% c("PET_Precipitation", "PET_NetRain", "Evapotranspiration", "Snow",
                 "WTD", "Export", "SoilVol")) {
    return(.plot_wb(WaterBalance = x$WaterBalance, Soil = x$Soil, input_soil = input$soil, 
                    type = type, dates = dates, 
                    xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                    summary.freq = summary.freq, ...))
  } 
  return(plot.pwb(x, type=type, cohorts = cohorts, bySpecies = bySpecies,
                  dates = dates, xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                  summary.freq = summary.freq, ...))
}

plot.pwb<-function(x, type="PlantTranspiration", cohorts = NULL, bySpecies = FALSE,
                   dates = NULL, subdaily = FALSE,
                   xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                   summary.freq = NULL,...) {
  
  if(subdaily) return(.plotsubdaily(x,type, cohorts, bySpecies, dates, 
                                    xlim, ylim, xlab, ylab))
  
  if(("spwb" %in% class(x)) || ("pwb" %in% class(x))) {
    input = x$spwbInput
  } else {
    input = x$growthInput
  }
  nlayers = length(input$soil$W)
  transpMode = input$control$transpirationMode
  
  if(is.null(cohorts))  cohorts = row.names(input$cohorts)
  else if(length(cohorts)==1) {
    if(cohorts=="T") {
      cohorts = row.names(input$cohorts)
      cohorts = cohorts[substr(cohorts,1,1)=="T"]
    } else if(cohorts=="S") {
      cohorts = row.names(input$cohorts)
      cohorts = cohorts[substr(cohorts,1,1)=="S"]
    }
  }

  spnames = as.character(input$cohorts[cohorts,"Name"])
  
  PlantsLAIlive = x$Plants$LAIlive[,cohorts, drop=FALSE]
  PlantsLAI = x$Plants$LAI[,cohorts, drop=FALSE]
  
  TYPES = .getDailyPWBPlotTypes(transpMode)
  type = match.arg(type,TYPES)  
  if(is.null(xlab)) xlab = ""  
  
  if(type %in% c("SoilPsi", "SoilTheta", "SoilRWC", "PlantExtraction", "HydraulicRedistribution")) {
    return(.plot_soil(Soil = x$Soil, input_soil = input$soil, input_control = input$control,
                    type = type, dates = dates, 
                    xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                    summary.freq = summary.freq, ...))
  }
  else if(type %in% c("LAI", "GroundIrradiance")) {
    return(.plot_stand(Stand = x$Stand,
                      type = type, dates = dates, 
                      xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                      summary.freq = summary.freq, ...))
  }
  else if(type %in% c("FPAR", "AbsorbedSWRFraction")) {
    OM = x$Plants[[type]][,cohorts,drop=FALSE]
    return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                          type, bySpecies = bySpecies, dates = dates, 
                          xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                          summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("StemPLC", "StemRWC", "LeafRWC")) {
    OM = x$Plants[[type]][,cohorts,drop=FALSE]*100
    return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("PlantLAI","PlantLAIlive","PlantTranspiration","PlantNetPhotosynthesis", "PlantGrossPhotosynthesis",
                      "PlantAbsorbedSWR","PlantNetLWR")) {
    subtype = substr(type,6,nchar(type))
    OM = x$Plants[[subtype]][,cohorts,drop=FALSE]
    return(.plot_plant_om_sum(OM, spnames,
                          type, bySpecies = bySpecies, dates = dates, 
                          xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                          summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("SoilPlantConductance","PlantPsi", "LeafPsiMin",
                      "LeafPsiMax", "StemPsi", "RootPsi", "PlantStress",
                      "PlantWaterBalance", "LFMC")) {
    if(type=="SoilPlantConductance") OM = x$Plants[["dEdP"]][,cohorts,drop=FALSE]
    else OM = x$Plants[[type]][,cohorts,drop=FALSE]
    return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("LeafPsiMin_SL", "LeafPsiMax_SL", "GSWMin_SL", "GSWMax_SL", "TempMin_SL", "TempMax_SL")) {
    subType = strsplit(type,"_")[[1]][1]
    OM = x$SunlitLeaves[[subType]][,cohorts,drop=FALSE]
    return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("LeafPsiMin_SH", "LeafPsiMax_SH", "GSWMin_SH", "GSWMax_SH", "TempMin_SH", "TempMax_SH")) {
    subType = strsplit(type,"_")[[1]][1]
    OM = x$ShadeLeaves[[subType]][,cohorts,drop=FALSE]
    return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("TranspirationPerLeaf","GrossPhotosynthesisPerLeaf", "NetPhotosynthesisPerLeaf",
                      "AbsorbedSWRPerLeaf", "NetLWRPerLeaf")) {
    subtype = substr(type, 1, nchar(type)-7)
    df = x$Plants[[subtype]][,cohorts,drop=FALSE]
    df = df/PlantsLAI
    df[PlantsLAI==0] = NA
    return(.plot_plant_om(df, PlantsLAIlive, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type == "LeafPsiRange") {
    OM1 = x$Plants$LeafPsiMax[,cohorts,drop=FALSE]
    OM2 = x$Plants$LeafPsiMin[,cohorts,drop=FALSE]
    if(bySpecies) {
      OM1 = .averageByLAISpecies(OM1, PlantsLAIlive, spnames)
      OM2 = .averageByLAISpecies(OM2, PlantsLAIlive, spnames)
    } 
    if(!is.null(dates)) {
      OM1 = OM1[row.names(OM1) %in% as.character(dates),,drop = FALSE]
      OM2 = OM2[row.names(OM2) %in% as.character(dates),,drop = FALSE]
    }
    if(!is.null(summary.freq)) {
      OM1 = .temporalSummary(OM1, summary.freq, mean, na.rm=TRUE)
      OM2 = .temporalSummary(OM2, summary.freq, mean, na.rm=TRUE)
    }
    return(.multiple_dynamics_range(as.matrix(OM1), as.matrix(OM2),  xlab = xlab, ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantWUE") {
    OM = x$Plants$GrossPhotosynthesis/x$Plants$Transpiration
    OM[x$Plants$Transpiration==0] = 0
    OM = OM[,cohorts,drop=FALSE]
    return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                          type, bySpecies = bySpecies, dates = dates, 
                          xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                          summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("Temperature","TemperatureRange", "AirTemperature",
                      "CanopyTemperature", "SoilTemperature")) {
    return(.plot_temperature(x$Temperature, type,  
                                dates = dates, 
                                xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                                summary.freq = summary.freq, ...))
  }
  else if(type %in% c("CanopyEnergyBalance", "SoilEnergyBalance")) {
    return(.plot_energybalance(x$EnergyBalance, type,  
                               dates = dates, 
                               xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                               summary.freq = summary.freq, ...))
  } 
}
plot.growth<-function(x, type="PET_Precipitation", cohorts = NULL, bySpecies = FALSE, 
                      dates = NULL, subdaily = FALSE, 
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                      summary.freq = NULL, ...) {
  
  if(subdaily) return(.plotsubdaily(x,type, cohorts, bySpecies, dates, 
                                    xlim, ylim, xlab, ylab))
  
  # Get common elements
  input = x$growthInput
  nlayers = length(input$soil$W)
  transpMode = input$control$transpirationMode
  
  TYPES_GROWTH = .getDailyGROWTHPlotTypes(transpMode)
  TYPES_SWB = .getDailySPWBPlotTypes(transpMode)  

  type = match.arg(type,TYPES_GROWTH)  
  
  if(is.null(cohorts))  cohorts = row.names(input$cohorts)
  else if(length(cohorts)==1) {
    if(cohorts=="T") {
      cohorts = row.names(input$cohorts)
      cohorts = cohorts[substr(cohorts,1,1)=="T"]
    } else if(cohorts=="S") {
      cohorts = row.names(input$cohorts)
      cohorts = cohorts[substr(cohorts,1,1)=="S"]
    }
  }

  spnames = as.character(input$cohorts[cohorts,"Name"])
  PlantsLAIlive = x$Plants$LAIlive[,cohorts, drop=FALSE]
  PlantsLAI = x$Plants$LAI[,cohorts, drop=FALSE]
  
  if(type %in% TYPES_SWB) {
    return(plot.spwb(x,type, cohorts, bySpecies, dates, subdaily, xlim, ylim, xlab, ylab, 
              summary.freq, ...))
  } 
  else if(type == "BiomassBalance") {
    return(.plot_biomass(BiomassBalance= x$BiomassBalance,
                       type = type, dates = dates, 
                       xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                       summary.freq = summary.freq, ...))
  }
  else if(type == "CarbonBalance") {
    return(.plot_carbon(CarbonBalance= x$CarbonBalance,
                         type = type, dates = dates, 
                         xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                         summary.freq = summary.freq, ...))
  }
  else if(type %in% c("GrossPhotosynthesis", "MaintenanceRespiration",  "GrowthCosts", 
                      "LabileCarbonBalance", 
                      "SugarLeaf","StarchLeaf","SugarSapwood","StarchSapwood", "SugarTransport", "RootExudation")) {
      OM = x$LabileCarbonBalance[[type]][,cohorts,drop=FALSE]
  } 
  else if(type == "PhotosynthesisMaintenanceRatio") {
    OM1 = x$LabileCarbonBalance[["GrossPhotosynthesis"]][,cohorts,drop=FALSE]
    OM2 = x$LabileCarbonBalance[["MaintenanceRespiration"]][,cohorts,drop=FALSE]
    OM = OM1/OM2
  } 
  else if(type %in% c("StructuralBiomassBalance","LabileBiomassBalance", "PlantBiomassBalance",
                      "MortalityBiomassLoss","CohortBiomassBalance")) {
    OM = x$PlantBiomassBalance[[type]][,cohorts,drop=FALSE]
  }
  else if(type %in% c("SapwoodBiomass", "LeafBiomass", "FineRootBiomass", 
                      "SapwoodArea", "LeafArea", "FineRootArea", 
                      "DBH", "Height", "HuberValue", "RootAreaLeafArea")) {
    OM = x$PlantStructure[[type]][,cohorts,drop=FALSE]
  } 
  else if(type %in% c("SAgrowth", "LAgrowth", "FRAgrowth", "StarvationRate", "DessicationRate", "MortalityRate")) {
    OM = x$GrowthMortality[[type]][,cohorts,drop=FALSE]
  } 
  return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                        type, bySpecies = bySpecies, dates = dates, 
                        xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                        summary.freq = summary.freq, ...))
}

plot.fordyn<-function(x, type="StandBasalArea",  
                      cohorts = NULL, bySpecies = FALSE, dates = NULL, 
                      xlim = NULL, ylim=NULL, xlab = NULL, ylab=NULL, 
                      summary.freq = NULL,...) {
  
  input_control = x$GrowthResults[[1]]$growthInput$control
  TYPES_GROWTH = .getDailyGROWTHPlotTypes(input_control$transpirationMode)
  TYPES_GROWTH_UNIQUE = .getUniqueDailyGROWTHPlotTypes(input_control$transpirationMode)
  TYPES_FORDYN_UNIQUE = .getUniqueFORDYNPlotTypes(input_control$transpirationMode)  
  
  type = match.arg(type,c(TYPES_GROWTH, TYPES_FORDYN_UNIQUE))  
  
  if(type %in% TYPES_GROWTH) {
    
    if(is.null(cohorts))  cohorts = row.names(x$GrowthResults[[1]]$growthInput$cohorts)
    PlantsLAI = summary(x, freq = "days", output = "Plants$LAI")[,cohorts,drop=FALSE]
    PlantsLAIlive = summary(x, freq = "days", output = "Plants$LAIlive")[,cohorts,drop=FALSE]
    spnames = as.character(x$GrowthResults[[1]]$growthInput$cohorts[cohorts,"Name"])
    input_soil = x$GrowthResults[[1]]$growthInput$soil
    
    if(type %in% c("PET_Precipitation", "PET_NetRain", "Evapotranspiration" ,"Snow",
                   "WTD", "Export", "SoilVol")) {
      input_soil = x$GrowthResults[[1]]$growthInput$soil
      WaterBalance = summary(x, freq = "days", output = "WaterBalance")
      Soil = summary(x, freq = "days", output = "Soil")
      return(.plot_wb(WaterBalance = WaterBalance, Soil = Soil, input_soil = input_soil,
                      type=type, dates = dates, ylim = ylim, xlab = xlab, ylab = ylab,
                      summary.freq = summary.freq,...))
    }
    else if(type %in% c("SoilPsi", "SoilTheta", "SoilRWC", "PlantExtraction", "HydraulicRedistribution")) {
      Soil = summary(x, freq = "days", output = "Soil")
      return(.plot_soil(Soil = Soil, input_soil = input_soil, input_control = input_control,
                        type=type, dates = dates, 
                        xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                        summary.freq = summary.freq, ...))
    }
    else if(type %in% c("LAI", "GroundIrradiance")) {
      Stand = summary(x, freq = "days", output = "Stand")
      return(.plot_stand(Stand = Stand,
                         type = type, dates = dates, 
                         xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab,
                         summary.freq = summary.freq,...))
    }
    else if(type %in% c("StemPLC", "StemRWC", "LeafRWC", "StemSympRWC", "LeafSympRWC")) {
      OM = summary(x, freq = "days", output = paste0("Plants$",type))[,cohorts,drop=FALSE]*100
      return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                            type, bySpecies = bySpecies, dates = dates, 
                            xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                            summary.freq = summary.freq, ...))
    } 
    else if(type %in% c("PlantLAI","PlantLAIlive","PlantTranspiration","PlantNetPhotosynthesis", "PlantGrossPhotosynthesis",
                        "PlantAbsorbedSWR","PlantNetLWR")) {
      subtype = substr(type,6,nchar(type))
      OM = summary(x, freq = "days", output = paste0("Plants$",subtype))[,cohorts,drop=FALSE]
      return(.plot_plant_om_sum(OM, spnames,
                                type, bySpecies = bySpecies, dates = dates, 
                                xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                                summary.freq = summary.freq, ...))
    } 
    else if(type=="PlantWUE") {
      GP = summary(x, freq = "days", output = "Plants$GrossPhotosynthesis")[,cohorts,drop=FALSE]
      E = summary(x, freq = "days", output = "Plants$Transpiration")[,cohorts,drop=FALSE]
      OM = GP/E
      OM[E==0] = 0
      return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                            type, bySpecies = bySpecies, dates = dates, 
                            xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                            summary.freq = summary.freq, ...))
    } 
    else if(type %in% c("SoilPlantConductance","PlantPsi", "LeafPsiMin",
                        "LeafPsiMax", "StemPsi", "RootPsi", "PlantStress",
                        "PlantWaterBalance", "FPAR", "AbsorbedSWRFraction")) {
      if(type=="SoilPlantConductance") OM = summary(x, freq = "days", output = "Plants$dEdP")[,cohorts,drop=FALSE]
      else OM = summary(x, freq = "days", output = paste0("Plants$",type))[,cohorts,drop=FALSE]
      return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                            type, bySpecies = bySpecies, dates = dates, 
                            xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                            summary.freq = summary.freq, ...))
    } 
    else if(type == "LeafPsiRange") {
      OM1 = summary(x, freq = "days", output = "Plants$LeafPsiMax")[,cohorts,drop=FALSE]
      OM2 = summary(x, freq = "days", output = "Plants$LeafPsiMin")[,cohorts,drop=FALSE]
      if(bySpecies) {
        OM1 = .averageByLAISpecies(OM1, PlantsLAIlive, spnames)
        OM2 = .averageByLAISpecies(OM2, PlantsLAIlive, spnames)
      } 
      if(!is.null(dates)) {
        OM1 = OM1[row.names(OM1) %in% as.character(dates),,drop = FALSE]
        OM2 = OM2[row.names(OM2) %in% as.character(dates),,drop = FALSE]
      }
      if(!is.null(summary.freq)) {
        OM1 = .temporalSummary(OM1, summary.freq, mean, na.rm=TRUE)
        OM2 = .temporalSummary(OM2, summary.freq, mean, na.rm=TRUE)
      }
      return(.multiple_dynamics_range(as.matrix(OM1), as.matrix(OM2),  xlab = xlab, ylab = ylab, ylim = ylim))
    }
    else if(type %in% c("TranspirationPerLeaf","GrossPhotosynthesisPerLeaf", "NetPhotosynthesisPerLeaf",
                        "AbsorbedSWRPerLeaf", "NetLWRPerLeaf")) {
      subtype = substr(type, 1, nchar(type)-7)
      OM = summary(x, freq = "days", output = paste0("Plants$",subtype))[,cohorts,drop=FALSE]
      OM = OM/PlantsLAI
      OM[PlantsLAI==0] = NA
      return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                            type, bySpecies = bySpecies, dates = dates, 
                            xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                            summary.freq = summary.freq, ...))
    } 
    else if(type %in% c("Temperature","TemperatureRange", "AirTemperature",
                        "CanopyTemperature", "SoilTemperature")) {
      OM = summary(x, freq = "days", output = "Temperature")
      return(.plot_temperature(OM, type,  
                               dates = dates, 
                               xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                               summary.freq = summary.freq, ...))
    }
    else if(type %in% c("CanopyEnergyBalance", "SoilEnergyBalance")) {
      OM = summary(x, freq = "days", output = "EnergyBalance")
      return(.plot_energybalance(OM, type,  
                                 dates = dates, 
                                 xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                                 summary.freq = summary.freq, ...))
    } 
    else if(type %in% TYPES_GROWTH_UNIQUE) {
      if(type=="BiomassBalance") {
        OM = summary(x, freq = "days", output = "BiomassBalance")
        return(.plot_biomass(OM, type,  
                             dates = dates, 
                             xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                             summary.freq = summary.freq, ...))
      }
      if(type=="CarbonBalance") {
        OM = summary(x, freq = "days", output = "CarbonBalance")
        return(.plot_carbon(OM, type,  
                             dates = dates, 
                             xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                             summary.freq = summary.freq, ...))
      }
      if(type %in% c("GrossPhotosynthesis", "MaintenanceRespiration",  "GrowthCosts", 
                     "LabileCarbonBalance", 
                     "SugarLeaf","StarchLeaf","SugarSapwood","StarchSapwood", "SugarTransport", "RootExudation")) {
        OM = summary(x, freq = "days", output = paste0("LabileCarbonBalance$",type))[,cohorts,drop=FALSE]
      } 
      if(type =="PhotosynthesisMaintenanceRatio") {
        OM1 = summary(x, freq = "days", output = "LabileCarbonBalance$GrossPhotosynthesis")[,cohorts,drop=FALSE]
        OM2 = summary(x, freq = "days", output = "LabileCarbonBalance$MaintenanceRespiration")[,cohorts,drop=FALSE]
        OM = OM1/OM2
      }
      else if(type %in% c("StructuralBiomassBalance", "LabileBiomassBalance", "PlantBiomassBalance", 
                          "MortalityBiomassLoss", "CohortBiomassBalance")) {
        OM = summary(x, freq = "days", output = paste0("PlantBiomassBalance$",type))[,cohorts,drop=FALSE]
      } 
      else if(type %in% c("LeafBiomass", "SapwoodBiomass", "FineRootBiomass", "SapwoodArea", "LeafArea", "FineRootArea","DBH", "Height",
                          "HuberValue", "RootAreaLeafArea")) {
        OM = summary(x, freq = "days", output = paste0("PlantStructure$",type))[,cohorts,drop=FALSE]
      } 
      else if(type %in% c("SAgrowth", "LAgrowth", "FRAgrowth", "StarvationRate", "MortalityRate", "DessicationRate")) {
        OM = summary(x, freq = "days", output = paste0("GrowthMortality$",type))[,cohorts,drop=FALSE]
      } 
      return(.plot_plant_om(OM, PlantsLAIlive, spnames,
                            type, bySpecies = bySpecies, dates = dates, 
                            xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                            summary.freq = summary.freq, ...))
    }
  }


  ## FORDYN PLOT

  if(is.null(ylab)) ylab = .getYLab(type)
  if(is.null(xlab)) xlab = "Step"

  if(type %in% c("NumTreeSpecies", "NumTreeCohorts", "NumShrubSpecies", "NumShrubCohorts"))  {
    out = x$StandSummary
    var = type
  }
  else if(type == "StandDensity")  {
    out = x$StandSummary
    var = "TreeDensityLive"
  }
  else if(type == "StandBasalArea")  {
    out = x$StandSummary
    var = "TreeBasalAreaLive"
  }
  else if(type == "StandShrubCover")  {
    out = x$StandSummary
    var = "ShrubCoverLive"
  }
  else if(type %in% c("DominantTreeHeight", "DominantTreeDiameter", "QuadraticMeanTreeDiameter", "HartBeckingIndex")) {
    out = x$StandSummary
    var = type
  }
  else if(type == "SpeciesDensity")  {
    out = x$SpeciesSummary
    var = "TreeDensityLive"
  }
  else if(type == "SpeciesBasalArea")  {
    out = x$SpeciesSummary
    var = "TreeBasalAreaLive"
  } 
  else if(type == "SpeciesShrubCover")  {
    out = x$SpeciesSummary
    var = "ShrubCoverLive"
  } 
  else if(type == "NumCohortsSpecies")  {
    out = x$SpeciesSummary
    var = "NumCohorts"
  } 
  else if(type == "CohortDensity")  {
    out = x$CohortSummary
    var = "TreeDensityLive"
  }
  else if(type == "CohortBasalArea")  {
    out = x$CohortSummary
    var = "TreeBasalAreaLive"
  }
  else if(type == "CohortShrubCover")  {
    out = x$CohortSummary
    var = "ShrubCoverLive"
  }
  
  df = data.frame(Step = out[["Step"]], y = out[[var]])
  if(type %in% c("SpeciesDensity", "SpeciesBasalArea", "SpeciesShrubCover", "NumCohortsSpecies")) df$group = as.character(out[["Name"]])
  else if(type %in% c("CohortBasalArea", "CohortDensity", "CohortShrubCover")) {
    df$group = paste0(as.character(out[["Cohort"]]), " (", as.character(out[["Name"]]),")")
  }
  df = df[!is.na(df$y),]
  g<-ggplot(df, aes_string(x="Step", y="y"))
  if("group" %in% names(df)) {
    g <- g+ geom_line(aes_string(col="group"))+
      scale_color_discrete(name="")
  } else {
    g <- g+ geom_line()
  }
  if(!is.null(ylim)) g <- g+ylim(ylim)
  g<-g+theme_bw()+ylab(ylab)+xlab(xlab)
  return(g)
}

