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
  spnames = as.character(input$cohorts[cohorts,"Name"])
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
  else if(type=="LAI") {
    return(.plot_stand(Stand = x$Stand,
                      type = type, dates = dates, 
                      xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                      summary.freq = summary.freq, ...))
  }
  else if(type %in% c("StemPLC", "StemRWC", "LeafRWC", "StemSympRWC", "LeafSympRWC")) {
    OM = x$Plants[[type]][,cohorts,drop=FALSE]*100
    return(.plot_plant_om(OM, PlantsLAI, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("PlantLAI","PlantTranspiration","PlantNetPhotosynthesis", "PlantGrossPhotosynthesis",
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
                      "PlantWaterBalance")) {
    if(type=="SoilPlantConductance") OM = Plants[["dEdP"]][,cohorts,drop=FALSE]
    else OM = x$Plants[[type]][,cohorts,drop=FALSE]
    return(.plot_plant_om(OM, PlantsLAI, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("LeafPsiMin_SL", "LeafPsiMax_SL", "GSWMin_SL", "GSWMax_SL", "TempMin_SL", "TempMax_SL")) {
    subType = strsplit(type,"_")[[1]][1]
    OM = x$SunlitLeaves[[subType]][,cohorts,drop=FALSE]
    return(.plot_plant_om(OM, PlantsLAI, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type %in% c("LeafPsiMin_SH", "LeafPsiMax_SH", "GSWMin_SH", "GSWMax_SH", "TempMin_SH", "TempMax_SH")) {
    subType = strsplit(type,"_")[[1]][1]
    OM = x$ShadeLeaves[[subType]][,cohorts,drop=FALSE]
    return(.plot_plant_om(OM, PlantsLAI, spnames,
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
    return(.plot_plant_om(df, PlantsLAI, spnames,
                   type, bySpecies = bySpecies, dates = dates, 
                   xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                   summary.freq = summary.freq, ...))
  } 
  else if(type == "LeafPsiRange") {
    OM1 = x$Plants$LeafPsiMax[,cohorts,drop=FALSE]
    OM2 = x$Plants$LeafPsiMin[,cohorts,drop=FALSE]
    if(bySpecies) {
      OM1 = .averageByLAISpecies(OM1, PlantsLAI, spnames)
      OM2 = .averageByLAISpecies(OM2, PlantsLAI, spnames)
    } 
    if(!is.null(dates)) {
      OM1 = OM1[row.names(OM1) %in% as.character(dates),]
      OM2 = OM2[row.names(OM2) %in% as.character(dates),]
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
    return(.plot_plant_om(OM, PlantsLAI, spnames,
                          type, bySpecies = bySpecies, dates = dates, 
                          xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                          summary.freq = summary.freq, ...))
  } 
  else if(type=="Temperature") {
    if(is.null(ylab)) ylab = "Temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Above-canopy"]] = x$Temperature$Tatm_mean
    df[["Inside-canopy"]] = x$Temperature$Tcan_mean
    df[["Soil"]] = x$Temperature$Tsoil_mean
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  }
  else if(type=="TemperatureRange") {
    if(is.null(ylab)) ylab = "Temperature (Celsius)"
    df1 = data.frame(row.names=row.names(x$Temperature))
    df1[["Above-canopy"]] = x$Temperature$Tatm_min
    df1[["Inside-canopy"]] = x$Temperature$Tcan_min
    df1[["Soil"]] = x$Temperature$Tsoil_min
    df2 = data.frame(row.names=row.names(x$Temperature))
    df2[["Above-canopy"]] = x$Temperature$Tatm_max
    df2[["Inside-canopy"]] = x$Temperature$Tcan_max
    df2[["Soil"]] = x$Temperature$Tsoil_max
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
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tatm_mean
    df[["Minimum"]] = x$Temperature$Tatm_min
    df[["Maximum"]] = x$Temperature$Tatm_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="CanopyTemperature") {
    if(is.null(ylab)) ylab = "Canopy temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tcan_mean
    df[["Minimum"]] = x$Temperature$Tcan_min
    df[["Maximum"]] = x$Temperature$Tcan_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="SoilTemperature") {
    if(is.null(ylab)) ylab = "Soil temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tsoil_mean
    df[["Minimum"]] = x$Temperature$Tsoil_min
    df[["Maximum"]] = x$Temperature$Tsoil_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="CanopyEnergyBalance") {
    if(is.null(ylab)) ylab = expression(MJ%.%m^{-2})    
    df = data.frame(row.names=row.names(x$EnergyBalance))
    df[["Balance"]] = x$EnergyBalance$Ebalcan
    df[["SWR abs."]] = x$EnergyBalance$SWRcan 
    df[["Net LWR"]] = x$EnergyBalance$LWRcan
    df[["Latent heat"]] = -x$EnergyBalance$LEcan
    df[["Convection can./atm."]] = -x$EnergyBalance$Hcan
    df[["Convection soil/can."]] = -x$EnergyBalance$Hcansoil
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    g<-.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim)
    return(g)
  } 
  else if(type=="SoilEnergyBalance") {
    if(is.null(ylab)) ylab = expression(MJ%.%m^{-2})    
    df = data.frame(row.names=row.names(x$EnergyBalance))
    df[["Balance"]] = x$EnergyBalance$Ebalsoil
    df[["SWR abs."]] = x$EnergyBalance$SWRsoil
    df[["Net LWR"]] = x$EnergyBalance$LWRsoil
    df[["Convection soil/can."]] = x$EnergyBalance$Hcansoil
    df[["Latent heat"]] = -x$EnergyBalance$LEsoil
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(!is.null(summary.freq)) df = .temporalSummary(df, summary.freq, mean, na.rm=TRUE)
    g<-.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim)
    return(g)
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
  spnames = as.character(input$cohorts[cohorts,"Name"])
  PlantsLAI = x$Plants$LAI[,cohorts, drop=FALSE]
  
  if(type %in% TYPES_SWB) {
    plot.spwb(x,type, cohorts, bySpecies, dates, subdaily, xlim, ylim, xlab, ylab, 
              summary.freq, ...)
  } 
  else if(type %in% c("GrossPhotosynthesis", "MaintenanceRespiration",  "GrowthCosts", 
                      "CarbonBalance", 
                      "SugarLeaf","StarchLeaf","SugarSapwood","StarchSapwood", "SugarTransport", "RootExudation",
                      "LeafPI0", "StemPI0")) {
      OM = x$PlantCarbonBalance[[type]][,cohorts,drop=FALSE]
  } 
  else if(type %in% c("SapwoodArea", "LeafArea", "FineRootArea", 
                      "SapwoodBiomass", "LeafBiomass", "FineRootBiomass",
                      "LabileBiomass", "TotalLivingBiomass")) {
    OM = x$PlantStructure[[type]][,cohorts,drop=FALSE]
  } 
  else if(type %in% c("SAgrowth", "LAgrowth", "FRAgrowth")) {
    OM = x$PlantGrowth[[type]][,cohorts,drop=FALSE]
  } 
  else if(type=="HuberValue") {
    OM = x$PlantStructure[["SapwoodArea"]][,cohorts,drop=FALSE] / x$PlantStructure[["LeafArea"]][,cohorts,drop=FALSE]
  } 
  else if(type=="RootAreaLeafArea") {
    OM = x$PlantStructure[["FineRootArea"]][,cohorts,drop=FALSE] / x$PlantStructure[["LeafArea"]][,cohorts,drop=FALSE]
  } 
  return(.plot_plant_om(OM, PlantsLAI, spnames,
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
  TYPES_FORDYN_UNIQUE =   c("StandBasalArea", "StandLAI", "StandDensity",
              "SpeciesBasalArea", "SpeciesLAI", "SpeciesDensity",
              "CohortBasalArea", "CohortLAI", "CohortDensity")
  type = match.arg(type,c(TYPES_GROWTH, TYPES_FORDYN))  
  
  if(type %in% TYPES_GROWTH) {
    
    PlantsLAI = summary(x, freq = "days", output = "Plants$LAI")
    if(is.null(cohorts))  cohorts = row.names(x$GrowthResults[[1]]$growthInput$cohorts)
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
    else if(type=="LAI") {
      Stand = summary(x, freq = "days", output = "Stand")
      return(.plot_stand(Stand = Stand,
                         type = type, dates = dates, 
                         xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab,
                         summary.freq = summary.freq,...))
    }
    else if(type %in% TYPES_GROWTH_UNIQUE) {
      if(type %in% c("GrossPhotosynthesis", "MaintenanceRespiration",  "GrowthCosts", 
                     "CarbonBalance", 
                     "SugarLeaf","StarchLeaf","SugarSapwood","StarchSapwood", "SugarTransport", "RootExudation",
                     "LeafPI0", "StemPI0")) {
        OM = summary(x, freq = "days", output = paste0("PlantCarbonBalance$",type))[,cohorts,drop=FALSE]
      } 
      else if(type %in% c("SapwoodArea", "LeafArea", "FineRootArea", 
                          "SapwoodBiomass", "LeafBiomass", "FineRootBiomass",
                          "LabileBiomass", "TotalLivingBiomass")) {
        OM = summary(x, freq = "days", output = paste0("PlantStructure$",type))[,cohorts,drop=FALSE]
      } 
      else if(type %in% c("SAgrowth", "LAgrowth", "FRAgrowth")) {
        OM = summary(x, freq = "days", output = paste0("PlantGrowth$",type))[,cohorts,drop=FALSE]
      } 
      else if(type %in% c("HuberValue")) {
        SA = summary(x, freq = "days", output = "PlantStructure$SapwoodArea")[,cohorts,drop=FALSE]
        LA = summary(x, freq = "days", output = "PlantStructure$LeafArea")[,cohorts,drop=FALSE]
        OM = SA/LA
      } 
      return(.plot_plant_om(OM, PlantsLAI, spnames,
                            type, bySpecies = bySpecies, dates = dates, 
                            xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                            summary.freq = summary.freq, ...))
    }
  }


  ## FORDYN PLOT
  i_type = which(TYPES_FORDYN_UNIQUE %in% type)
  
  vars = rep(c("TreeBasalAreaLive", "LeafAreaIndex", "TreeDensityLive"),3)
  tables = c(rep("StandSummary",3),rep("SpeciesSummary",3),rep("CohortSummary",3))
  
  if(is.null(ylab)) ylab = .getYLab(type)
  if(is.null(xlab)) xlab = "Step"
  
  out = x[[tables[i_type]]]
  df = data.frame(Step = out[["Step"]], 
                  y = out[[vars[i_type]]])
  if(tables[i_type]=="SpeciesSummary") df$group = as.character(out[["Species"]])
  else if(tables[i_type]=="CohortSummary") df$group = as.character(out[["Cohort"]])
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

