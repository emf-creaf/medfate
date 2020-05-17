
.averageBySpecies<-function(OM, spnames) {
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
plot.spwb<-function(x, type="PET_Precipitation", cohorts = NULL, bySpecies = FALSE,
                    dates = NULL, subdaily = FALSE, 
                    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,...) {
  
  if(subdaily) return(.plotsubdaily(x,type, cohorts, bySpecies, dates, 
                                    xlim, ylim, xlab, ylab))
  
  if(("spwb" %in% class(x)) || ("pwb" %in% class(x))) {
    input = x$spwbInput
  } else {
    input = x$growthInput
  }
  
  soilInput = x$soilInput
  WaterBalance = x$WaterBalance
  Soil = x$Soil
  Stand = x$Stand
  Plants = x$Plants
  
  nlayers = length(soilInput$W)
  
  transpMode = input$control$transpirationMode

  TYPES = .getDailySPWBPlotTypes(transpMode)  

  type = match.arg(type,TYPES)  
  if(is.null(xlab)) xlab = ""
  if(type=="PET_Precipitation") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})
    if(is.null(ylab)) ylab = expression(L%.%m^{-2}) 
    df = data.frame(row.names=row.names(WaterBalance))
    df[["PET"]] = WaterBalance$PET
    df[["Precipitation"]] = WaterBalance$Precipitation
    df[["Snow"]] = WaterBalance$Snow
    df[["Date"]] = as.Date(row.names(WaterBalance))
    if(!is.null(dates)) df = df[df$Date %in% dates,]
    Date = df$Date
    Snow = df$Snow
    Precipitation = df$Precipitation
    g<-ggplot(df)+
      geom_area(aes(x=Date, y=Precipitation, fill="Precipitation"))+
      geom_area(aes(x=Date, y=Snow, fill="Snow"))+
      geom_path(aes(x=Date, y=PET, col="PET"))+
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
    Date = df$Date
    PET = df$PET
    NetRain = df$NetRain
    g<-ggplot(df)+
      geom_area( aes(x=Date, y=NetRain, fill="NetRain"))+
      geom_path(aes(x=Date, y=PET, col="PET"))+
      scale_fill_manual(name="", values=c("NetRain"="black"))+
      scale_color_manual(name="", values=c("PET"="gray"))+
      ylab(ylab)+xlab(xlab)+
      theme_bw()
    return(g)
  } 
  else if(type=="Snow") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})  
    df = data.frame(row.names=row.names(WaterBalance))
    df[["Snow"]] = WaterBalance$Snow
    df[["Snowpack"]] = Soil$SWE
    df[["Date"]] = as.Date(row.names(WaterBalance))
    if(!is.null(dates)) df = df[df$Date %in% dates,]
    Date = df$Date
    Snow = df$Swow
    Snowpack = df$Snowpack
    g<-ggplot(df)+
      geom_area( aes(x=Date, y=Snow, fill="Snow"))+
      geom_path(aes(x=Date, y=Snowpack, col="Snowpack"))+
      scale_fill_manual(name="", values=c("Snow"="black"))+
      scale_color_manual(name="", values=c("Snowpack"="gray"))+
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
    return(.multiple_dynamics(as.matrix(df), ylab=ylab, xlab=xlab, ylim = ylim))
  } 
  else if(type=="WTD") {
    if(is.null(ylab)) ylab = expression(paste("Water table depth  (mm)"))
    xv = Soil$WTD
    names(xv) = row.names(Soil)
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
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
    Date = df$Date
    Export = df$Export
    DeepDrainage = df$DeepDrainage
    Runoff = df$Runoff
    g<-ggplot(df)+
      geom_line( aes(x=Date, y=Export, col="Export"))+
      geom_line( aes(x=Date, y=DeepDrainage, col="Deep drainage"))+
      geom_line( aes(x=Date, y=Runoff, col="Runoff"))+
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
    return(.multiple_dynamics(as.matrix(MLM), ylab = ylab, ylim = ylim,
                              xlab=xlab, labels = c("Total", paste("Layer", 1:nlayers))))
  } 
  else {
    plot.pwb(x, type=type, cohorts = cohorts, bySpecies = bySpecies,
             dates = dates, xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  }
}

plot.pwb<-function(x, type="PlantTranspiration", cohorts = NULL, bySpecies = FALSE,
                   dates = NULL, subdaily = FALSE,
                   xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...) {
  
  if(subdaily) return(.plotsubdaily(x,type, cohorts, bySpecies, dates, 
                                    xlim, ylim, xlab, ylab))
  
  if(("spwb" %in% class(x)) || ("pwb" %in% class(x))) {
    input = x$spwbInput
  } else {
    input = x$growthInput
  }
  soilInput = x$soilInput
  
  WaterBalance = x$WaterBalance
  Soil = x$Soil
  Stand = x$Stand
  Plants = x$Plants
  SunlitLeaves = x$SunlitLeaves
  ShadeLeaves = x$ShadeLeaves
  
  nlayers = length(soilInput$W)
  
  transpMode = input$control$transpirationMode
  
  if(is.null(cohorts))  cohorts = row.names(input$cohorts)
  spnames = as.character(input$cohorts[cohorts,"Name"])
  PlantsLAI = Plants$LAI[,cohorts, drop=FALSE]
  
  TYPES = .getDailyPWBPlotTypes(transpMode)
  type = match.arg(type,TYPES)  
  if(is.null(xlab)) xlab = ""  
  
  if(type=="SoilPsi") {
    PsiM = Soil[,paste("psi",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Soil water potential (MPa)"    
    if(!is.null(dates)) PsiM = PsiM[row.names(PsiM) %in% as.character(dates),]
    return(.multiple_dynamics(as.matrix(PsiM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="SoilTheta") {
    WM = Soil[,paste("W",1:nlayers,sep=".")]
    if(!is.null(dates)) WM = WM[row.names(WM) %in% as.character(dates),]
    theta_FC = soil_thetaFC(soilInput, model = input$control$soilFunctions)
    WM = 100*sweep(WM, 2,theta_FC, "*")
    if(is.null(ylab)) ylab = "Soil moisture (% volume)"
    return(.multiple_dynamics(as.matrix(WM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="SoilRWC") {
    WM = Soil[,paste("W",1:nlayers,sep=".")]
    if(!is.null(dates)) WM = WM[row.names(WM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = "Soil moisture (% field capacity)"
    return(.multiple_dynamics(as.matrix(WM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="LAI") {
    if(is.null(ylab)) ylab = expression(paste("Leaf Area Index   ",(m^{2}%.%m^{-2})))
    df = Stand[,c("LAI", "LAIexpanded", "LAIdead")]
    names(df)<-c("Total (live+dead)", "Live unfolded","Dead standing")
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    return(.multiple_dynamics(as.matrix(df), ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantExtraction") {
    extrBal = Soil[,paste("PlantExt",1:nlayers,sep=".")]
    if(!is.null(dates)) extrBal = extrBal[row.names(extrBal) %in% as.character(dates),]
    if(is.null(ylab)) ylab = "Extraction from soil layer (mm)"    
    g<-.multiple_dynamics(as.matrix(extrBal),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers))
    g<-g+geom_abline(slope=0, intercept=0, col="gray")
    return(g)
  } 
  else if(type=="HydraulicRedistribution") {
    hydrIn = Soil[,paste("HydraulicInput",1:nlayers,sep=".")]
    if(!is.null(dates)) hydrIn = hydrIn[row.names(hydrIn) %in% as.character(dates),]
    if(is.null(ylab)) ylab = "Hydraulic input (mm)"    
    return(.multiple_dynamics(as.matrix(hydrIn),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="PlantLAI") {
    OM = PlantsLAI
    if(bySpecies) OM = .averageBySpecies(OM, spnames)
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = expression(paste("Leaf Area Index   ",(m^{2}%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="SoilPlantConductance") {
    OM = Plants$dEdP[,cohorts,drop=FALSE]
    if(bySpecies) OM = .averageByLAISpecies(OM, PlantsLAI, spnames)
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = expression(paste("Average soil-plant conductance ",(mmol%.%m^{-2}%.%s^{-1})))
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("StemPLC", "StemRWC", "LeafRWC", "StemSympRWC", "LeafSympRWC")) {
    OM = Plants[[type]][,cohorts,drop=FALSE]*100
    if(bySpecies) OM = .averageByLAISpecies(OM, PlantsLAI, spnames)
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("PlantPsi", "LeafPsiMin","LeafPsiMax", "StemPsi", "RootPsi", "PlantStress")) {
    OM = Plants[[type]][,cohorts,drop=FALSE]
    if(bySpecies) OM = .averageByLAISpecies(OM, PlantsLAI, spnames)
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type == "LeafPsiRange") {
    OM1 = Plants$LeafPsiMax[,cohorts,drop=FALSE]
    OM2 = Plants$LeafPsiMin[,cohorts,drop=FALSE]
    if(bySpecies) {
      OM1 = .averageByLAISpecies(OM1, PlantsLAI, spnames)
      OM2 = .averageByLAISpecies(OM2, PlantsLAI, spnames)
    } 
    if(!is.null(dates)) {
      OM1 = OM1[row.names(OM1) %in% as.character(dates),]
      OM2 = OM2[row.names(OM2) %in% as.character(dates),]
    }
    return(.multiple_dynamics_range(as.matrix(OM1), as.matrix(OM2),  xlab = xlab, ylab = ylab, ylim = ylim))
  }
  else if(type %in% c("LeafPsiMin_SL", "LeafPsiMax_SL", "GW_SL")) {
    subType = strsplit(type,"_")[[1]][1]
    OM = SunlitLeaves[[subType]][,cohorts,drop=FALSE]
    if(bySpecies) OM = .averageByLAISpecies(OM, PlantsLAI, spnames)
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("LeafPsiMin_SH", "LeafPsiMax_SH", "GW_SH")) {
    subType = strsplit(type,"_")[[1]][1]
    OM = ShadeLeaves[[subType]][,cohorts,drop=FALSE]
    if(bySpecies) OM = .averageByLAISpecies(OM, PlantsLAI, spnames)
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("PlantTranspiration","PlantNetPhotosynthesis", "PlantGrossPhotosynthesis","PlantPhotosynthesis",
                      "PlantAbsorbedSWR","PlantAbsorbedLWR")) {
    subtype = substr(type,6,nchar(type))
    OM = Plants[[subtype]][,cohorts,drop=FALSE]
    if(bySpecies) OM = .averageBySpecies(OM, spnames)
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantWaterBalance") {
    OM = Plants$PlantWaterBalance[,cohorts,drop=FALSE]
    if(bySpecies) OM = .averageBySpecies(OM, spnames)
    if(is.null(ylab)) ylab = expression(paste("Plant water balance   ",(L%.%m^{-2})))
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("TranspirationPerLeaf","PhotosynthesisPerLeaf", "GrossPhotosynthesisPerLeaf", "NetPhotosynthesisPerLeaf",
                      "AbsorbedSWRPerLeaf", "AbsorbedLWRPerLeaf")) {
    subtype = substr(type, 1, nchar(type)-7)
    df = Plants[[subtype]][,cohorts,drop=FALSE]
    df = df/PlantsLAI
    df[PlantsLAI==0] = NA
    if(bySpecies) df = .averageByLAISpecies(df, PlantsLAI, spnames)
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantWUE") {
    if("Photosynthesis" %in% names(Plants)) OM = Plants$Photosynthesis/Plants$Transpiration
    else OM = Plants$GrossPhotosynthesis/Plants$Transpiration
    OM[Plants$Transpiration==0] = 0
    OM = OM[,cohorts,drop=FALSE]
    if(bySpecies) OM = .averageBySpecies(OM, spnames)
    if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
    if(is.null(ylab)) ylab = expression(paste("Plant daily WUE   ",(g*C%.%L^{-1})))
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="Temperature") {
    if(is.null(ylab)) ylab = "Temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Above-canopy"]] = x$Temperature$Tatm_mean
    df[["Inside-canopy"]] = x$Temperature$Tcan_mean
    df[["Soil"]] = x$Temperature$Tsoil_mean
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="AirTemperature") {
    if(is.null(ylab)) ylab = "Above-canopy temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tatm_mean
    df[["Minimum"]] = x$Temperature$Tatm_min
    df[["Maximum"]] = x$Temperature$Tatm_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="CanopyTemperature") {
    if(is.null(ylab)) ylab = "Canopy temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tcan_mean
    df[["Minimum"]] = x$Temperature$Tcan_min
    df[["Maximum"]] = x$Temperature$Tcan_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="SoilTemperature") {
    if(is.null(ylab)) ylab = "Soil temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tsoil_mean
    df[["Minimum"]] = x$Temperature$Tsoil_min
    df[["Maximum"]] = x$Temperature$Tsoil_max
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="CanopyEnergyBalance") {
    if(is.null(ylab)) ylab = expression(MJ%.%m^{-2})    
    df = data.frame(row.names=row.names(x$EnergyBalance))
    df[["Balance"]] = x$EnergyBalance$Ebalcan
    df[["SWR abs. from atm."]] = x$EnergyBalance$SWRcanin 
    df[["LWR abs. from atm."]] = x$EnergyBalance$LWRcanin
    df[["LWR abs. from soil"]] = x$EnergyBalance$LWRsoilcan
    df[["LWR emmited"]] = -x$EnergyBalance$LWRcanout
    df[["Latent heat"]] = -x$EnergyBalance$LEcan
    df[["Convection can./atm."]] = -x$EnergyBalance$Hcan
    df[["Convection soil/can."]] = -x$EnergyBalance$Hcansoil
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    g<-.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim)
    return(g)
    # lines(dates, x$EnergyBalance$Ebalcan, col="black",...)
    # lines(dates, x$EnergyBalance$SWRcanin, col="red",...)
    # lines(dates, x$EnergyBalance$LWRcanin, col="brown",...)
    # lines(dates, -x$EnergyBalance$LWRcanout, col="blue",...)
    # lines(dates, x$EnergyBalance$LWRsoilcan, col="orange",...)
    # lines(dates, -x$EnergyBalance$LEcan, col="green",...)
    # lines(dates, -x$EnergyBalance$Hcan, col="gray",...)
    # lines(dates, -x$EnergyBalance$Hcansoil, col="dark gray",...)
    # legend("topright", bty="n", col=c("red","brown","orange", "blue","green", "gray", "dark gray", "black"), lty=1,
    #        legend=c("SWR abs. from atm.","LWR abs. from atm.","LWR abs. from soil","LWR emmited", "Latent heat (L)",
    #                 "Convection can./atm.","Convection soil/can.", "Balance"),...)        
  } 
  else if(type=="SoilEnergyBalance") {
    if(is.null(ylab)) ylab = expression(MJ%.%m^{-2})    
    df = data.frame(row.names=row.names(x$EnergyBalance))
    df[["Balance"]] = x$EnergyBalance$Ebalsoil
    df[["SWR abs. from atm."]] = x$EnergyBalance$SWRsoilin
    df[["LWR abs. from atm."]] = x$EnergyBalance$LWRsoilin
    df[["LWR abs. from canopy"]] = x$EnergyBalance$LWRcanout
    df[["LWR emmited"]] = -x$EnergyBalance$LWRsoilout
    df[["Convection soil/can."]] = x$EnergyBalance$Hcansoil
    df[["Latent heat"]] = -x$EnergyBalance$LEsoil
    if(!is.null(dates)) df = df[row.names(df) %in% as.character(dates),]
    g<-.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim)
    return(g)
    # lines(dates, x$EnergyBalance$Ebalsoil, col="black",...)
    # lines(dates, x$EnergyBalance$SWRsoilin, col="red",...)
    # lines(dates, x$EnergyBalance$LWRsoilin, col="brown",...)
    # lines(dates, x$EnergyBalance$LWRcanout, col="orange",...)
    # lines(dates, -x$EnergyBalance$LEsoil, col="green",...)
    # lines(dates, -x$EnergyBalance$LWRsoilout, col="blue",...)
    # lines(dates, x$EnergyBalance$Hcansoil, col="gray",...)
    # legend("topright", bty="n", col=c("red","brown","orange", "blue", "green", "gray", "black"), lty=1,
    #        legend=c("SWR abs. from atm.","LWR abs. from atm.", "LWR abs. from canopy","LWR emmited","Latent heat (L)",  "Convection soil/can.", "Balance"),...)        
  }
}
plot.growth<-function(x, type="PET_Precipitation", cohorts = NULL, bySpecies = FALSE, 
                      dates = NULL, subdaily = FALSE, 
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...) {
  
  if(subdaily) return(.plotsubdaily(x,type, cohorts, bySpecies, dates, 
                                    xlim, ylim, xlab, ylab))
  
  # Get common elements
  input = x$growthInput
  soilInput = x$soilInput
  PlantGrowth = x$PlantGrowth
  Plants = x$Plants
  PCB = x$PlantCarbonBalance
  
  nlayers = length(soilInput$W)
  
  transpMode = input$control$transpirationMode
  
  TYPES_GROWTH = .getDailyGROWTHPlotTypes(transpMode)
  TYPES_SWB = .getDailySPWBPlotTypes(transpMode)  

  type = match.arg(type,TYPES_GROWTH)  
  
  if(is.null(cohorts))  cohorts = row.names(input$cohorts)
  spnames = as.character(input$cohorts[cohorts,"Name"])
  PlantsLAI = Plants$LAI[,cohorts, drop=FALSE]
  
  if(type %in% TYPES_SWB) {
    plot.spwb(x,type, cohorts, bySpecies, dates, subdaily, xlim, ylim, xlab, ylab, ...)
  } 
  else if(type %in% c("GrossPhotosynthesis", "MaintenanceRespiration",  "GrowthRespiration", "CarbonBalance", 
                      "SugarLeaf","StarchLeaf","SugarSapwood","StarchSapwood", "SugarTransport", "LeafPI0", "StemPI0")) {
      OM = PCB[[type]][,cohorts,drop=FALSE]
      if(bySpecies) OM = .averageBySpecies(OM, spnames)
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab=.getYLab(type)
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type %in% c("SapwoodArea", "LeafArea","SAgrowth", "LAgrowth", "HuberValue")) {
      OM = PlantGrowth[[type]][,cohorts,drop=FALSE]
      if(bySpecies) OM = .averageBySpecies(OM, spnames)
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = .getYLab(type)
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
}

.averageSubdailyBySpecies<-function(OM, spnames) {
  if(ncol(OM)>2) OM = cbind(OM[,1],t(apply(OM[,-1],1, tapply, spnames, sum, na.rm=T)))
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
  else if(type=="LeafAbsorbedLWR") {
    mSu = extractSubdaily(x, "SunlitLeaves$Abs_LWR", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$Abs_LWR", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Absorbed LWR per leaf area ",(W%.%m^{-2})))
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
    mSu = extractSubdaily(x, "SunlitLeaves$GW", dates)[,c("datetime", cohorts), drop=FALSE]
    mSh = extractSubdaily(x, "ShadeLeaves$GW", dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=expression(paste("Stomatal conductance ", (mmol%.%m^{-2}%.%s^{-1})))
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
  else if(type %in% c("GrossPhotosynthesis", "MaintenanceRespiration", "GrowthRespiration", "CarbonBalance",
                      "SugarLeaf","StarchLeaf","SugarSapwood","StarchSapwood", "SugarTransport")) {
    m = extractSubdaily(x, type, dates)[,c("datetime", cohorts), drop=FALSE]
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_dynamics_subdaily(m,  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
}