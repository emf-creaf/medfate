plot.spwb<-function(x, type="PET_Precipitation", bySpecies = FALSE,
                    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...) {
  dates = as.Date(rownames(x$WaterBalance))
  
  input = x$spwbInput
  soilInput = x$soilInput
  
  WaterBalance = x$WaterBalance
  Soil = x$Soil
  Stand = x$Stand
  nlayers = length(soilInput$W)
  
  transpMode = input$control$transpirationMode
  
  TYPES = c("PET_Precipitation","PET_NetRain","Snow","Evapotranspiration",
            "SoilPsi","SoilTheta","SoilRWC","SoilVol", 
            "Export", "LAI", "WTD",
            "PlantExtraction","PlantLAI",
            "PlantStress", "PlantPsi","PlantPhotosynthesis", "PlantTranspiration", "PlantWUE",
            "PlantPhotosynthesisPerLeaf","PlantTranspirationPerLeaf")
  if(transpMode=="Sperry") {
    TYPES = c("PET_Precipitation","PET_NetRain","Snow","Evapotranspiration",
              "SoilPsi","SoilTheta", "SoilRWC", "SoilVol", 
              "Export", "LAI", "WTD",
              "PlantExtraction","HydraulicRedistribution",
              "PlantLAI",
              "SoilPlantConductance","PlantStress", 
              "PlantPhotosynthesis", "PlantTranspiration","PlantWUE",
              "PlantPhotosynthesisPerLeaf","PlantTranspirationPerLeaf", 
              "LeafPsiMin", "LeafPsiMax", "LeafPsiMin_SL", "LeafPsiMax_SL", "LeafPsiMin_SH", "LeafPsiMax_SH",
              "StemPsi","RootPsi","StemPLC", "StemRWC", "LeafRWC", 
              "PlantWaterBalance",
              "PlantAbsorbedSWR", "PlantAbsorbedSWRPerLeaf",
              "PlantAbsorbedLWR", "PlantAbsorbedLWRPerLeaf",
              "AirTemperature","SoilTemperature", "CanopyTemperature",
              "CanopyEnergyBalance", "SoilEnergyBalance")
  } 
  type = match.arg(type,TYPES)  
  if(is.null(xlab)) xlab = ""
  if(type=="PET_Precipitation") {
    if(is.null(ylab)) ylab = expression(L%.%m^{-2})
    if(is.null(ylab)) ylab = expression(L%.%m^{-2}) 
    df = data.frame(row.names=row.names(x$WaterBalance))
    df[["PET"]] = x$WaterBalance$PET
    df[["Precipitation"]] = x$WaterBalance$Precipitation
    df[["Snow"]] = x$WaterBalance$Snow
    df[["Date"]] = as.Date(row.names(x$WaterBalance))
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
    df = data.frame(row.names=row.names(x$WaterBalance))
    df[["PET"]] = x$WaterBalance$PET
    df[["NetRain"]] = x$WaterBalance$NetRain
    df[["Date"]] = as.Date(row.names(x$WaterBalance))
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
    df = data.frame(row.names=row.names(x$WaterBalance))
    df[["Snow"]] = x$WaterBalance$Snow
    df[["Snowpack"]] = x$Soil$SWE
    df[["Date"]] = as.Date(row.names(x$WaterBalance))
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
    df = data.frame(row.names=row.names(x$WaterBalance))
    df[["Total evapotranspiration"]] = x$WaterBalance$Evapotranspiration
    df[["Interception evaporation"]] = x$WaterBalance$Interception
    df[["Plant transpiration"]] = x$WaterBalance$Transpiration
    df[["Bare soil evaporation"]] = x$WaterBalance$SoilEvaporation
    return(.multiple_dynamics(as.matrix(df), ylab=ylab, xlab=xlab, ylim = ylim))
  } 
  else if(type=="LAI") {
    if(is.null(ylab)) ylab = expression(paste("Leaf Area Index   ",(m^{2}%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(Stand[,c("LAIcell", "LAIcelldead")]), ylab = ylab, ylim = ylim))
  } 
  else if(type=="WTD") {
    if(is.null(ylab)) ylab = expression(paste("Water table depth  (mm)"))
    xv = Soil$WTD
    names(xv) = row.names(Soil)
    return(.single_dynamics(xv, ylab = ylab, ylim = ylim))
  } 
  else if(type=="Export") {
    if(is.null(ylab)) ylab =  expression(L%.%m^{-2})    
    df = data.frame(row.names=row.names(x$WaterBalance))
    df[["Export"]] = x$WaterBalance$DeepDrainage + x$WaterBalance$Runoff
    df[["DeepDrainage"]] = x$WaterBalance$DeepDrainage
    df[["Runoff"]] = x$WaterBalance$Runoff 
    df[["Date"]]= as.Date(row.names(x$WaterBalance))
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
    return(.multiple_dynamics(as.matrix(MLM), ylab = ylab, ylim = ylim,
                              xlab=xlab, labels = c("Total", paste("Layer", 1:nlayers))))
  } 
  else {
    plot.pwb(x, type=type, bySpecies = bySpecies,
             xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  }
}

plot.pwb<-function(x, type="PlantTranspiration", bySpecies = FALSE,
                   xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...) {
  dates = as.Date(rownames(x$WaterBalance))
  
  input = x$spwbInput
  soilInput = x$soilInput
  
  WaterBalance = x$WaterBalance
  Soil = x$Soil
  Stand = x$Stand
  nlayers = length(soilInput$W)
  
  transpMode = input$control$transpirationMode
  
  TYPES = c("SoilPsi","SoilTheta", "SoilRWC",
            "PlantExtraction","PlantLAI",
            "PlantStress", "PlantPsi","PlantPhotosynthesis", "PlantTranspiration", "PlantWUE",
            "PlantPhotosynthesisPerLeaf","PlantTranspirationPerLeaf")
  if(transpMode=="Sperry") {
    TYPES = c("SoilPsi","SoilTheta", "SoilRWC",
              "PlantExtraction","HydraulicRedistribution",
              "PlantLAI",
              "SoilPlantConductance","PlantStress", 
              "PlantPhotosynthesis", "PlantTranspiration","PlantWUE",
              "PlantPhotosynthesisPerLeaf","PlantTranspirationPerLeaf", 
              "LeafPsiMin", "LeafPsiMax", "LeafPsiMin_SL", "LeafPsiMax_SL", "LeafPsiMin_SH", "LeafPsiMax_SH",
              "StemPsi","RootPsi","StemPLC", "StemRWC", "LeafRWC", 
              "PlantWaterBalance",
              "PlantAbsorbedSWR", "PlantAbsorbedSWRPerLeaf",
              "PlantAbsorbedLWR", "PlantAbsorbedLWRPerLeaf",
              "AirTemperature","SoilTemperature", "CanopyTemperature",
              "CanopyEnergyBalance", "SoilEnergyBalance")
  } 
  type = match.arg(type,TYPES)  
  if(is.null(xlab)) xlab = ""  
  
  if(type=="SoilPsi") {
    PsiM = Soil[,paste("psi",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Soil water potential (MPa)"    
    return(.multiple_dynamics(as.matrix(PsiM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="SoilTheta") {
    WM = Soil[,paste("W",1:nlayers,sep=".")]
    theta_FC = soil_thetaFC(soilInput, model = input$control$soilFunctions)
    WM = 100*sweep(WM, 2,theta_FC, "*")
    if(is.null(ylab)) ylab = "Soil moisture (% volume)"
    return(.multiple_dynamics(as.matrix(WM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="SoilRWC") {
    WM = Soil[,paste("W",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Soil moisture (% field capacity)"
    return(.multiple_dynamics(as.matrix(WM),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="PlantExtraction") {
    extrBal = Soil[,paste("PlantExt",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Extraction from soil layer (mm)"    
    g<-.multiple_dynamics(as.matrix(extrBal),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers))
    g<-g+geom_abline(slope=0, intercept=0, col="gray")
    return(g)
  } 
  else if(type=="HydraulicRedistribution") {
    hydrIn = Soil[,paste("HydraulicInput",1:nlayers,sep=".")]
    if(is.null(ylab)) ylab = "Hydraulic input (mm)"    
    return(.multiple_dynamics(as.matrix(hydrIn),  xlab = xlab, ylab = ylab, ylim = ylim,
                              labels = paste("Layer", 1:nlayers)))
  } 
  else if(type=="PlantLAI") {
    OM = x$PlantLAI
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
    if(is.null(ylab)) ylab = expression(paste("Leaf Area Index   ",(m^{2}%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="SoilPlantConductance") {
    OM = x$dEdP
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = expression(paste("Average soil-plant conductance ",(mmol%.%m^{-2}%.%s^{-1})))
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantStress") {
    OM = x$PlantStress
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Drought stress [0-1]"
    if(is.null(ylim)) ylim = c(0,1)
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="StemPLC") {
    OM = x$StemPLC*100
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Percent loss conductance in stem [%]"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="StemRWC") {
    OM = x$StemRWC*100
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Relative water content in stem [%]"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafRWC") {
    OM = x$LeafRWC*100
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Relative water content in leaf [%]"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantPsi") {
    OM = x$PlantPsi
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Plant water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsiMin") {
    OM = x$LeafPsiMin
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Minimum (midday) leaf water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsiMax") {
    OM = x$LeafPsiMax
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Maximum (predawn) leaf water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsiMin_SL") {
    OM = x$LeafPsiMin_SL
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Minimum (midday) sunlit leaf water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsiMax_SL") {
    OM = x$LeafPsiMax_SL
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Maximum (predawn) sunlit leaf water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsiMin_SH") {
    OM = x$LeafPsiMin_SH
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Minimum (midday) shade leaf water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="LeafPsiMax_SH") {
    OM = x$LeafPsiMax_SH
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Maximum (predawn) shade leaf water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  }
  else if(type=="StemPsi") {
    OM = x$StemPsi
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Midday stem water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="RootPsi") {
    OM = x$RootPsi
    if(bySpecies) {
      lai1 = t(apply(x$PlantLAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(x$PlantLAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
    if(is.null(ylab)) ylab = "Midday root crown water potential (MPa)"
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantTranspiration") {
    OM = x$PlantTranspiration
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant transpiration   ",(L%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantTranspirationPerLeaf") {
    df = x$PlantTranspiration
    if(bySpecies) {
      m1 = apply(df,1, tapply, input$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$PlantLAI,1,tapply, input$cohorts$Name, sum, na.rm=T)
      df = t(m1/lai1)
    } else {
      df = df/x$PlantLAI
      df[x$PlantLAI==0] = NA
    }
    if(is.null(ylab)) ylab = expression(paste("Plant transpiration per leaf area  ",(L%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantWaterBalance") {
    OM = x$PlantWaterBalance
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant water balance   ",(L%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantPhotosynthesis") {
    df = x$PlantPhotosynthesis
    if(bySpecies) {
      df = t(apply(df,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant photosynthesis   ",(g*C%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantPhotosynthesisPerLeaf") {
    df = x$PlantPhotosynthesis
    if(bySpecies) {
      m1 = apply(df,1, tapply, input$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$PlantLAI,1,tapply, input$cohorts$Name, sum, na.rm=T)
      df = t(m1/lai1)
    } else {
      df = df/x$PlantLAI
      df[x$PlantLAI==0] = NA
    }
    if(is.null(ylab)) ylab = expression(paste("Plant photosynthesis per leaf area   ",(g*C%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantWUE") {
    OM = x$PlantPhotosynthesis/x$PlantTranspiration
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant daily WUE   ",(g*C%.%L^{-1})))
    return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantAbsorbedSWR") {
    df = x$PlantAbsorbedSWR
    if(bySpecies) {
      df = t(apply(df,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant absorbed SWR  ",(MJ%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantAbsorbedSWRPerLeaf") {
    df = x$PlantAbsorbedSWR
    if(bySpecies) {
      m1 = apply(df,1, tapply, input$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$PlantLAI,1,tapply, input$cohorts$Name, sum, na.rm=T)
      df = t(m1/lai1)
    } else {
      df = df/x$PlantLAI
      df[x$PlantLAI==0] = NA
    }
    if(is.null(ylab)) ylab = expression(paste("Plant absorbed SWR per leaf area  ",(MJ%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantAbsorbedLWR") {
    df = x$PlantAbsorbedLWR
    if(bySpecies) {
      df = t(apply(df,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
    if(is.null(ylab)) ylab = expression(paste("Plant absorbed LWR  ",(MJ%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="PlantAbsorbedLWRPerLeaf") {
    df = x$PlantAbsorbedLWR
    if(bySpecies) {
      m1 = apply(df,1, tapply, input$cohorts$Name, sum, na.rm=T)
      lai1 = apply(x$PlantLAI,1,tapply, input$cohorts$Name, sum, na.rm=T)
      df = t(m1/lai1)
    } else {
      df = df/x$PlantLAI
      df[x$PlantLAI==0] = NA
    }
    if(is.null(ylab)) ylab = expression(paste("Plant absorbed LWR per leaf area  ",(MJ%.%m^{-2})))
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab = ylab, ylim = ylim))
  } 
  else if(type=="AirTemperature") {
    if(is.null(ylab)) ylab = "Above-canopy temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tatm_mean
    df[["Minimum"]] = x$Temperature$Tatm_min
    df[["Maximum"]] = x$Temperature$Tatm_max
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="CanopyTemperature") {
    if(is.null(ylab)) ylab = "Canopy temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tcan_mean
    df[["Minimum"]] = x$Temperature$Tcan_min
    df[["Maximum"]] = x$Temperature$Tcan_max
    return(.multiple_dynamics(as.matrix(df),  xlab = xlab, ylab=ylab, ylim = ylim))
  } 
  else if(type=="SoilTemperature") {
    if(is.null(ylab)) ylab = "Soil temperature (Celsius)"
    df = data.frame(row.names=row.names(x$Temperature))
    df[["Mean"]] = x$Temperature$Tsoil_mean
    df[["Minimum"]] = x$Temperature$Tsoil_min
    df[["Maximum"]] = x$Temperature$Tsoil_max
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
