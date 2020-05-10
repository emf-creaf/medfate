.getDailyPWBPlotTypes<-function(transpirationMode = "Granier") {
  if(transpirationMode=="Granier") {
    TYPES = c("SoilPsi","SoilTheta", "SoilRWC", "LAI",
              "PlantExtraction","PlantLAI", 
              "PlantStress", "PlantPsi","StemPLC",
              "PlantPhotosynthesis", "PlantTranspiration", "PlantWUE",
              "PhotosynthesisPerLeaf","TranspirationPerLeaf")
  } else {
    TYPES = c("SoilPsi","SoilTheta", "SoilRWC", "LAI",
              "PlantExtraction","HydraulicRedistribution",
              "PlantLAI",
              "SoilPlantConductance","PlantStress", 
              "PlantNetPhotosynthesis", "PlantGrossPhotosynthesis", "PlantTranspiration","PlantWUE",
              "NetPhotosynthesisPerLeaf","GrossPhotosynthesisPerLeaf","TranspirationPerLeaf", 
              "GW_SL", "GW_SH", "LeafPsiRange",
              "LeafPsiMin", "LeafPsiMax", "LeafPsiMin_SL", "LeafPsiMax_SL", "LeafPsiMin_SH", "LeafPsiMax_SH",
              "StemPsi","RootPsi","StemPLC", "StemRWC", "LeafRWC", "StemSympRWC", "LeafSympRWC", 
              "PlantWaterBalance",
              "PlantAbsorbedSWR", "AbsorbedSWRPerLeaf",
              "PlantAbsorbedLWR", "AbsorbedLWRPerLeaf",
              "Temperature","AirTemperature","SoilTemperature", "CanopyTemperature",
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
            "LeafAbsorbedSWR","LeafAbsorbedLWR",
            "LeafCi", "LeafIntrinsicWUE",
            "LeafVPD","LeafStomatalConductance", "LeafTemperature",
            "Temperature","CanopyEnergyBalance", "SoilEnergyBalance", 
            "PlantWaterBalance", "WaterBalancePerLeaf")
  return(TYPES)
}
.getDailyGROWTHPlotTypes<-function(transpirationMode = "Granier"){
  TYPES = c("GrossPhotosynthesis","MaintenanceRespiration","GrowthRespiration", "CarbonBalance",
            "SugarTransport", "LeafPI0", "StemPI0",
            "SugarLeaf", "SugarSapwood", "StarchLeaf", "StarchSapwood","SugarTransport",
            "SapwoodArea", "LeafArea","SAgrowth", "LAgrowth", "HuberValue",
            .getDailySPWBPlotTypes(transpirationMode))
  return(TYPES)
}
.getSubdailyGROWTHPlotTypes<-function(){
  TYPES = c("GrossPhotosynthesis","MaintenanceRespiration","GrowthRespiration", "CarbonBalance",
            "SugarTransport",
            "SugarLeaf", "SugarSapwood", "StarchLeaf", "StarchSapwood","SugarTransport",
            .getSubdailySPWBPlotTypes())
  return(TYPES)
}
.getYLab<-function(type) {
  ylab="Unknown"
  if(type=="PlantTranspiration") ylab = expression(paste("Plant transpiration ",(L%.%m^{-2})))
  else if(type=="PlantPhotosynthesis") ylab = expression(paste("Plant photosynthesis ",(g*C%.%m^{-2})))
  else if(type=="PlantGrossPhotosynthesis") ylab = expression(paste("Plant gross photosynthesis ",(g*C%.%m^{-2})))
  else if(type=="PlantNetPhotosynthesis") ylab = expression(paste("Plant net photosynthesis ",(g*C%.%m^{-2})))
  else if(type=="PlantAbsorbedSWR") ylab = expression(paste("Plant absorbed SWR ",(MJ%.%m^{-2})))
  else if(type=="PlantAbsorbedLWR") ylab = expression(paste("Plant absorbed LWR ",(MJ%.%m^{-2})))
  else if(type=="TranspirationPerLeaf") ylab = expression(paste("Transpiration per leaf area ",(L%.%m^{-2})))
  else if(type=="PhotosynthesisPerLeaf") ylab = expression(paste("Photosynthesis per leaf area ",(g*C%.%m^{-2})))
  else if(type=="GrossPhotosynthesisPerLeaf") ylab = expression(paste("Gross photosynthesis per leaf area ",(g*C%.%m^{-2})))
  else if(type=="NetPhotosynthesisPerLeaf") ylab = expression(paste("Net photosynthesis per leaf area ",(g*C%.%m^{-2})))
  else if(type=="AbsorbedSWRPerLeaf") ylab = expression(paste("Absorbed SWR per leaf area ",(MJ%.%m^{-2})))
  else if(type=="AbsorbedLWRPerLeaf") ylab = expression(paste("Absorbed LWR per leaf area ",(MJ%.%m^{-2})))
  else if(type=="GrossPhotosynthesis") ylab=expression(paste("Gross photosynthesis ", (gGluc%.%gdry^{-1})))
  else if(type=="MaintenanceRespiration") ylab=expression(paste("Maintenance respiration ", (gGluc%.%gdry^{-1})))
  else if(type=="GrowthRespiration") ylab=expression(paste("Growth respiration ", (gGluc%.%gdry^{-1})))
  else if(type=="CarbonBalance") ylab=expression(paste("Carbon balance ", (gGluc%.%gdry^{-1})))
  else if(type=="SugarLeaf") ylab=expression(paste("Leaf sugar concentration  ", (mol%.%L^{-1})))
  else if(type=="StarchLeaf") ylab=expression(paste("Leaf starch concentration  ", (mol%.%L^{-1})))
  else if(type=="SugarSapwood") ylab=expression(paste("Sapwood sugar concentration  ", (mol%.%L^{-1})))
  else if(type=="StarchSapwood") ylab=expression(paste("Sapwood starch concentration  ", (mol%.%L^{-1})))
  else if(type=="SugarTransport") ylab=expression(paste("Floem sugar transport rate ", (mmol%.%s^{-1})))
  else if(type=="SapwoodArea")  ylab = expression(paste("Sapwood area  ",(cm^2)))
  else if(type=="LeafArea")  ylab = expression(paste("Leaf area  ",(m^2)))
  else if(type=="HuberValue")  ylab = expression(paste("Huber value  ",(cm^2 %.% m^{-2})))
  else if(type=="SAgrowth") ylab = expression(paste("Sapwood area growth rate ",(cm^2 %.% cm^{-2} %.% d^{-1})))
  else if(type=="LAgrowth") ylab = expression(paste("Leaf area growth rate ",(m^2 %.% cm^{-2} %.% d^{-1})))
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
  else if(type=="LeafPsiMax_SL") ylab = "Maximum (midday) sunlit leaf water potential (MPa)"
  else if(type=="LeafPsiMin_SH") ylab = "Minimum (midday) shade leaf water potential (MPa)"
  else if(type=="LeafPsiMax_SH") ylab = "Maximum (midday) shade leaf water potential (MPa)"
  else if(type=="GW_SH") ylab = expression(paste("Shade leaf stomatal conductance ",(mmol%.%m^{-2}%.%s^{-1})))
  else if(type=="GW_SL") ylab = expression(paste("Sunlit leaf stomatal conductance ",(mmol%.%m^{-2}%.%s^{-1})))
  return(ylab)
}

.multiple_x<-function(x, y, xlab = "", ylab=NULL, xlim = NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x)
  df = data.frame("X" = as.vector(x), 
                  "Y" = y,
                  "Cohort" = gl(length(colnames(x)), nrow(x), labels=labels))
  g<-ggplot(df, aes_string(x="X", y="Y"))+
    geom_path(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(xlim)) g <- g+xlim(xlim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}
.multiple_y<-function(x, y, xlab = "", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(y)
  df = data.frame("Y" = as.vector(y), "X" = x,
                  "Cohort" = gl(length(colnames(y)), nrow(y), labels=labels))
  g<-ggplot(df, aes_string(x="X", y="Y"))+
    geom_path(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.multiple_dynamics<-function(x, xlab = "", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x)
  df = data.frame("Y" = as.vector(x), 
                  "Date" = as.Date(rownames(x)),
                  "Cohort" = gl(length(colnames(x)), nrow(x), labels=labels))
  g<-ggplot(df, aes_string(x="Date", y="Y"))+
    geom_line(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}
.multiple_dynamics_range<-function(x1,x2, xlab = "", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x1)
  df = data.frame("Y1" = as.vector(x1),
                  "Y2" = as.vector(x2),
                  "Date" = as.Date(rownames(x1)),
                  "Cohort" = gl(length(colnames(x1)), nrow(x1), labels=labels))
  g<-ggplot(df, aes_string(x="Date"))+
    geom_ribbon(aes_string(ymin = "Y2", ymax="Y1",fill="Cohort"), alpha = 0.5)+
    scale_fill_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}


.multiple_dynamics_subdaily<-function(x, xlab = "", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x)[-1]
  df = data.frame("Y" = as.vector(as.matrix(x[,-1])), 
                  "DateTime" = as.POSIXct(x$datetime),
                  "Cohort" = gl(length(names(x)[-1]), nrow(x), labels=labels))
  g<-ggplot(df, aes_string(x="DateTime", y="Y"))+
    geom_line(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}
.multiple_dynamics_subdaily_sunlit_shade<-function(x_sl, x_sh, xlab = "", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x_sl)[-1]
  df_sl = data.frame("Y" = as.vector(as.matrix(x_sl[,-1])), 
                     "DateTime" = as.POSIXct(x_sl$datetime),
                     "Cohort" = gl(length(colnames(x_sl)[-1]), nrow(x_sl), labels=labels),
                     "LeafType" = "Sunlit", stringsAsFactors = F)
  df_sh = data.frame("Y" = as.vector(as.matrix(x_sh[,-1])), 
                     "DateTime" = as.POSIXct(x_sh$datetime),
                     "Cohort" = gl(length(colnames(x_sh)[-1]), nrow(x_sh), labels=labels),
                     "LeafType" = "Shade", stringsAsFactors = F)
  df = as.data.frame(rbind(df_sl, df_sh), stringsAsFactors = F)
  df$LeafType = factor(df$LeafType, levels =c("Sunlit","Shade"))
  g<-ggplot(df, aes_string(x="DateTime", y="Y"))+
    geom_line(aes_string(col="Cohort", linetype = "Cohort"))+
    facet_wrap(~LeafType, ncol=1)+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}
.single_dynamics<-function(x, xlab="", ylab=NULL, ylim = NULL) {
  df = data.frame("Y" = x, 
                  "Date" = as.Date(names(x)))
  g<-ggplot(df, aes_string(x = "Date", y= "Y"))+
    geom_line()+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.single_subday_dynamics <-function(x, xlab="Time step", ylab=NULL, ylim = NULL) {
  df = data.frame("Y" = x, 
                  "TimeStep" = 1:length(x))
  g<-ggplot(df, aes_string(x = "TimeStep", y= "Y"))+
    geom_line()+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.multiple_subday_dynamics<-function(x, xlab = "Time step", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x)
  df = data.frame("Y" = as.vector(x), 
                  "TimeStep" = as.numeric(rownames(x)),
                  "Cohort" = gl(length(colnames(x)), nrow(x), labels=labels))
  g<-ggplot(df, aes_string(x="TimeStep", 
                           y="Y"))+
    geom_line(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.multiple_subday_dynamics_sunlit_shade<-function(x_sl, x_sh, xlab = "Time step", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x_sl)
  df_sl = data.frame("Y" = as.vector(x_sl), 
                     "TimeStep" = as.numeric(rownames(x_sl)),
                     "Cohort" = gl(length(colnames(x_sl)), nrow(x_sl), labels=labels),
                     "LeafType" = "Sunlit", stringsAsFactors = F)
  df_sh = data.frame("Y" = as.vector(x_sh), 
                     "TimeStep" = as.numeric(rownames(x_sh)),
                     "Cohort" = gl(length(colnames(x_sh)), nrow(x_sh), labels=labels),
                     "LeafType" = "Shade", stringsAsFactors = F)
  df = as.data.frame(rbind(df_sl, df_sh), stringsAsFactors = F)
  df$LeafType = factor(df$LeafType, levels =c("Sunlit","Shade"))
  g<-ggplot(df, aes_string(x="TimeStep", y="Y"))+
    geom_line(aes_string(col="Cohort", linetype = "Cohort"))+
    facet_wrap(~LeafType, ncol=2)+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}