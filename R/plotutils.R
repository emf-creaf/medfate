.getDailySPWBPlotTypes<-function(transpirationMode = "Granier") {
   if(transpirationMode=="Granier") {
     TYPES = c("PET_Precipitation","PET_NetRain","Snow","Evapotranspiration",
               "SoilPsi","SoilTheta","SoilRWC","SoilVol", 
               "Export", "LAI", "WTD",
               "PlantExtraction","PlantLAI",
               "PlantStress", "PlantPsi","PlantPhotosynthesis", "PlantTranspiration", "PlantWUE",
               "PhotosynthesisPerLeaf","TranspirationPerLeaf")
   } else {
     TYPES = c("PET_Precipitation","PET_NetRain","Snow","Evapotranspiration",
               "SoilPsi","SoilTheta", "SoilRWC", "SoilVol", 
               "Export", "LAI", "WTD",
               "PlantExtraction","HydraulicRedistribution",
               "PlantLAI",
               "SoilPlantConductance","PlantStress", 
               "PlantNetPhotosynthesis", "PlantGrossPhotosynthesis", "PlantTranspiration","PlantWUE",
               "NetPhotosynthesisPerLeaf","GrossPhotosynthesisPerLeaf","TranspirationPerLeaf", 
               "LeafPsiMin", "LeafPsiMax", 
               "GW_SL", "GW_SH","LeafPsiMin_SL", "LeafPsiMax_SL", "LeafPsiMin_SH", "LeafPsiMax_SH",
               "StemPsi","RootPsi","StemPLC", "StemRWC", "LeafRWC", 
               "PlantWaterBalance",
               "PlantAbsorbedSWR", "AbsorbedSWRPerLeaf",
               "PlantAbsorbedLWR", "AbsorbedLWRPerLeaf",
               "Temperature","AirTemperature","SoilTemperature", "CanopyTemperature",
               "CanopyEnergyBalance", "SoilEnergyBalance")
   }
  return(TYPES)
}

.getSubdailySPWBPlotTypes<-function(){
  TYPES = c("LeafPsi","LeafPsiAverage","RootPsi", "StemPsi", 
            "StemPLC","StemRWC", "LeafRWC",
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
            "SugarTransport",
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
  if(type=="GrossPhotosynthesis") ylab=expression(paste("Gross photosynthesis  ", (gGluc%.%gdry^{-1})))
  else if(type=="MaintenanceRespiration") ylab=expression(paste("Maintenance respiration  ", (gGluc%.%gdry^{-1})))
  else if(type=="GrowthRespiration") ylab=expression(paste("Growth respiration  ", (gGluc%.%gdry^{-1})))
  else if(type=="CarbonBalance") ylab=expression(paste("Carbon balance  ", (gGluc%.%gdry^{-1})))
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