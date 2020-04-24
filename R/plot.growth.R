plot.growth<-function(x, type="PET_Precipitation", bySpecies = FALSE, 
                      dates = NULL, subdaily = FALSE, 
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...) {

  if(subdaily) return(.plotsubdaily(x,type, bySpecies, dates, 
                                    xlim, ylim, xlab, ylab))
  
  # Get common elements
  input = x$spwbInput
  soilInput = x$soilInput
  WaterBalance = x$WaterBalance
  Soil = x$Soil
  Stand = x$Stand
  Plants = x$Plants
  
  nlayers = length(soilInput$W)
  
  transpMode = input$control$transpirationMode
  
  TYPES_SWB = .getDailySPWBPlotTypes(transpMode)  
  TYPES = c(TYPES_SWB, "PlantRespiration", "PlantRespirationPerLeaf", "PlantRespirationPerIndividual",
            "PlantCBalance", "PlantCBalancePerLeaf","PlantCBalancePerIndividual",
            "PlantCstorageFast", "PlantCstorageSlow", "PlantSAgrowth", "PlantRelativeSAgrowth", "PlantSA",
            "PlantLAIlive","PlantLAIdead")

  type = match.arg(type,TYPES)  

  if(type %in% TYPES_SWB) {
    x$PlantLAI = x$PlantLAIlive
    plot.spwb(x,type, bySpecies, xlim, ylim, xlab, ylab, ...)
  } else {
    if(is.null(xlab)) xlab = ""
    if(type=="PlantRespiration") {
      OM = x$PlantRespiration
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Plant respiration ",(g*C%.%m^{-2})))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
    else if(type=="PlantRespirationPerLeaf") {
      OM = x$PlantRespiration
      if(bySpecies) {
        m1 = apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T)
        lai1 = apply(x$PlantLAIlive,1,tapply, input$cohorts$Name, sum, na.rm=T)
        OM = t(m1/lai1)
      } else {
        OM = OM/x$PlantLAIlive
        OM[x$PlantLAIlive==0] = NA
      }
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Plant respiration per leaf area ",(g*C%.%m^{-2})))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
    else if(type=="PlantRespirationPerIndividual") {
      OM = x$PlantRespiration
      if(bySpecies) {
        m1 = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
        dens = as.numeric(tapply(input$above$N/10000, input$cohorts$Name, sum, na.rm=T))
        OM = sweep(m1, 2, dens, "/")
      } else {
        OM = sweep(OM, 2, input$above$N/10000, "/")
      }
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Plant respiration per individual ",(g*C%.%ind^{-1})))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
    else if(type=="PlantCBalance") {
      OM = x$PlantPhotosynthesis-x$PlantRespiration
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Plant C balance ",(g*C%.%m^{-2})))
      g<-.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim)+
        geom_hline(yintercept = 0, col="gray")
      return(g)
    } 
    else if(type=="PlantCBalancePerLeaf") {
      OM = (x$PlantPhotosynthesis-x$PlantRespiration)
      if(bySpecies) {
        m1 = apply(OM,1, tapply, x$cohorts$Name, sum, na.rm=T)
        lai1 = apply(x$PlantLAIlive,1,tapply, x$cohorts$Name, sum, na.rm=T)
        OM = t(m1/lai1)
      } else {
        OM = OM/x$PlantLAIlive
        OM[x$PlantLAIlive==0] = NA
      }
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Plant C balance per leaf area ",(g*C%.%m^{-2})))
      g<-.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim)+
        geom_hline(yintercept = 0, col="gray")
      return(g)
    }
    else if(type=="PlantCBalancePerIndividual") {
      OM = (x$PlantPhotosynthesis-x$PlantRespiration)
      if(bySpecies) {
        m1 = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
        dens = as.numeric(tapply(input$above$N/10000, input$cohorts$Name, sum, na.rm=T))
        OM = sweep(m1, 2, dens, "/")
      } else {
        OM = sweep(OM, 2, input$above$N/10000, "/")
      }
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Plant C balance per individual ",(g*C%.%ind^{-1})))
      g<-.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim)+
        geom_hline(yintercept = 0, col="gray")
      return(g)
    }
    else if(type=="PlantCstorageFast") {
      OM = x$PlantCstorageFast
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, mean, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Fast C storage pool [0-1]"))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
    else if(type=="PlantCstorageSlow") {
      OM = x$PlantCstorageSlow
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, mean, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Slow C storage pool  [0-1]"))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
    else if(type=="PlantSAgrowth") {
      OM = x$PlantSAgrowth
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, mean, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Sapwood area growth ",(cm^2)))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
    else if(type=="PlantRelativeSAgrowth") {
      OM = x$PlantSAgrowth/x$PlantSA
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, mean, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Relative sapwood area growth ",(cm^2%.%cm^{-2})))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
    else if(type=="PlantSA") {
      OM = x$PlantSA
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, mean, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Sapwood area  ",(cm^2)))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    }
    else if(type=="PlantLAIlive") {
      OM = x$PlantLAIlive
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Live Leaf Area Index   ",(m^{2}%.%m^{-2})))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
    else if(type=="PlantLAIdead") {
      OM = x$PlantLAIdead
      if(bySpecies) {
        OM = t(apply(OM,1, tapply, x$cohorts$Name, sum, na.rm=T))
      } 
      if(!is.null(dates)) OM = OM[row.names(OM) %in% as.character(dates),]
      if(is.null(ylab)) ylab = expression(paste("Dead Leaf Area Index   ",(m^{2}%.%m^{-2})))
      return(.multiple_dynamics(as.matrix(OM),  xlab = xlab, ylab = ylab, ylim = ylim))
    } 
  }
}
