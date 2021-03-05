summary.spwb<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){  
  dates = as.Date(rownames(object$WaterBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  input = object$spwbInput
  if(output=="WaterBalance") OM = object$WaterBalance
  else if(output=="Soil") OM = object$Soil
  else if(output=="Stand") OM = object$Stand
  else if(output=="PlantPhotosynthesis") OM = object$PlantPhotosynthesis
  else if(output=="PlantGrossPhotosynthesis") OM = object$PlantGrossPhotosynthesis
  else if(output=="PlantNetPhotosynthesis") OM = object$PlantNetPhotosynthesis
  else if(output=="PlantTranspiration") OM = object$PlantTranspiration
  else if(output=="EnergyBalance") OM = object$EnergyBalance
  else if(output=="Temperature") OM = object$Temperature
  else if(output=="TemperatureLayers") OM = object$TemperatureLayers
  else if(output=="PlantLAI") {
    OM = object$Plants$LAI
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
    } 
  } 
  else if(output=="subdaily") {
    stop("Cannot summarize subdaily output")
  }
  else {
    OM = object$Plants[[output]]
    if(bySpecies) {
      lai1 = t(apply(object$Plants$LAI,1, tapply, input$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(object$Plants$LAI * OM,1, tapply, input$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
  }
  
  if(ncol(OM)==length(date.factor) && nrow(OM)==1) OM = t(OM)
  
  #Perform summary at the desired temporal scale
  M <- apply(OM,2,tapply, INDEX=date.factor, FUN)
  
  if(is.vector(M)) {
    M = t(as.matrix(M))
    rownames(M) <- levels(date.factor)
  }
  # if(sum(is.na(M[nrow(M), drop=FALSE]))==ncol(M)) M = M[-nrow(M), drop=FALSE] #Remove empty row
  ncases = table(date.factor)
  M = M[ncases>0, ,drop = FALSE]
  
  return(M)
}

summary.pwb<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){  
  summary.spwb(object = object, freq = freq, output = output, FUN = FUN, bySpecies = bySpecies, ...)
}