summary.spwb<-function(object, freq="years", output="DailyBalance", FUN=sum, bySpecies = FALSE, ...){  
  dates = as.Date(rownames(object$DailyBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  if(output=="DailyBalance") OM = object$DailyBalance
  else if(output=="SoilWaterBalance") OM = object$SoilWaterBalance
  else if(output=="TemperatureBalance") OM = object$TemperatureBalance
  else if(output=="PlantLAI") {
    OM = object$PlantLAI
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, object$cohorts$Name, sum, na.rm=T))
    } 
  } else if(output=="PlantPsi") {
    OM = object$PlantPsi
    if(bySpecies) {
      lai1 = t(apply(object$PlantLAI,1, tapply, object$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(object$PlantLAI * OM,1, tapply, object$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
  }
  else if(output=="PlantStress") {
    OM = object$PlantStress
    if(bySpecies) {
      lai1 = t(apply(object$PlantLAI,1, tapply, object$cohorts$Name, sum, na.rm=T))
      m1 = t(apply(object$PlantLAI * OM,1, tapply, object$cohorts$Name, sum, na.rm=T))
      OM = m1/lai1
      OM[lai1==0] = NA
    } 
  }
  else if(output=="PlantTranspiration") {
    OM = object$PlantTranspiration
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, object$cohorts$Name, sum, na.rm=T))
    } 
  }
  else if(output=="PlantPhotosynthesis") {
    OM = object$PlantPhotosynthesis
    if(bySpecies) {
      OM = t(apply(OM,1, tapply, object$cohorts$Name, sum, na.rm=T))
    } 
  }

  if(ncol(OM)==length(date.factor) && nrow(OM)==1) OM = t(OM)
  
  #Perform summary at the desired temporal scale
  M <- apply(OM,2,tapply, INDEX=date.factor, FUN)
  if(sum(is.na(M[nrow(M), drop=FALSE]))==ncol(M)) M = M[-nrow(M), drop=FALSE] #Remove empty row
  return(M)
}