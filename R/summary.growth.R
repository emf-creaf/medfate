summary.growth<-function(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...){  
  dates = as.Date(rownames(object$WaterBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  OM = object[[output]]

  #Aggregate results by species
  if(bySpecies) {
    lai1 = t(apply(object$PlantLAI,1, tapply, object$cohorts$Name, sum))
    m1 = t(apply(object$PlantLAI * OM,1, tapply, object$cohorts$Name, sum))
    OM = m1/lai1
    OM[lai1==0] = NA
  } 
  if(ncol(OM)==length(date.factor) && nrow(OM)==1) OM = t(OM)
  
  #Perform summary at the desired temporal scale
  M <- apply(OM,2,tapply, INDEX=date.factor, FUN)
  if(sum(is.na(M[nrow(M), drop=FALSE]))==ncol(M)) M = M[-nrow(M), drop=FALSE] #Remove empty row
  return(M)
}