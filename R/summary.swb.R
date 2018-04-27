summary.swb<-function(object, freq="years", output="DailyBalance", FUN=sum, bySpecies = FALSE, ...){  
  dates = as.Date(rownames(object$DailyBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  speciesOutput = FALSE
  if(output=="DailyBalance") OM = object$DailyBalance
  else if(output=="SoilWaterBalance") OM = object$SoilWaterBalance
  else if(output=="TemperatureBalance") OM = object$TemperatureBalance
  else if(output=="PlantLAI") {
    OM = object$PlantLAI
    speciesOutput = TRUE
  } else if(output=="PlantPsi") {
    OM = object$PlantPsi
    speciesOutput = TRUE
  }
  else if(output=="PlantStress") {
    OM = object$PlantStress
    speciesOutput = TRUE
  }
  else if(output=="PlantTranspiration") {
    OM = object$PlantTranspiration
    speciesOutput = TRUE
  }
  else if(output=="PlantPhotosynthesis") {
    OM = object$PlantPhotosynthesis
    speciesOutput = TRUE
  }

  #Aggregate results by species
  if(speciesOutput && bySpecies) {
    lai1 = t(apply(object$PlantLAI,1, tapply, object$cohorts$Name, sum))
    m1 = t(apply(object$PlantLAI * OM,1, tapply, object$cohorts$Name, sum))
    OM = m1/lai1
    OM[lai1==0] = NA
  } 
  
  #Perform summary at the desired temporal scale
  ncoh = ncol(OM)
  M = data.frame(matrix(0,nrow=length(levels(date.factor)), ncol=ncol(OM)))
  row.names(M)<-levels(date.factor)
  if(class(OM)=="matrix") names(M)<-colnames(OM)
  else if(class(OM)=="data.frame") names(M)<-names(OM)
  for(c in 1:ncol(OM)){
    M[,c] <-tapply(OM[,c], INDEX=date.factor,FUN=FUN)
  }
  if(sum(is.na(M[nrow(M),]))==ncol(M)) M = M[-nrow(M),] #Remove empty row
  return(M)
}