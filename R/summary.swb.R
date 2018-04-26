summary.swb<-function(object, freq="years", output="DailyBalance", FUN=sum,...){  
  dates = as.Date(rownames(object$DailyBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  if(output=="DailyBalance") OM = object$DailyBalance
  else if(output=="PlantLAI") OM = object$PlantLAI
  else if(output=="PlantPsi") OM = object$PlantPsi
  else if(output=="PlantStress") OM = object$PlantStress
  else if(output=="PlantTranspiration") OM = object$PlantTranspiration
  else if(output=="PlantPhotosynthesis") OM = object$PlantPhotosynthesis
  else if(output=="SoilWaterBalance") OM = object$SoilWaterBalance
  else if(output=="TemperatureBalance") OM = object$TemperatureBalance
  
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