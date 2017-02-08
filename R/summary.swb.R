summary.swb<-function(object, freq="years", output="DailyBalance", FUN=sum,...){  
  dates = as.Date(rownames(object$DailyBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  if(output=="DailyBalance") OM = object$DailyBalance
  else if(output=="PlantStress") OM = object$PlantStress
  else if(output=="PlantTranspiration") OM = object$PlantTranspiration
  else if(output=="SoilWaterBalance") OM = object$SoilWaterBalance
  
  ncoh = ncol(OM)
  M = data.frame(matrix(0,nrow=length(levels(date.factor)), ncol=ncol(OM)))
  rownames(M)<-levels(date.factor)
  colnames(M)<-names(OM)
  for(c in 1:ncol(OM)){
    M[,c] <-tapply(OM[,c], INDEX=date.factor,FUN=FUN)
  }
  if(sum(is.na(M[nrow(M),]))==ncol(M)) M = M[-nrow(M),] #Remove empty row
  return(M)
}