summary.growth<-function(object, freq="years", output="DailyBalance", FUN=sum, bySpecies = FALSE, ...){  
  dates = as.Date(rownames(object$DailyBalance))
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