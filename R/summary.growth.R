summary.growth<-function(object, freq="years", output="DailyBalance", FUN=sum,...){  
  dates = as.Date(rownames(object$DailyBalance))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  OM = object[[output]]

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