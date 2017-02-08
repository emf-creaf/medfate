translateSpeciesCodes<-function(x, SFIcodes) {  
  lsfi = strsplit(SFIcodes,"[,./]")
  sfiNumCod = unique(as.numeric(unlist(lsfi)))
  repVect = rep(NA,max(sfiNumCod))
  for(i in 1:length(lsfi)) {
    cv = as.numeric(lsfi[[i]])
    for(ch in cv) {
      repVect[ch] = (i-1) #Species indices start from 0 in medfate
    }
  }
  return(repVect[as.numeric(x$Especie)])
}