spwb.stress<-function(x, index = "NDD", freq = "years", bySpecies = FALSE) {
  index = match.arg(index,c("NDD","DI", "ADS", "MDS","WSI"))  
  dates = as.Date(rownames(x$PlantStress))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  
  transpMode = x$spwbInput$control$transpirationMode
  
  ndd<-function(dds) { # De Caceres et al (2015)
    return(sum(dds>0.5))
  }
  di<-function(dds) { # De Caceres et al (2015)
    return(sum(pmax(0,(dds-0.5)/0.5))/length(dds))
  }
  wsi<-function(lwp) { # Myers (1988)
    c = max(lwp, na.rm=T)
    return(abs(sum(lwp-c, na.rm=T)))
  }
  if(index=="NDD") {
    M <- apply(x$PlantStress,2,tapply, INDEX=date.factor, ndd)
  } else if(index=="DI") {
    M <- apply(x$PlantStress,2,tapply, INDEX=date.factor, di)
  } else if(index=="ADS") {
    M <- apply(x$PlantStress,2,tapply, INDEX=date.factor, function(x) {return(mean(x, na.rm=T))})
  } else if(index=="MDS") {
    M <- apply(x$PlantStress,2,tapply, INDEX=date.factor, function(x) {return(max(x, na.rm=T))})
  } else if(index=="WSI") {
    if(transpMode=="Simple") {
      M <- apply(x$PlantPsi,2,tapply, INDEX=date.factor, wsi)
    } else {
      M <- apply(x$LeafPsi,2,tapply, INDEX=date.factor, wsi)
    }
  }
  ncases = table(date.factor)
  M = M[ncases>0, ,drop = FALSE]
  if(bySpecies) {
    cohlai = apply(x$PlantLAI,2,max, na.rm=T)
    cohsp = as.character(x$spwbInput$cohorts$Name)
    lai1 = tapply(cohlai, cohsp, sum, na.rm=T)
    m1 = t(apply(sweep(M,2,cohlai,"*"),1,tapply, cohsp, sum, na.rm=T))
    M = sweep(m1,2,lai1,"/")
  }
  return(M)
}
