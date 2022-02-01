.droughtstress_sim<-function(x, index = "NDD", freq = "years", bySpecies = FALSE) {
  dates = as.Date(rownames(x$Plants$PlantStress))
  ndaysTotal = length(dates)
  date.factor = cut(dates, breaks=freq)
  if(inherits(x, c("spwb","pwb"))) {
    transpMode = x$spwbInput$control$transpirationMode
  } else {
    transpMode = x$growthInput$control$transpirationMode
  }
  if("spwbInput" %in% names(x)) input = x$spwbInput
  else input = x$growthInput
  
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
    M <- apply(x$Plants$PlantStress,2,tapply, INDEX=date.factor, ndd)
  } else if(index=="DI") {
    M <- apply(x$Plants$PlantStress,2,tapply, INDEX=date.factor, di)
  } else if(index=="ADS") {
    M <- apply(x$Plants$PlantStress,2,tapply, INDEX=date.factor, function(x) {return(mean(x, na.rm=T))})
  } else if(index=="MDS") {
    M <- apply(x$Plants$PlantStress,2,tapply, INDEX=date.factor, function(x) {return(max(x, na.rm=T))})
  } else if(index=="WSI") {
    if(transpMode=="Granier") {
      M <- apply(x$Plants$PlantPsi,2,tapply, INDEX=date.factor, wsi)
    } else {
      M <- apply(x$Plants$LeafPsi,2,tapply, INDEX=date.factor, wsi)
    }
  }
  if(is.vector(M)) {
    M = t(as.matrix(M))
    rownames(M) <- levels(date.factor)
  }
  if(bySpecies) {
    cohlai = apply(x$Plants$LAI,2,max, na.rm=T)
    cohsp = as.character(input$cohorts$Name)
    lai1 = tapply(cohlai, cohsp, sum, na.rm=T)
    m1 = t(apply(sweep(M,2,cohlai,"*"),1,tapply, cohsp, sum, na.rm=T))
    if(length(lai1)==1) {
      m1 = t(m1)
    } 
    M = sweep(m1,2,lai1,"/")
    colnames(M) = names(lai1)
  }
  ncases = table(date.factor)
  M = M[ncases>0, ,drop = FALSE]
  return(M)
}
droughtStress<-function(x, index = "NDD", freq = "years", bySpecies = FALSE, draw = TRUE) {
  if(!inherits(x,c("spwb","pwb","growth", "fordyn"))) {
    stop("'x' should be of class 'spwb', 'pwb', 'growth' or 'fordyn'")
  }
  index = match.arg(index,c("NDD","DI", "ADS", "MDS","WSI"))  
  
  if(inherits(x, c("spwb","pwb", "growth"))) {
     M = .droughtstress_sim(x = x, index = index, freq = freq, bySpecies = bySpecies)
  } else {
    vec<-vector("list", length(x$GrowthResults))
    for(i in 1:length(x$GrowthResults)) {
      vec[[i]] <- .droughtstress_sim(x = x$GrowthResults[[i]], index = index, freq = freq, bySpecies = bySpecies)
    }
    M = .mergeVectorOfMatrices(vec)
  }
  if(draw) {
    if(index=="NDD") {
      g<-.multiple_dynamics(M, ylab = "Number of days with DS > 0.5")
    }
    else if(index=="DI") {
      g<-.multiple_dynamics(M, ylab = "Average drought stress over 0.5")
    }
    else if(index=="ADS") {
      g<-.multiple_dynamics(M, ylab = "Average drought stress")
    }
    else if(index=="MDS") {
      g<-.multiple_dynamics(M, ylab = "Maximum drought stress")
    }
    else if(index=="WSI") {
      g<-.multiple_dynamics(M, ylab = "Water stress integral")
    }
    return(g)
  } else {
    return(M)
  }
}
