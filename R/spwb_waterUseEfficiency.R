spwb_waterUseEfficiency<-function(x, type = "An/E", leaves = "average", freq="days") {
  if(!("spwb" %in% class(x)) && !("pwb" %in% class(x))) {
    stop("'x' should be of class 'spwb' or 'pwb'")
  }
  type = match.arg(type, c("iWUE", "An/E"))
  if(type=="iWUE") {
    if(x$spwbInput$control$transpirationMode != "Sperry") {
      stop("iWUE can only be calculated with transpirationMode = 'Sperry'")
    }
    if(!x$spwbInput$control$subdailyResults) {
      stop("iWUE can only be calculated with subdailyResults = TRUE")
    }
    leaves = match.arg(leaves, c("average", "sunlit", "shade"))
    sd = x$subdaily
    ndays = length(sd)
    coh = x$spwbInput$cohorts
    ncoh = nrow(coh)
    dates = as.Date(names(sd))
    iWUEdays = matrix(NA, nrow=ndays, ncol=ncoh)
    rownames(iWUEdays)= as.character(dates)
    colnames(iWUEdays) = rownames(coh)
    for(i in 1:ndays) {
      sl = sd[[i]]$PlantsInst$SunlitLeaves
      sh = sd[[i]]$PlantsInst$ShadeLeaves
      sl_lai = sd[[i]]$SunlitLeaves$LAI
      sh_lai = sd[[i]]$ShadeLeaves$LAI
      if(leaves =="sunlit") {
        iwueinst = sl$An/sl$GW
        iwueinst[iwueinst<0] = 0
        iWUEdays[i,] = rowSums(iwueinst*sl$An)/rowSums(sl$An) #Photosynthesis-weighted iWUE
      }
      if(leaves =="shade") {
        iwueinst = sh$An/sh$GW
        iwueinst[iwueinst<0] = 0
        iWUEdays[i,] = rowSums(iwueinst*sh$An)/rowSums(sh$An) #Photosynthesis-weighted iWUE
      }
      else {
        iwueinst_sl = sl$An/sl$GW
        iwueinst_sh = sh$An/sh$GW
        iwueinst = ((iwueinst_sl*sl_lai) + (iwueinst_sh*sh_lai))/(sl_lai+sh_lai)
        iwueinst[iwueinst<0] = 0
        An_tot = sh$An+sl$An
        iWUEdays[i,] = rowSums(iwueinst*An_tot)/rowSums(An_tot) #Photosynthesis-weighted iWUE
      }
    }
    if(freq=="days") {
      return(iWUEdays)
    } else {
      
      date.factor = cut(dates, breaks=freq)
      
      Andays = x$PlantPhotosynthesis
      #Perform summary at the desired temporal scale (weighted average of iWUE with An as weights)
      sumiWUEAn <- apply(iWUEdays*Andays,2,tapply, INDEX=date.factor, sum)
      sumAn <- apply(Andays,2,tapply, INDEX=date.factor, sum)
      M = sumiWUEAn/sumAn
      if(is.vector(M)) {
        M = t(as.matrix(M))
        rownames(M) <- levels(date.factor)
      }
      ncases = table(date.factor)
      M = M[ncases>0, ,drop = FALSE]
      return(M)
    }
  }
  else if(type =="An/E") {
    if(freq=="days") {
      return(x$PlantPhotosynthesis/x$PlantTranspiration)
    } else {
      pt = summary(x, freq=freq, output="PlantTranspiration", FUN=sum)
      pp = summary(x, freq=freq, output="PlantPhotosynthesis", FUN=sum)
      return(pp/pt)
    }
  }
}