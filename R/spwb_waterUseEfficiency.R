spwb_waterUseEfficiency<-function(x, type = "Plant An/E", leaves = "average", freq="days", draw = TRUE) {
  if(!("spwb" %in% class(x)) && !("pwb" %in% class(x))) {
    stop("'x' should be of class 'spwb' or 'pwb'")
  }
  type = match.arg(type, c("Leaf Ci", "Leaf iWUE", "Plant An/E", "Stand An/E"))
  if(type=="Leaf iWUE") {
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
        iWUEdays[i,] = rowSums(iwueinst*sl$An, na.rm=T)/rowSums(sl$An, na.rm=T) #Photosynthesis-weighted iWUE
      }
      else if(leaves =="shade") {
        iwueinst = sh$An/sh$GW
        iwueinst[iwueinst<0] = 0
        iWUEdays[i,] = rowSums(iwueinst*sh$An, na.rm=T)/rowSums(sh$An, na.rm=T) #Photosynthesis-weighted iWUE
      }
      else {
        iwueinst_sl = sl$An/sl$GW
        iwueinst_sh = sh$An/sh$GW
        iwueinst = ((iwueinst_sl*sl_lai) + (iwueinst_sh*sh_lai))/(sl_lai+sh_lai)
        iwueinst[iwueinst<0] = 0
        An_tot = sh$An+sl$An
        iWUEdays[i,] = rowSums(iwueinst*An_tot, na.rm=T)/rowSums(An_tot, na.rm=T) #Photosynthesis-weighted iWUE
      }
    }
    if(freq=="days") {
      res = iWUEdays
    } else {
      
      date.factor = cut(dates, breaks=freq)
      
      Andays = x$PlantPhotosynthesis
      #Perform summary at the desired temporal scale (weighted average of iWUE with An as weights)
      sumiWUEAn <- apply(iWUEdays*Andays,2,tapply, INDEX=date.factor, sum, na.rm=T)
      sumAn <- apply(Andays,2,tapply, INDEX=date.factor, sum, na.rm=T)
      M = sumiWUEAn/sumAn
      if(is.vector(M)) {
        M = t(as.matrix(M))
        rownames(M) <- levels(date.factor)
      }
      ncases = table(date.factor)
      M = M[ncases>0, ,drop = FALSE]
      res = M
    }
  }
  else if(type=="Leaf Ci") {
    if(x$spwbInput$control$transpirationMode != "Sperry") {
      stop("Ci can only be calculated with transpirationMode = 'Sperry'")
    }
    if(!x$spwbInput$control$subdailyResults) {
      stop("Ci can only be calculated with subdailyResults = TRUE")
    }
    leaves = match.arg(leaves, c("average", "sunlit", "shade"))
    sd = x$subdaily
    ndays = length(sd)
    coh = x$spwbInput$cohorts
    ncoh = nrow(coh)
    dates = as.Date(names(sd))
    Cidays = matrix(NA, nrow=ndays, ncol=ncoh)
    rownames(Cidays)= as.character(dates)
    colnames(Cidays) = rownames(coh)
    for(i in 1:ndays) {
      sl = sd[[i]]$PlantsInst$SunlitLeaves
      sh = sd[[i]]$PlantsInst$ShadeLeaves
      sl_lai = sd[[i]]$SunlitLeaves$LAI
      sh_lai = sd[[i]]$ShadeLeaves$LAI
      if(leaves =="sunlit") {
        ciinst = sl$Ci
        Cidays[i,] = rowSums(ciinst*sl$An, na.rm=T)/rowSums(sl$An, na.rm=T) #Photosynthesis-weighted iWUE
      }
      else if(leaves =="shade") {
        ciinst = sh$Ci
        Cidays[i,] = rowSums(ciinst*sh$An, na.rm=T)/rowSums(sh$An, na.rm=T) #Photosynthesis-weighted iWUE
      }
      else {
        ciinst_sl = sl$Ci
        ciinst_sh = sh$Ci
        ciinst = ((ciinst_sl*sl_lai) + (ciinst_sh*sh_lai))/(sl_lai+sh_lai)
        An_tot = sh$An+sl$An
        Cidays[i,] = rowSums(ciinst*An_tot, na.rm=T)/rowSums(An_tot, na.rm=T) #Photosynthesis-weighted iWUE
      }
    }
    if(freq=="days") {
      res = Cidays
    } else {
      
      date.factor = cut(dates, breaks=freq)
      
      Andays = x$PlantPhotosynthesis
      #Perform summary at the desired temporal scale (weighted average of iWUE with An as weights)
      sumCiAn <- apply(Cidays*Andays,2,tapply, INDEX=date.factor, sum, na.rm=T)
      sumAn <- apply(Andays,2,tapply, INDEX=date.factor, sum, na.rm=T)
      M = sumCiAn/sumAn
      if(is.vector(M)) {
        M = t(as.matrix(M))
        rownames(M) <- levels(date.factor)
      }
      ncases = table(date.factor)
      M = M[ncases>0, ,drop = FALSE]
      res = M
    }
  }
  else if(type =="Plant An/E") {
    if(freq=="days") {
      res = x$PlantPhotosynthesis/x$PlantTranspiration
    } else {
      pt = summary(x, freq=freq, output="PlantTranspiration", FUN=sum, na.rm=T)
      pp = summary(x, freq=freq, output="PlantPhotosynthesis", FUN=sum, na.rm=T)
      res = pp/pt
    }
  }
  else if(type =="Stand An/E") {
    if(freq=="days") {
      res = rowSums(x$PlantPhotosynthesis)/rowSums(x$PlantTranspiration)
    } else {
      pt = summary(x, freq=freq, output="PlantTranspiration", FUN=sum, na.rm=T)
      pp = summary(x, freq=freq, output="PlantPhotosynthesis", FUN=sum, na.rm=T)
      res = rowSums(pp)/rowSums(pt)
    }
  }
  if(!draw) {
    return(res)
  } else {
    if(type=="Stand An/E") {
      g<-.single_dynamics(res, ylab = "Stand An/E (gC/L)")
    } 
    else if(type=="Plant An/E") {
      g<-.multiple_dynamics(res, ylab = "Plant An/E (gC/L)")
    }
    else if(type %in% c("Leaf iWUE", "Leaf Ci")) {
      g<-.multiple_dynamics(res, ylab = type)
    }
    return(g)
  }
}