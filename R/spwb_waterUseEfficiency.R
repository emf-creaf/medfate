spwb_waterUseEfficiency<-function(x, type = "Plant An/E", leaves = "average", freq="days", draw = TRUE,
                                  ylim=NULL) {
  if(!("spwb" %in% class(x)) && !("pwb" %in% class(x))) {
    stop("'x' should be of class 'spwb' or 'pwb'")
  }
  type = match.arg(type, c("Leaf Ci", "Leaf iWUE", "Plant An/E", "Stand An/E"))
  
  if("Photosynthesis" %in% names(x$Plants)) {
    Andays = x$Plants$Photosynthesis
  } else {
    Andays = x$Plants$NetPhotosynthesis
  }
  
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
      sl = sd[[i]]$SunlitLeavesInst
      sh = sd[[i]]$ShadeLeavesInst
      sl_lai = sd[[i]]$SunlitLeaves$LAI
      sh_lai = sd[[i]]$ShadeLeaves$LAI
      an_sl = sl$An
      an_sl[an_sl<0] = 0
      an_sh = sh$An
      an_sh[an_sh<0] = 0
      if(leaves =="sunlit") {
        iwueinst = an_sl/sl$GW
        iwueinst[iwueinst<0] = 0
        iWUEdays[i,] = rowSums(iwueinst*an_sl, na.rm=T)/rowSums(an_sl, na.rm=T) #Photosynthesis-weighted iWUE
      }
      else if(leaves =="shade") {
        iwueinst = an_sh/sh$GW
        iwueinst[iwueinst<0] = 0
        iWUEdays[i,] = rowSums(iwueinst*an_sh, na.rm=T)/rowSums(an_sh, na.rm=T) #Photosynthesis-weighted iWUE
      }
      else {
        iwueinst_sl = an_sl/sl$GW
        iwueinst_sh = an_sh/sh$GW
        iwueinst = ((iwueinst_sl*sl_lai) + (iwueinst_sh*sh_lai))/(sl_lai+sh_lai)
        iwueinst[iwueinst<0] = 0
        an_tot = an_sl + an_sh
        iWUEdays[i,] = rowSums(iwueinst*an_tot, na.rm=T)/rowSums(an_tot, na.rm=T) #Photosynthesis-weighted iWUE
      }
    }
    if(freq=="days") {
      res = iWUEdays
    } else {
      
      date.factor = cut(dates, breaks=freq)
      
      Andays[Andays<0] = 0
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
      sl = sd[[i]]$SunlitLeavesInst
      sh = sd[[i]]$ShadeLeavesInst
      sl_lai = sd[[i]]$SunlitLeaves$LAI
      sh_lai = sd[[i]]$ShadeLeaves$LAI
      an_sl = sl$An
      an_sl[an_sl<0] = 0
      an_sh = sh$An
      an_sh[an_sh<0] = 0
      if(leaves =="sunlit") {
        ciinst = sl$Ci
        Cidays[i,] = rowSums(ciinst*an_sl, na.rm=T)/rowSums(an_sl, na.rm=T) #Photosynthesis-weighted iWUE
      }
      else if(leaves =="shade") {
        ciinst = sh$Ci
        Cidays[i,] = rowSums(ciinst*an_sh, na.rm=T)/rowSums(an_sh, na.rm=T) #Photosynthesis-weighted iWUE
      }
      else {
        ciinst_sl = sl$Ci
        ciinst_sh = sh$Ci
        ciinst = ((ciinst_sl*sl_lai) + (ciinst_sh*sh_lai))/(sl_lai+sh_lai)
        an_tot = an_sh+an_sl
        Cidays[i,] = rowSums(ciinst*an_tot, na.rm=T)/rowSums(an_tot, na.rm=T) #Photosynthesis-weighted iWUE
      }
    }
    if(freq=="days") {
      res = Cidays
    } else {
      
      date.factor = cut(dates, breaks=freq)
      
      Andays[Andays<0] = 0
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
    if("Photosynthesis" %in% names(x$Plants)) x$Plants$Photosynthesis[x$Plants$Photosynthesis<0] = 0
    else x$Plants$NetPhotosynthesis[x$Plants$NetPhotosynthesis<0] = 0
    x$Plants$Transpiration[x$Plants$Transpiration<0] = 0
    if(freq=="days") {
      if("Photosynthesis" %in% names(x$Plants)) res = x$Plants$Photosynthesis/x$Plants$Transpiration
      else res = x$Plants$NetPhotosynthesis/x$Plants$Transpiration
    } else {
      pt = summary(x, freq=freq, output="Transpiration", FUN=sum, na.rm=T)
      pp = summary(x, freq=freq, output="Photosynthesis", FUN=sum, na.rm=T)
      res = pp/pt
    }
  }
  else if(type =="Stand An/E") {
    if("Photosynthesis" %in% names(x$Plants)) x$Plants$Photosynthesis[x$Plants$Photosynthesis<0] = 0
    else x$Plants$NetPhotosynthesis[x$Plants$NetPhotosynthesis<0] = 0
    x$Plants$Transpiration[x$Plants$Transpiration<0] = 0
    if(freq=="days") {
      if("Photosynthesis" %in% names(x$Plants)) res = rowSums(x$Plants$Photosynthesis, na.rm=T)/rowSums(x$Plants$Transpiration, na.rm=T)
      else res = rowSums(x$Plants$NetPhotosynthesis, na.rm=T)/rowSums(x$Plants$Transpiration, na.rm=T)
    } else {
      pt = summary(x, freq=freq, output="Transpiration", FUN=sum, na.rm=T)
      pp = summary(x, freq=freq, output="Photosynthesis", FUN=sum, na.rm=T)
      res = rowSums(pp, na.rm=T)/rowSums(pt, na.rm=T)
    }
  }
  if(!draw) {
    return(res)
  } else {
    if(type=="Stand An/E") {
      g<-.single_dynamics(res, ylab = "Stand An/E (gC/L)", ylim = ylim)
    } 
    else if(type=="Plant An/E") {
      g<-.multiple_dynamics(res, ylab = "Plant An/E (gC/L)", ylim = ylim)
    }
    else if(type %in% c("Leaf iWUE", "Leaf Ci")) {
      if(type=="Leaf iWUE" && leaves == "sunlit") ylab = expression(paste("Sunlit leaf iWUE  (",mu%.%mol%.%mol^{-1},")"))
      if(type=="Leaf iWUE" && leaves == "shade") ylab = expression(paste("Shade leaf iWUE  (",mu%.%mol%.%mol^{-1},")"))
      if(type=="Leaf iWUE" && leaves == "average") ylab = expression(paste("Average leaf iWUE  (",mu%.%mol%.%mol^{-1},")"))
      if(type=="Leaf Ci" && leaves == "sunlit") ylab = expression(paste("Sunlit leaf Ci  ",(ppm)))
      if(type=="Leaf Ci" && leaves == "shade") ylab = expression(paste("Shade leaf Ci  ",(ppm)))
      if(type=="Leaf Ci" && leaves == "average") ylab = expression(paste("Average leaf Ci  ",(ppm)))
      g<-.multiple_dynamics(res, ylab = ylab, ylim = ylim)
    }
    return(g)
  }
}