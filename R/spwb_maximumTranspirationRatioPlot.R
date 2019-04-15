spwb_maximumTranspirationRatioPlot<-function(x, soil,  meteo, latitude, elevation, ndays = 100, 
                                        LAI_seq = c(0.1,0.25, seq(0.5, 8, by=0.5))) {
  meteo = meteo[meteo$Precipitation==0, ] #Exclude days with precipitation
  meteo = meteo[!is.na(meteo$PET), ] #Exclude days without PET
  PET = meteo$PET
  PET_cut = cut(PET,breaks = seq(0, max(PET), length.out = 20))
  nPET = table(PET_cut)
  PETw = as.numeric(1/nPET[PET_cut])
  PETw = PETw/sum(PETw)
  s = sample(1:nrow(meteo),ndays,replace = TRUE, prob = PETw)
  s = s[order(PET[s])]
  meteo = meteo[s,]
  PET = meteo$PET
  
  ncoh = nrow(x$above)
  ndays = nrow(meteo)
  nlai = length(LAI_seq)
  
 
  pwb_list = vector("list", ncoh)
  Tmax_list = vector("list", ncoh)
  names(Tmax_list) = row.names(x$above)

  for(i in 1:ncoh) {
    cat(paste0("\n Cohort: ", row.names(x$above)[i],"\n"))
    
    xIni = x
    xIni$above$LAI_live[-i] = 0 #Set LAI of other cohors to zero
    xIni$above$LAI_expanded[-i] = 0
    s = soil
    xIni$control$cavitationRefill = TRUE
    xIni$control$verbose = FALSE
    Tmax = matrix(NA, nrow=ndays, ncol = nlai)  
    colnames(Tmax) = LAI_seq
    rownames(Tmax) = row.names(meteo)

    pwb_res = vector("list", nlai)
    pb = txtProgressBar(0, nlai, style=3)
    for(j in 1:nlai) {
      LAIComp = LAI_seq[j]
      xlai = xIni
      xlai$above$LAI_live = xlai$above$LAI_live*LAIComp/sum(xlai$above$LAI_live)
      xlai$above$LAI_expanded = xlai$above$LAI_live
      spwb_resetInputs(xlai, s)
      W = matrix(1, nrow=ndays, ncol = length(s$W))
      pwb_res[[j]] = pwb(xlai,s,meteo,W, latitude = latitude, elevation = elevation)
      Tmax[,j] = pwb_res[[j]]$WaterBalance$Transpiration
      setTxtProgressBar(pb, j)
    }
    Tmax_list[[i]] = Tmax
    pwb_list[[i]] = pwb_res
  }
  
  Tmax_all = matrix(NA, nrow=ndays, ncol = nlai)  
  colnames(Tmax_all) = LAI_seq
  rownames(Tmax_all) = row.names(meteo)
  
  xIni = x
  s = soil
  xIni$control$cavitationRefill = TRUE
  xIni$control$verbose = FALSE
  
  cat(paste0("\nAll cohorts. \n"))
  pwb_res = vector("list", nlai)
  pb = txtProgressBar(0, nlai, style=3)
  for(j in 1:nlai) {
    LAIComp = LAI_seq[j]
    xlai = xIni
    xlai$above$LAI_live = xlai$above$LAI_live*LAIComp/sum(xlai$above$LAI_live)
    xlai$above$LAI_expanded = xlai$above$LAI_live
    spwb_resetInputs(xlai, s)
    W = matrix(1, nrow=ndays, ncol = length(s$W))
    pwb_res[[j]] = pwb(xlai,s,meteo,W, latitude = latitude, elevation = elevation)
    Tmax_all[,j] = pwb_res[[j]]$WaterBalance$Transpiration
    setTxtProgressBar(pb, j)
  }
  
  TmaxPETGranier = -0.006*(LAI_seq^2)+0.134*LAI_seq+0.036
  plot(LAI_seq, TmaxPETGranier, type="l", col="gray", lwd=2, 
       xlab = "Stand's Leaf Area Index", ylab = "Tmax/PET", ylim=c(0,1))
  for(i in 1:length(Tmax_list)) {
    Tmax = Tmax_list[[i]]
    TmaxRatio = sweep(Tmax,1,PET,"/")
    m = apply(TmaxRatio, 2, median, na.rm=T)
    lines(LAI_seq, m, col=i+1, lwd=1, lty=i)
  }
  TmaxRatio = sweep(Tmax_all,1,PET,"/")
  m = apply(TmaxRatio, 2, median, na.rm=T)
  lines(LAI_seq, m, col="black", lwd=2, lty=1)
  legend("topleft", legend=c("Granier's empirical equation", "Whole stand"), lwd=2, lty=1, col=c("gray", "black"),bty="n", cex =1)
  legend("bottomright", legend=c(row.names(x$cohorts)), lwd=1, lty=1:ncoh, col=2:(ncoh+1),bty="n", cex=1)
}