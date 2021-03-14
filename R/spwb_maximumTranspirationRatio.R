spwb_maximumTranspirationRatioModels<-function(x, meteo, latitude, elevation, slope, aspect,
                                         ndays = 100, 
                                         LAI_seq = c(0.1,0.25, seq(0.5, 7, by=0.5)),
                                         draw = TRUE) {
  
  #Exclude days with precipitation
  meteo = meteo[meteo$Precipitation==0, ] 
  
  #Calculate PET using penman
  if("PET" %in% names(meteo)) meteo$PET = NULL
  PET <- numeric(nrow(meteo))
  
  cat(paste0("\n Calculating PET...\n"))
  for (i in 1:length(meteo[['MinTemperature']])) {
    PET[i] <- meteoland::penman(
      latrad = latitude*pi/180,
      elevation = elevation,
      slorad = slope*pi/180,
      asprad = aspect*pi/180,
      J = meteoland::radiation_dateStringToJulianDays(row.names(meteo)[i]),
      Tmin = meteo[['MinTemperature']][i],
      Tmax = meteo[['MaxTemperature']][i],
      RHmin = meteo[['MinRelativeHumidity']][i],
      RHmax = meteo[['MaxRelativeHumidity']][i],
      R_s = meteo[['Radiation']][i],
      u = meteo[['WindSpeed']][i],
      z = 2
    )
  }
  
  #Exclude days without PET
  meteo = meteo[!is.na(PET), ]
  PET = PET[!is.na(PET)]
  
  #Subsample days from PET
  PET_cut = cut(PET,breaks = seq(0, max(PET), length.out = 20))
  nPET = table(PET_cut)
  PETw = as.numeric(1/nPET[PET_cut])
  PETw = PETw/sum(PETw)
  s = sample(1:nrow(meteo),ndays,replace = TRUE, prob = PETw)
  s = s[order(PET[s])]
  meteo = meteo[s,]
  PET =PET[s]
  
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
    xIni$control$cavitationRefill = "total"
    xIni$control$verbose = FALSE
    xIni$control$leafPhenology = FALSE
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
      resetInputs(xlai)
      W = matrix(1, nrow=ndays, ncol = length(xlai$soil$W))
      pwb_res[[j]] = pwb(xlai, meteo,W, latitude = latitude, elevation = elevation, slope = slope, aspect = aspect)
      Tmax[,j] = pwb_res[[j]]$WaterBalance$Transpiration
      setTxtProgressBar(pb, j)
    }
    Tmax_list[[i]] = Tmax
    pwb_list[[i]] = pwb_res
  }
  # Tmax_all = matrix(NA, nrow=ndays, ncol = nlai)  
  # colnames(Tmax_all) = LAI_seq
  # rownames(Tmax_all) = row.names(meteo)
  # 
  #   
  # xIni = x
  # xIni$control$cavitationRefill = "total"
  # xIni$control$verbose = FALSE
  
  # cat(paste0("\nAll cohorts. \n"))
  # pwb_res = vector("list", nlai)
  # pb = txtProgressBar(0, nlai, style=3)
  # for(j in 1:nlai) {
  #   LAIComp = LAI_seq[j]
  #   xlai = xIni
  #   xlai$above$LAI_live = xlai$above$LAI_live*LAIComp/sum(xlai$above$LAI_live)
  #   xlai$above$LAI_expanded = xlai$above$LAI_live
  #   resetInputs(xlai)
  #   W = matrix(1, nrow=ndays, ncol = length(xlai$soil$W))
  #   pwb_res[[j]] = pwb(xlai,meteo,W, latitude = latitude, elevation = elevation)
  #   Tmax_all[,j] = pwb_res[[j]]$WaterBalance$Transpiration
  #   setTxtProgressBar(pb, j)
  # }
  mods = vector("list", ncoh)
  names(mods) <- row.names(x$cohorts)
  for(i in 1:ncoh) {
    Tmax = Tmax_list[[i]] 
    TmaxRatio = sweep(Tmax,1,PET,"/")
    Tmaxratiovec = as.vector(TmaxRatio)
    laivec = as.numeric(as.character(gl(n=ncol(Tmax),k = nrow(Tmax), labels = colnames(Tmax))))
    df = data.frame(y=Tmaxratiovec, LAI = laivec)
    df = df[(df$y > 0.0) & (df$y < 1.0),]
    mods[[i]] <- glm(y ~ -1 + LAI + I(LAI^2), 
                     start = c(0.134,-0.006),
                     data =df, family=Gamma(link="identity"))
  }
  # TmaxRatio = sweep(Tmax_all,1,PET,"/")
  # m = apply(TmaxRatio, 2, median, na.rm=T)
  # lines(LAI_seq, m, col="black", lwd=2, lty=1)
  if(draw==TRUE) {
    TmaxPETGranier = -0.006*(LAI_seq^2)+0.134*LAI_seq
    plot(LAI_seq, TmaxPETGranier, type="l", col="gray", lwd=2, 
         xlab = "Stand's Leaf Area Index", ylab = "Tmax/PET", ylim=c(0,1))
    df2<-data.frame(LAI = LAI_seq)
    for(i in 1:ncoh) {
      m = mods[[i]]
      lines(df2$LAI, predict(m, newdata = df2), col=i, lwd=1, lty=i)
    }
    legend("topright", legend=c("Granier's equation"), lwd=2, lty=1, col=c("gray"),bty="n", cex =0.8)
    legend("topleft", legend=c(row.names(x$cohorts)), lwd=1, lty=1:ncoh, col=1:ncoh,bty="n", cex=0.8)
  }
  return(mods)
}
