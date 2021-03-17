transp_maximumTranspirationModel<-function(x, meteo, latitude, elevation, slope, aspect,
                                           LAI_seq = c(0.1,0.25, seq(0.5, 10, by=0.5)),
                                           draw = TRUE) {
  
  
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
  # PET_cut = cut(PET,breaks = seq(0, max(PET), length.out = 20))
  # nPET = table(PET_cut)
  # PETw = as.numeric(1/nPET[PET_cut])
  # PETw = PETw/sum(PETw)
  # s = sample(1:nrow(meteo),ndays,replace = TRUE, prob = PETw)
  # s = s[order(PET[s])]
  # meteo = meteo[s,]
  # PET =PET[s]
  
  ncoh = nrow(x$above)
  ndays = nrow(meteo)
  nlai = length(LAI_seq)
  
  cohnames <- row.names(x$cohorts)
  LAItotal <- sum(x$above$LAI_live)
  
  xIni = x
  xIni$control$modifyInput = FALSE
  xIni$control$unlimitedSoilWater = TRUE
  xIni$control$cavitationRefill = "total"
  xIni$control$verbose = FALSE
  Tmax = matrix(NA, nrow=ndays, ncol = nlai)  
  colnames(Tmax) = LAI_seq
  rownames(Tmax) = row.names(meteo)
  LAI = matrix(NA, nrow=ndays, ncol = nlai)  
  colnames(LAI) = LAI_seq
  rownames(LAI) = row.names(meteo)
  
  s_res = vector("list", nlai)
  pb = txtProgressBar(0, nlai, style=3)
  for(j in 1:nlai) {
    customParams = LAI_seq[j]*(x$above$LAI_live/LAItotal)
    names(customParams) = paste0(cohnames,"/LAI_live")
    xlai = modifyInputParams(xIni, customParams, FALSE)
    s_res[[j]] = spwb(xlai, meteo,
                      latitude = latitude, 
                      elevation = elevation, slope = slope, aspect = aspect)
    Tmax[,j] = s_res[[j]]$WaterBalance$Transpiration
    LAI[,j] = s_res[[j]]$Stand$LAI
    setTxtProgressBar(pb, j)
  }
  TmaxRatio = sweep(Tmax,1,PET,"/")
  Tmaxratiovec = as.vector(TmaxRatio)
  laivec = as.vector(LAI)
  df = data.frame(y=Tmaxratiovec, LAI = laivec, Prec = meteo$Precipitation)
  df = df[df$Prec==0,] #Exclude precipitation days
  df = df[!is.na(df$y),, drop=FALSE] # Exclude missing ratio
  df = df[(df$y > 0.0) & (df$y < 1.0),, drop=FALSE] # Exclude extreme ratio
  mod <- glm(y ~ -1 + LAI + I(LAI^2), 
                   start = c(0.134,-0.006),
                   data =df, family=Gamma(link="identity"))
  if(draw==TRUE) {
    TmaxPETGranier = -0.006*(LAI_seq^2)+0.134*LAI_seq
    plot(LAI_seq, TmaxPETGranier, type="l", col="gray", lwd=2, 
         xlab = "Stand's Leaf Area Index", ylab = "Tmax/PET", ylim=c(0,1))
    df2<-data.frame(LAI = LAI_seq)
    lines(df2$LAI, predict(mod, newdata = df2), col="black", lwd=2)
    legend("topright", legend=c("Granier's equation", "Forest stand"),
           lwd=2, lty=1, col=c("gray", "black"),bty="n", cex =0.8)
  }
  return(mod)
}
