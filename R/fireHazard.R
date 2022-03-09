.firehaz_sim<-function(forest, x, SpParams, standardConditions = FALSE, freq = "days", fun = "max") {
  slope = x$topography[["slope"]]
  weather = x$weather
  fmc = x$Plants$LFMC
  if(is.na(slope)) slope = 0.0

  dates = as.Date(row.names(weather))
  ndays = length(dates)
  
  fb_vec = vector("list", ndays)
  for(i in 1:ndays){
    windSpeed = weather$WindSpeed[i]
    fccs = fuel_FCCS(forest, SpParams, cohortFMC = fmc[i,])
    vp = meteoland::utils_averageDailyVP(Tmin = weather$MinTemperature[i], Tmax = weather$MaxTemperature[i],
                                         RHmin = weather$MinRelativeHumidity[i], RHmax = weather$MaxRelativeHumidity[i])
    D = max(0, meteoland::utils_saturationVP(weather$MaxTemperature[i]) - vp)
    fm_dead = 5.43 + 52.91*exp(-0.64*D) # Resco de Dios, V., A. W. Fellows, R. H. Nolan, M. M. Boer, R. a. Bradstock, F. Domingo, and M. L. Goulden. 2015. A semi-mechanistic model for predicting the moisture content of fine litter. Agricultural and Forest Meteorology 203:64–73.
    MdeadSI = rep(fm_dead, 5)
    if(standardConditions) {
      windSpeed = 11.0
      fccs$ActFMC = NA
      MdeadSI = c(6,6,6,6,6)
    }
    fb_vec[[i]] = unlist(fire_FCCS(fccs, slope = slope, windSpeedSI = windSpeed,
                                   MdeadSI = MdeadSI))
  }
  fbm = matrix(NA, nrow = ndays, ncol = length(fb_vec[[1]]))
  rownames(fbm) <- row.names(weather)
  colnames(fbm) <- names(fb_vec[[1]])
  for(i in 1:ndays){
    fbm[i,] <- fb_vec[[i]]
  }
  if(freq!="days") {
    date.factor = cut(dates, breaks=freq)
    fbm = apply(fbm,2,function(x){tapply(x,INDEX=date.factor, FUN = fun, na.rm=TRUE)})
  } 
  return(fbm)
}
fireHazard<-function(x, SpParams, forest = NULL, standardConditions = FALSE,
                     freq="days", fun = "max") {
  if(!inherits(x, c("spwb","pwb","growth", "fordyn"))) {
    stop("'x' should be of class 'spwb', 'pwb', 'growth' or 'fordyn'")
  }
  if(inherits(x, c("spwb","pwb", "growth"))) {
    if(is.null(forest)) stop("You must supply a 'forest' object when 'x' is of class 'spwb', 'pwb' or 'growth'")
    res = .firehaz_sim(forest, x, SpParams, standardConditions,
                       freq = freq, fun = fun)
  } else {
    vec<-vector("list", length(x$GrowthResults))
    for(i in 1:length(x$GrowthResults)) {
      vec[[i]] <- .firehaz_sim(x$ForestStructures[[i]], x$GrowthResults[[i]], SpParams, standardConditions,
                               freq = freq, fun = fun)
    }
    res = .mergeVectorOfMatrices(vec)
  }
  return(res)
}