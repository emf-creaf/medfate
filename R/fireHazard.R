.firehaz_sim<-function(forest, x, SpParams, standardConditions = FALSE) {
  slope = x$topography[["slope"]]
  weather = x$weather
  fmc = moisture_cohortFMC(x, SpParams)
  if(is.na(slope)) slope = 0.0

  ndays = nrow(weather)
  fb_vec = vector("list", ndays)
  for(i in 1:ndays){
    windSpeed = weather$WindSpeed[i]
    fccs = fuel_FCCS(forest, SpParams, cohortFMC = fmc[i,])
    if(standardConditions) {
      windSpeed = 11.0
      fccs$ActFMC = NA
    }
    fb_vec[[i]] = unlist(fire_FCCS(fccs, slope = slope, windSpeedSI = windSpeed))
  }
  fbm = matrix(NA, nrow = ndays, ncol = length(fb_vec[[1]]))
  rownames(fbm) <- row.names(weather)
  colnames(fbm) <- names(fb_vec[[1]])
  for(i in 1:ndays){
    fbm[i,] <- fb_vec[[i]]
  }
  return(fbm)
}
fireHazard<-function(x, SpParams, forest = NULL, standardConditions = FALSE) {
  if(!inherits(x, c("spwb","pwb","growth", "fordyn"))) {
    stop("'x' should be of class 'spwb', 'pwb', 'growth' or 'fordyn'")
  }
  if(inherits(x, c("spwb","pwb", "growth"))) {
    if(is.null(forest)) stop("You must supply a 'forest' object when 'x' is of class 'spwb', 'pwb' or 'growth'")
    res = .firehaz_sim(forest, x, SpParams, standardConditions)
  } else {
    vec<-vector("list", length(x$GrowthResults))
    for(i in 1:length(x$GrowthResults)) {
      vec[[i]] <- .firehaz_sim(x$ForestStructures[[i]], x$GrowthResults[[i]], SpParams, standardConditions)
    }
    res = .mergeVectorOfMatrices(vec)
  }
  return(res)
}