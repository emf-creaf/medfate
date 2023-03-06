.firehaz_sim<-function(forest, x, SpParams, standardConditions = FALSE, freq = "days", fun = "max") {
  slope = x$topography[["slope"]]
  weather = x$weather
  fmc = x$Plants$LFMC
  if(is.na(slope)) slope = 0.0

  dates = as.Date(row.names(weather))
  ndays = length(dates)
  
  #Calculate FCCS without FMC
  fccs = fuel_FCCS(forest, SpParams)
  cohLoading = plant_fuel(forest, SpParams)
  cohHeight = plant_height(forest, SpParams)
  cohHeight[is.na(cohHeight)] = 0
  cohCR = plant_crownRatio(forest,SpParams)
  cohCR[is.na(cohCR)] = 0
  
  fb_vec = vector("list", ndays)
  for(i in 1:ndays){
    windSpeed = weather$WindSpeed[i]
    
    # Estimate moisture of dead fine fuels (Resco de Dios et al. 2015)
    vp = meteoland::utils_averageDailyVP(Tmin = weather$MinTemperature[i], Tmax = weather$MaxTemperature[i],
                                         RHmin = weather$MinRelativeHumidity[i], RHmax = weather$MaxRelativeHumidity[i])
    D = max(0, meteoland::utils_saturationVP(weather$MaxTemperature[i]) - vp)
    fm_dead = 5.43 + 52.91*exp(-0.64*D) 

    #Calculate cohort canopy moisture to the average of canopy live and dead fuels, considering that a fraction of LAI is dead
    #proportionally to stem PLC (Ruffault et al. 2023)
    LFMC = x$Plants$LFMC[i,]
    PLC = x$Plants$StemPLC[i,]
    canopyFMC = (LFMC*(1.0 - PLC) + fm_dead*PLC)
    
    #Average canopy moisture in the crown and surface layers
    if(fccs$w[1]>0) fccs$ActFMC[1] = .layerFuelAverageParameter(200.0, 10000.0, canopyFMC, cohLoading, cohHeight, cohCR)
    else fccs$ActFMC[1] = NA
    if(fccs$w[2]>0.0) fccs$ActFMC[2] = .layerFuelAverageParameter(0.0, 200.0, canopyFMC, cohLoading, cohHeight, cohCR)
    else fccs$ActFMC[1] = NA

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

#' Fire hazard
#' 
#' Estimates potential fire behaviour at each daily step of a simulation
#' 
#' @param x An object of class \code{\link{spwb}}, \code{\link{spwb_day}}, \code{\link{pwb}}, \code{\link{growth}}, \code{\link{growth_day}} or \code{\link{fordyn}}.
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsDefinition}} and \code{\link{SpParamsMED}}).
#' @param forest An object of class \code{\link{forest}} (needed if \code{x} is not of class \code{\link{fordyn}}).
#' @param standardConditions A logical flag to indicate that standard fire weather conditions are to be used (instead of deriving fuel moisture and windspeed from \code{x}).
#' @param freq Frequency of summary statistics (see \code{\link{cut.Date}}).
#' @param fun Summary function (by default, maximum values).
#' 
#' @details Live fuel moisture of shrub and canopy layers is estimated from plant water status. 
#' Dead fuel moisture is estimated following Resco-de-Dios et al. (2015).
#' 
#' @return A matrix with fire behaviour variables (columns) for each simulated day (rows) or coarser time steps if summaries are requested.
#' 
#' @references 
#' 
#' Resco de Dios, V., A. W. Fellows, R. H. Nolan, M. M. Boer, R. A. Bradstock, F. Domingo, and M. L. Goulden. 2015. A semi-mechanistic model for predicting the moisture content of fine litter. Agricultural and Forest Meteorology 203:64–73.
#' 
#' Ruffault J, Limousin JM, Pimont F, Dupuy JL, De Cáceres M, Cochard H, Mouillot F, Blackman C, Torres-Ruiz JM, Parsons R, 
#' Moreno M, Delzon S, Jansen S, Olioso A, Choat B, Martin-StPaul N. 2023. Plant hydraulic modelling of leaf and canopy fuel moisture content reveals increasing vulnerability of a Mediterranean forest to wildfires under extreme drought. 
#' New Phytologist. (10.1111/nph.18614).
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{spwb}}, \code{\link{fuel_FCCS}}, \code{\link{fire_FCCS}}
#' 
#' @examples 
#' \donttest{
#' #Load example daily meteorological data
#' data(examplemeteo)
#' 
#' #Load example plot plant data
#' data(exampleforestMED)
#' 
#' #Default species parameterization
#' data(SpParamsMED)
#' 
#' #Initialize soil with default soil params (4 layers)
#' examplesoil <- soil(defaultSoilParams(4))
#' 
#' #Initialize control parameters
#' control <- defaultControl("Granier")
#' 
#' #Initialize input
#' x1 <- forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
#' 
#' #Call simulation function
#' S1 <- spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)
#' 
#' #Evaluate fire hazard
#' F1 <- fireHazard(S1, SpParamsMED, exampleforestMED)
#' }
fireHazard<-function(x, SpParams, forest = NULL, standardConditions = FALSE,
                     freq="days", fun = "max") {
  if(!inherits(x, c("spwb","pwb","growth", "fordyn", "spwb_day", "growth_day"))) {
    stop("'x' should be of class 'spwb', 'spwb_day', 'pwb', 'growth', 'growth_day' or 'fordyn'")
  }
  if(inherits(x, c("spwb","pwb", "growth"))) {
    if(is.null(forest)) stop("You must supply a 'forest' object when 'x' is of class 'spwb', 'spwb_day', 'pwb', 'growth' or 'growth_day'")
    res = .firehaz_sim(forest, x, SpParams, standardConditions,
                       freq = freq, fun = fun)
  } else if(inherits(x, c("spwb_day", "growth_day"))) {
    if(is.null(forest)) stop("You must supply a 'forest' object when 'x' is of class 'spwb', 'spwb_day', 'pwb', 'growth' or 'growth_day'")
    slope = x$topography[["slope"]]
    weather = x$weather
    
    # Estimate moisture of dead fine fuels (Resco de Dios et al. 2015)
    vp = meteoland::utils_averageDailyVP(Tmin = weather[["tmin"]], Tmax = weather[["tmax"]],
                                         RHmin = weather[["rhmin"]], RHmax = weather[["rhmax"]])
    D = max(0, meteoland::utils_saturationVP(weather[["tmax"]]) - vp)
    fm_dead = 5.43 + 52.91*exp(-0.64*D) 
    
    
    #Calculate cohort canopy moisture to the average of canopy live and dead fuels, considering that a fraction of LAI is dead
    #proportionally to stem PLC (Ruffault et al. 2023)
    canopyFMC= (x$Plants$LFMC*(1 - x$Plants$StemPLC) + fm_dead*x$Plants$StemPLC)
    
    if(is.na(slope)) slope = 0.0
    windSpeed = weather[["wind"]]
    fccs = fuel_FCCS(forest, SpParams, cohortFMC = canopyFMC)
    MdeadSI = rep(fm_dead, 5)
    if(standardConditions) {
      windSpeed = 11.0
      fccs$ActFMC = NA
      MdeadSI = c(6,6,6,6,6)
    }
    res = unlist(fire_FCCS(fccs, slope = slope, windSpeedSI = windSpeed, MdeadSI = MdeadSI))
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