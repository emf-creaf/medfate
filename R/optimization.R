optimization_evaluation_function<-function(x, soil,
                                           cohNames, parNames, 
                                           measuredData, type = "SWC", cohort = NULL, SpParams = NULL, 
                                           metric = "loglikelihood",
                                           ...) {
  l = list(...)
  if(!("meteo" %in% names(l))) stop("Please specify 'meteo' in '...'")
  if(!("latitude" %in% names(l))) stop("Please specify 'latitude' in '...'")
  
  if(!("elevation" %in% names(l))) l[["elevation"]] = NA
  if(!("aspect" %in% names(l))) l[["aspect"]] = NA
  if(!("slope" %in% names(l))) l[["slope"]] = NA
  if(inherits(x, "spwbInput")) model = "spwb"
  else model = "growth"
  
  custom = data.frame(Cohort = unique(cohNames), as.data.frame(matrix(NA, nrow = length(unique(cohNames)), ncol = length(parNames))))
  names(custom) = c("Cohort", parNames)
  x$control$verbose = FALSE
  
  yf<-function(v) {
    resetInputs(x, soil)
    for(i in 1:length(parNames)) custom[custom$Cohort==cohNames[i], parNames[i]] = v[i]
    x = modifyCohortParams(x, custom, soil)
    # print(x$below)
    if(model=="spwb") {
      S = spwb(x, soil, 
               meteo = l[["meteo"]], 
               latitude = l[["latitude"]], elevation = l[["elevation"]],
               slope  = l[["slope"]],aspect = l[["aspect"]])
    } 
    else {
      S = growth(x, soil, 
                 meteo = l[["meteo"]], 
                 latitude = l[["latitude"]], elevation = l[["elevation"]],
                 slope  = l[["slope"]], aspect = l[["aspect"]])
    }
    y = evaluation_metric(S, measuredData = measuredData, type=type, 
                           cohort=cohort, SpParams = SpParams, metric = metric)
    return(y)
  }
  return(yf)
}

optimization_prediction_function<-function(x, soil,
                                           cohNames, parNames, 
                                           summary_function,
                                           ...) {
  l = list(...)
  if(!("meteo" %in% names(l))) stop("Please specify 'meteo' in '...'")
  if(!("latitude" %in% names(l))) stop("Please specify 'latitude' in '...'")
  
  if(!("elevation" %in% names(l))) l[["elevation"]] = NA
  if(!("aspect" %in% names(l))) l[["aspect"]] = NA
  if(!("slope" %in% names(l))) l[["slope"]] = NA
  if(inherits(x, "spwbInput")) model = "spwb"
  else model = "growth"
  
  custom = data.frame(Cohort = unique(cohNames), as.data.frame(matrix(NA, nrow = length(unique(cohNames)), ncol = length(parNames))))
  names(custom) = c("Cohort", parNames)
  x$control$verbose = FALSE
  
  yf<-function(v) {
    resetInputs(x, soil)
    for(i in 1:length(parNames)) custom[custom$Cohort==cohNames[i], parNames[i]] = v[i]
    x = modifyCohortParams(x, custom, soil)
    # print(x$below)
    if(model=="spwb") {
      S = spwb(x, soil, 
               meteo = l[["meteo"]], 
               latitude = l[["latitude"]], elevation = l[["elevation"]],
               slope  = l[["slope"]],aspect = l[["aspect"]])
    } 
    else {
      S = growth(x, soil, 
                 meteo = l[["meteo"]], 
                 latitude = l[["latitude"]], elevation = l[["elevation"]],
                 slope  = l[["slope"]], aspect = l[["aspect"]])
    }
    y = summary_function(S)
    # print(ll)
    return(y)
  }
  return(yf)
}