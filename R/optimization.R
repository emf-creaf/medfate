multiple_runs<-function(parMatrix, x, soil,
                        meteo, latitude,
                        elevation = NA, slope = NA, aspect = NA, 
                        summary_function = NULL, args = NULL) {
  if(inherits(x, "spwbInput")) model = "spwb"
  else model = "growth"
  
  n = nrow(parMatrix)
  res = vector("list", n)
  parNames = colnames(parMatrix)
  for(r in 1:n) {
    resetInputs(x, soil)
    customParams = parMatrix[r,]
    names(customParams) <- parNames
    x = modifyInputParams(x, customParams, soil)
    x$control$verbose = FALSE
    S = do.call(model, list(x = x, soil = soil, 
                            meteo = meteo, 
                            latitude = latitude, elevation = elevation,
                            slope  = slope,aspect = aspect))
    if(!is.null(summary_function)) {
      res[[r]] = do.call(summary_function, c(list(S), args))
    } else {
      res[[r]] = S
    }
  }
  return(res)
}

optimization_function<-function(parNames, x, soil,  
                                meteo, latitude,
                                elevation = NA, slope = NA, aspect = NA, 
                                summary_function, args= NULL) {
  
  if(inherits(x, "spwbInput")) model = "spwb"
  else model = "growth"
  
  x$control$verbose = FALSE
  
  yf<-function(v) {
    resetInputs(x, soil)
    if(is.vector(v)) {
      customParams = v
      names(customParams) <- parNames
      x = modifyInputParams(x, customParams, soil)
      S = do.call(model, list(x = x, soil = soil, 
                              meteo = meteo, 
                              latitude = latitude, elevation = elevation,
                              slope  = slope,aspect = aspect))
      y = do.call(summary_function, c(list(S), args))
      return(y)
    } else if(is.matrix(v)) {
      colnames(v)<-parNames
      y<-multiple_runs(parMatrix = v, x = x, soil = soil, 
                       meteo = meteo,
                       latitude = latitude, elevation = elevation,
                       slope  = slope,aspect = aspect,
                       summary_function = summary_function, args = args)
      return(as.numeric(y))
    } else {
      stop("Wrong 'v' class")
    }
  }
  return(yf)
}

optimization_evaluation_function<-function(parNames, x, soil, 
                                           meteo, latitude,
                                           elevation = NA, slope = NA, aspect = NA, 
                                           measuredData, type = "SWC", cohort = NULL, SpParams = NULL, 
                                           metric = "loglikelihood") {
  sf<-function(S) {
    y = evaluation_metric(S, measuredData = measuredData, type=type, 
                           cohort=cohort, SpParams = SpParams, metric = metric)
    return(y)
  }
  return(optimization_function(parNames = parNames, x = x, soil = soil,
                               meteo = meteo, latitude = latitude,
                               elevation = elevation, slope = slope, aspect = aspect,
                               summary_function = sf))
}

