multiple_runs<-function(parMatrix, x,
                        meteo, latitude,
                        elevation = NA, slope = NA, aspect = NA, 
                        summary_function = NULL, args = NULL, 
                        verbose = TRUE) {
  if(inherits(x, "spwbInput")) model = "spwb"
  else model = "growth"
  
  n = nrow(parMatrix)
  res = vector("list", n)
  parNames = colnames(parMatrix)
  for(r in 1:n) {
    x_r = .cloneInput(x)
    customParams = parMatrix[r,]
    names(customParams) <- parNames
    x_r = modifyInputParams(x_r, customParams, FALSE)
    x_r$control$verbose = FALSE
    S = do.call(model, list(x = x_r,
                            meteo = meteo, 
                            latitude = latitude, elevation = elevation,
                            slope  = slope,aspect = aspect))
    if(!is.null(summary_function)) {
      res[[r]] = do.call(summary_function, c(list(S), args))
      if(verbose) cat(paste0(r, ". Parameter values = [", paste0(customParams, collapse=", "), "] f = ", res[[r]], "\n"))
    } else {
      res[[r]] = S
      if(verbose) cat(paste0(r, ". Parameter values = [", paste0(customParams, collapse=", "), "]\n"))
    }
  }
  return(res)
}

optimization_function<-function(parNames, x, 
                                meteo, latitude,
                                elevation = NA, slope = NA, aspect = NA, 
                                summary_function, args= NULL) {
  
  if(inherits(x, "spwbInput")) model = "spwb"
  else model = "growth"
  
  x$control$verbose = FALSE
  yf<-function(v, verbose = FALSE) {
    x_i <- .cloneInput(x)
    if(is.vector(v)) {
      customParams = v
      names(customParams) <- parNames
      x_i = modifyInputParams(x_i, customParams, FALSE)
      S = do.call(model, list(x = x_i, 
                              meteo = meteo, 
                              latitude = latitude, elevation = elevation,
                              slope  = slope,aspect = aspect))
      y = do.call(summary_function, c(list(S), args))
      if(verbose) cat(paste0("Parameter values = [", paste0(customParams, collapse=", "), "] f = ", y, "\n"))
      return(y)
    } else if(is.matrix(v)) {
      colnames(v)<-parNames
      y<-multiple_runs(parMatrix = v, x = x_i, 
                       meteo = meteo,
                       latitude = latitude, elevation = elevation,
                       slope  = slope,aspect = aspect,
                       summary_function = summary_function, args = args,
                       verbose = verbose)
      return(as.numeric(y))
    } else {
      stop("Wrong 'v' class")
    }
  }
  return(yf)
}

optimization_evaluation_function<-function(parNames, x, 
                                           meteo, latitude,
                                           elevation = NA, slope = NA, aspect = NA, 
                                           measuredData, type = "SWC", cohorts = NULL, 
                                           temporalResolution = "day", 
                                           metric = "loglikelihood") {
  sf<-function(S) {
    if(!is.null(cohorts)) { # If cohorts != NULL evaluate output for each cohort and average the result
      y = rep(NA, length(cohorts))
      for(i in 1:length(cohorts)) {
        if(paste0(type,"_", cohorts[i]) %in% names(measuredData)) {
          y[i] = evaluation_metric(S, measuredData = measuredData, type=type, 
                                   cohort=cohorts[i], 
                                   temporalResolution = temporalResolution, metric = metric)
        }
      }
      return(mean(y, na.rm=TRUE))
    } else { # If cohorts == NULL evaluate output for first cohort (or no referred to any)
      y = evaluation_metric(S, measuredData = measuredData, type=type, 
                            cohort=NULL, 
                            temporalResolution = temporalResolution, metric = metric)
    }
    return(y)
  }
  return(optimization_function(parNames = parNames, x = x,
                               meteo = meteo, latitude = latitude,
                               elevation = elevation, slope = slope, aspect = aspect,
                               summary_function = sf))
}

optimization_multicohort_function<-function(cohortParNames, cohortNames, x, 
                                            meteo, latitude, 
                                            otherParNames = NULL,
                                            elevation = NA, slope = NA, aspect = NA, 
                                            summary_function, args= NULL) {
  
  if(inherits(x, "spwbInput")) model = "spwb"
  else model = "growth"
  
  x$control$verbose = FALSE
  yf<-function(v, verbose = FALSE) {
    x_i <- .cloneInput(x)
    if(is.vector(v)) {
      for(j in 1:length(cohortNames)) {
        customParams = v
        names(customParams) <- c(paste0(cohortNames[j],"/",cohortParNames), otherParNames)
        x_i = modifyInputParams(x_i, customParams, FALSE)
      }
      S = do.call(model, list(x = x_i, 
                              meteo = meteo, 
                              latitude = latitude, elevation = elevation,
                              slope  = slope,aspect = aspect))
      y = do.call(summary_function, c(list(S), args))
      if(verbose) cat(paste0("Parameter values = [", paste0(customParams, collapse=", "), "] f = ", y, "\n"))
      return(y)
    } else if(is.matrix(v)) {
      colnames(v) <- c(paste0(cohortNames[j],"/",cohortParNames), otherParNames)
      y<-multiple_runs(parMatrix = v, x = x_i, 
                       meteo = meteo,
                       latitude = latitude, elevation = elevation,
                       slope  = slope,aspect = aspect,
                       summary_function = summary_function, args = args,
                       verbose = verbose)
      return(as.numeric(y))
    } else {
      stop("Wrong 'v' class")
    }
  }
  return(yf)
}

optimization_evaluation_multicohort_function<-function(cohortParNames, cohortNames, x, 
                                                       meteo, latitude,
                                                       otherParNames = NULL,
                                                       elevation = NA, slope = NA, aspect = NA, 
                                                       measuredData, type = "SWC", cohorts = cohortNames,
                                                       temporalResolution = "day", 
                                                       metric = "loglikelihood") {
  sf<-function(S) {
    if(!is.null(cohorts)) { # If cohorts != NULL evaluate output for each cohort and average the result
      y = rep(NA, length(cohorts))
      for(i in 1:length(cohorts)) {
        if(paste0(type,"_", cohorts[i]) %in% names(measuredData)) {
          y[i] = evaluation_metric(S, measuredData = measuredData, type=type, 
                                   cohort=cohorts[i], 
                                   temporalResolution = temporalResolution, metric = metric)
        }
      }
      return(mean(y, na.rm=TRUE))
    } else { # If cohorts == NULL evaluate output for first cohort (or no referred to any)
      y = evaluation_metric(S, measuredData = measuredData, type=type, 
                            cohort=NULL, 
                            temporalResolution = temporalResolution, metric = metric)
    }
    return(y)
  }
  return(optimization_multicohort_function(cohortParNames = cohortParNames, 
                                           cohortNames = cohortNames, 
                                           x = x,
                                           meteo = meteo, latitude = latitude,
                                           otherParNames = otherParNames,
                                           elevation = elevation, slope = slope, aspect = aspect,
                                           summary_function = sf))
}

