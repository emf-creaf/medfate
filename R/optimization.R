#' Multiple model runs and function factories for optimization routines
#'
#' Function factories to generate functions to be used in model calibration, uncertainty or sensitivity analysis.
#'
#' @param parMatrix A matrix of parameter values with runs in rows and parameters in columns. Column names should follow parameter modification naming rules (see examples and naming rules in \code{\link{modifyInputParams}}).
#' @param parNames A string vector of parameter names (see examples and naming rules in \code{\link{modifyInputParams}}).
#' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
#' @param meteo,latitude,elevation,slope,aspect Additional parameters to simulation functions \code{\link{spwb}} or \code{\link{growth}}.
#' @param verbose A flag to indicate extra console output.
#' @param measuredData A data frame with observed/measured values. Dates should be in row names, whereas columns should be named according to the type of output to be evaluated (see details).
#' @param type A string with the kind of model output to be evaluated. Accepted values are \code{"SWC"} (soil moisture content), \code{"REW"} relative extractable water, \code{"ETR"} (total evapotranspiration), \code{"E"} (transpiration per leaf area), \code{"LFMC"} (live fuel moisture content) and \code{"WP"} (plant water potentials).
#' @param cohorts A string or a vector of strings with the cohorts to be compared (e.g. "T1_68"). If several cohort names are provided, the function \code{optimization_cohorts_function} evaluates the performance for each one and provides the mean value. If \code{NULL} results for the first cohort will be evaluated.
#' @param temporalResolution A string to indicate the temporal resolution of the model evaluation, which can be "day", "week", "month" or "year". Observed and modelled values are aggregated temporally (using either means or sums) before comparison.
#' @param metric An evaluation metric (see \code{\link{evaluation_metric}}).
#' @param summary_function A function whose input is the result of \code{\link{spwb}} or \code{\link{growth}}. The function must return a numeric scalar in the case of \code{optimization_function}, but is not restricted in the case of \code{multiple_runs}.
#' @param args A list of additional arguments of \code{optimization_function}.
#' @param cohortParNames A string vector of vegetation parameter names for cohorts (e.g. 'Z95' or 'psiExtract').
#' @param cohortNames A string vector of cohort names. All cohorts will be given the same parameter values for each parameter in 'cohortParNames'.
#' @param otherParNames A string vector of parameter names (see examples and naming rules in \code{\link{modifyInputParams}}) for non-vegetation parameters (i.e. control parameters and soil parameters).
#' 
#' @details 
#' See \code{\link{evaluation}} for details regarding how to specify measured data.
#' 
#' Functions produced by these function factories should be useful for sensitivity analyses using package 'sensitivity'.
#' 
#' Parameter naming (i.e. \code{parNames}) should follow the rules specified in section details of \code{\link{modifyInputParams}}. 
#' The exception to the naming rules applies when multiple cohorts are to be modified to the same values with functions 
#' \code{optimization_multicohort_function} and \code{optimization_evaluation_multicohort_function}. 
#' Then, only a vector of parameter names is supplied for \code{cohortParNames}.
#' 
#' @return 
#' Function \code{multiple_runs} returns a list, whose elements are either the result of calling simulation models 
#' or the result of calling \code{summary_function} afterwards. 
#' 
#' Function \code{optimization_function} returns a function whose parameters are parameter values 
#' and whose return is a prediction scalar (e.g. total transpiration).
#' 
#' Function \code{optimization_evaluation_function} returns a function whose parameters are parameter values 
#' and whose return is an evaluation metric (e.g. loglikelihood of the data observations given model predictions). 
#' If evaluation data contains information for different cohorts (e.g. plant water potentials or transpiration rates) 
#' then the evaluation is performed for each cohort and the metrics are averaged. 
#' 
#' Function \code{optimization_multicohorts_function} returns a function whose parameters are parameter values 
#' and whose return is a prediction scalar (e.g. total transpiration). The difference with \code{optimization_function} 
#' is that multiple cohorts are set to the same parameter values.
#' 
#' Function \code{optimization_evaluation_multicohort_function} returns a function whose parameters are parameter values 
#' and whose return is an evaluation metric (e.g. loglikelihood of the data observations given model predictions). 
#' If evaluation data contains information for different cohorts (e.g. plant water potentials or transpiration rates) 
#' then the evaluation is performed for each cohort and the metrics are averaged. The difference with \code{optimization_evaluation_function} 
#' is that multiple cohorts are set to the same parameter values.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{evaluation_metric}}, \code{\link{modifyInputParams}}, \code{\link{spwb}}, \code{\link{growth}}
#' 
#' @examples 
#' \dontrun{
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
#' examplesoil = soil(defaultSoilParams(4))
#' 
#' #Initialize control parameters
#' control = defaultControl("Granier")
#' 
#' #Initialize input
#' x1 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
#' 
#' # Cohort name for Pinus halepensis
#' PH_coh = paste0("T1_", SpParamsMED$SpIndex[SpParamsMED$Name=="Pinus halepensis"])
#' PH_coh 
#' 
#' #Parameter names of interest
#' parNames = c(paste0(PH_coh,"/Z50"), paste0(PH_coh,"/Z95"))
#' 
#' #Specify parameter matrix
#' parMatrix <- cbind(c(200,300), c(500,1000))
#' colnames(parMatrix) <- parNames
#' 
#' #Define a summary function as the total transpiration over the simulated period
#' sf<-function(x) {sum(x$WaterBalance$Transpiration, na.rm=TRUE)}
#' 
#' #Perform two runs and evaluate the summary function
#' multiple_runs(parMatrix, 
#'               x1, examplemeteo, latitude = 42, elevation = 100,
#'               summary_function = sf)
#' 
#' #Load observed data (in this case the same simulation results with some added error)  
#' # Generate a prediction function for total transpiration over the simulated period
#' # as a function of parameters "Z50" and "Z95" for Pinus halepensis cohort 
#' of<-optimization_function(parNames = parNames,
#'                           x = x1,
#'                           meteo = examplemeteo, 
#'                           latitude = 41.82592, elevation = 100,
#'                           summary_function = sf)
#' 
#' # Evaluate for the values of the parameter matrix
#' of(parMatrix[1, ])
#' of(parMatrix)
#' 
#' 
#' # Generate a loglikelihood function for soil water content
#' # as a function of parameters "Z50" and "Z95" for Pinus halepensis cohort 
#' data(exampleobs)
#' oef<-optimization_evaluation_function(parNames = parNames,
#'                                       x = x1,
#'                                       meteo = examplemeteo, latitude = 41.82592, elevation = 100,
#'                                       measuredData = exampleobs, type = "SWC", 
#'                                       metric = "loglikelihood")
#' 
#' # Loglikelihood for the values of the parameter matrix
#' oef(parMatrix[1, ])
#' oef(parMatrix)
#' }
#' 
#' @name optimization
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

#' @rdname optimization
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

#' @rdname optimization
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

#' @rdname optimization
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

#' @rdname optimization
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

