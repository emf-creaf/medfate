#' Modify parameters
#' 
#' Routines to modify species parameter table or model input objects
#' 
#' @param x A model input object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
#' @param SpParams A species parameter data frame, typically \code{\link{SpParamsMED}}.
#' @param customParams A data frame or a named vector with new parameter values (see details).
#' @param subsetSpecies A logical flag to indicate that the output data frame should include only those species mentioned in \code{customParams}.
#' @param verbose A logical flag to indicate that messages should be printed on the console.
#' 
#' @details When calling function \code{modifySpParams}, \code{customParams} should be a data frame with as many rows as species 
#' and as many columns as parameters to modify, plus a column called 'SpIndex' to match species between the two tables.
#' 
#' When calling \code{modifyCohortParams}, \code{customParams} can be a data frame with as many rows as cohorts 
#' and as many columns as parameters to modify, plus a column called 'Cohort' which will be matched with the cohort names 
#' given by \code{\link{spwbInput}} or \code{\link{growthInput}}. 
#' Alternatively, \code{customParams} can be a named list or named numeric vector as for \code{modifyInputParams}.
#' 
#' When calling \code{modifyInputParams}, \code{customParams} must be either a named list or a named numeric vector. 
#' Cohort parameters are specified using the syntax "<cohortName>/<paramName>" for names (e.g. "T2_176/Z50" to modify parameter 'Z50' of cohort 'T2_176'). 
#' Soil layer parameters are specified using the syntax "<paramName>@#layer" for names, where #layer is the layer index (e.g. "rfc@1" will modify the rock fragment content of soil layer 1). 
#' Control parameters are specified using either "<paramName>" (e.g "phloemConductanceFactor") or "<paramName>$<subParamName>" (e.g "maximumRelativeGrowthRates$leaf"). 
#' It may seem unnecessary to modify soil or control parameters via a function, but \code{modifyInputParams} is called from optimization functions (see \code{\link{optimization}}).
#'  
#' @return Function \code{modifySpParams} returns a modified species parameter data frame. 
#' 
#' Functions \code{modifyCohortParams} and \code{modifyInputParams} return a modified \code{\link{spwbInput}} or \code{\link{growthInput}} object. Note that modifications may affect other parameters beyond those indicated in \code{customParams}, as a result of parameter dependencies (see examples).
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{spwbInput}}, \code{\link{SpParamsMED}}, \code{\link{optimization}}
#' 
#' @examples 
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
#' # Modify Z50 and Z95 of Pinus halepensis cohort 
#' customParams <- c(200,2000)
#' names(customParams) <- paste0(PH_coh,c("/Z50", "/Z95"))
#' x1m <- modifyInputParams(x1, customParams)
#' 
#' # Inspect original and modified objects 
#' x1$below
#' x1m$below
#' 
#' # Inspect dependencies: fine root distribution across soil layers
#' x1$belowLayers$V
#' x1m$belowLayers$V
#' 
#' # Modify rock fragment content and sand proportion of soil layer 1
#' x1s <- modifyInputParams(x1, c("rfc@1" = 5, "sand@1" = 10))
#' 
#' # Inspect original and modified soils 
#' x1$soil
#' x1s$soil
#' 
#' # When modifying growth input objects dependencies increase
#' x1 = forest2growthInput(exampleforestMED,examplesoil, SpParamsMED, control)
#' customParams <- c(2000,2)
#' names(customParams) <- paste0(PH_coh,c("/Al2As", "/LAI_live"))
#' x1m <- modifyInputParams(x1, customParams)
#' 
#' @name modifyParams

#' @rdname modifyParams
modifySpParams<-function(SpParams, customParams, subsetSpecies = TRUE) {
  
  # check if customParams exists, if not return SpParams without modification
  if (is.null(customParams)) {
    return(SpParams)
  }
  
  # get the names of the custom params and the SpParams
  custom_par <- names(customParams)
  if("SpIndex" %in% custom_par) custom_par = custom_par[-which(custom_par=="SpIndex")] # remove SpIndex from
  
  sp_par <- names(SpParams)
  
  # iterate between the custom params exisiting in SpParams
  for (param in custom_par) {
    
    # check if the param exists in SpParams
    if (param %in% sp_par) {
      
      # iterate by species, in case same variable has different values by sp
      for (sp in customParams[['SpIndex']]) {
        val <- customParams[which(customParams[['SpIndex']] == sp), param][1]
        if(!is.na(val)) {
          SpParams[which(SpParams[['SpIndex']] == sp), param] <- val
        }
      }
    }
  }
  # subset species
  if(subsetSpecies) {
    SpParams = SpParams[SpParams$SpIndex %in% customParams[['SpIndex']],]
  }
  # return the new SpParams
  return(SpParams)
}

#' @rdname modifyParams
modifyCohortParams<-function(x, customParams, verbose = TRUE) {
  
  # check if customParams exists, if not return x without modification
  if (is.null(customParams)) {
    return(x)
  }
  # check class of customParams
  if((!inherits(customParams, "data.frame")) && (!inherits(customParams, "numeric")) && (!inherits(customParams, "list"))) {
    stop("'customParams' must be a named numeric vector, a list or a data frame")
  }
  
  # get the names of the custom params and the input tables
  custom <- names(customParams)
  above_par <- names(x[['above']])
  below_par <- c("Z50","Z95")
  pheno_par <- names(x[['paramsPhenology']])
  base_par <- names(x[['paramsInterception']])
  transp_par <- names(x[['paramsTranspiration']])
  anatomy_par <- names(x[['paramsAnatomy']])
  waterstorage_par <- names(x[['paramsWaterStorage']])
  growth_par <- names(x[['paramsGrowth']])
  allom_par <- names(x[['paramsAllometries']])
  
  # clone object to modify
  x = .cloneInput(x)
  
  x_coh = row.names(x$cohorts)
  
  modifyParameterValue<-function(param, icoh, val) {
    if(!is.na(val)) {
      if(param %in% above_par) .modifyInputParam(x, "above", param, icoh - 1, val, verbose)
      if(param %in% below_par) .modifyInputParam(x, "below", param, icoh - 1, val, verbose) 
      if(param %in% base_par) .modifyInputParam(x, "paramsInterception", param, icoh - 1, val, verbose)
      if(param %in% transp_par) .modifyInputParam(x, "paramsTranspiration", param, icoh - 1, val, verbose) 
      if(param %in% anatomy_par)  .modifyInputParam(x, "paramsAnatomy", param, icoh - 1, val, verbose)
      if(!is.null(growth_par)) if (param %in% growth_par) .modifyInputParam(x, "paramsGrowth", param, icoh - 1, val, verbose)
      if(!is.null(allom_par)) if (param %in% allom_par) .modifyInputParam(x, "paramsAllometries", param, icoh - 1, val, verbose)
      if(!is.null(waterstorage_par)) if (param %in% waterstorage_par) .modifyInputParam(x, "paramsWaterStorage", param, icoh - 1, val, verbose)
    }
    return(x)
  }
  # iterate between the custom params
  if(inherits(customParams, "data.frame")) {
    if(sum(!(customParams[['Cohort']] %in% x_coh))>0) stop("Cohort names do not match between 'x' and 'customParams'")
    for (param in custom) {
      for (coh in customParams[['Cohort']]) {
        val <- customParams[customParams[['Cohort']] == coh, param]
        icoh <- which(x_coh==coh)
        modifyParameterValue(param, icoh, val)
      }
    }
  } 
  else { # assume a vector
    s <- strsplit(custom, "/")
    for(i in 1:length(s)) {
      paramCoh <- s[[i]]
      val <- customParams[[i]]
      coh <- paramCoh[[1]]
      if(!(coh %in% x_coh)) stop(paste0("Cohort '", coh,"' not found in 'x'"))
      icoh <- which(x_coh==coh)
      param <- paramCoh[[2]]
      # cat(paste0(coh, " ", param, " ", val,"\n"))
      modifyParameterValue(param, icoh, val)
    }
  }
  # return the modified input
  return(x)
}

#' @rdname modifyParams
modifyInputParams<-function(x, customParams, verbose = TRUE) {
  # check class of customParams
  if((!inherits(customParams, "numeric")) && (!inherits(customParams, "list"))) {
    stop("'customParams' must be a named numeric vector or a named list")
  }
  cn = names(customParams)
  isSoilParam = unlist(lapply(strsplit(cn, "@"), length))==2 #detect soil params
  isCohParam = unlist(lapply(strsplit(cn, "/"), length))==2 #detect cohort params
  customCohortParams = customParams[isCohParam]
  customSoilParams = customParams[isSoilParam]
  customControlParams = customParams[(!isCohParam) & (!isSoilParam)]
  # Modify control params
  if(length(customControlParams)>0) {
    for(i in 1:length(customControlParams)) {
      parName = names(customControlParams)[i]
      s = strsplit(parName, "\\$")[[1]]
      if(length(s)==1) {
        if(s %in% names(x[["control"]])) {
          x[["control"]][[s]] = customControlParams[[i]]
        } else {
          stop(paste0("Wrong control parameter name ", parName))
        }
      } else if(length(s)==2) {
        if((s[1] %in% names(x[["control"]])) && (s[2] %in% names(x[["control"]][[s[1]]]))) {
          x[["control"]][[s[1]]][[s[2]]] = customControlParams[[i]]
        } else {
          stop(paste0("Wrong control parameter name ", parName))
        }
      } else {
        stop(paste0("Wrong control parameter name ", parName))
      }
    }
  }
  # Modify cohort params
  if(length(customCohortParams)>0) {
    x = modifyCohortParams(x, customCohortParams, verbose)
  } else {
    x = .cloneInput(x)
  }
  # Modify soil layer params
  if(length(customSoilParams)>0) {
    s <- strsplit(names(customSoilParams), "@")
    for(i in 1:length(s)) {
      paramLayer <- s[[i]]
      val <- customSoilParams[[i]]
      param <- paramLayer[[1]]
      layer <- as.numeric(paramLayer[[2]])
      if(!(layer %in% 1:length(x$soil$dVec))) stop(paste0("Soil layer '", layer,"' not found in 'x'"))
      .modifySoilLayerParam(x$soil, param, layer-1, val)
    }
    .updateBelow(x)
  }
  return(x)
}