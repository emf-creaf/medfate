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