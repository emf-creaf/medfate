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

modifyCohortParams<-function(x, customParams) {
  
  # check if customParams exists, if not return swbInput without modification
  if (is.null(customParams)) {
    return(x)
  }
  
  # get the names of the custom params and the input tables
  custom <- names(customParams)
  above_par <- names(x[['above']])
  base_par <- names(x[['paramsBase']])
  transp_par <- names(x[['paramsTransp']])
  anatomy_par <- names(x[['paramsAnatomy']])
  waterstorage_par <- names(x[['paramsWaterStorage']])
  
  # iterate between the custom params
  for (param in custom) {
    for (coh in customParams[['Cohort']]) {
      val <- customParams[customParams[['Cohort']] == coh, param]
      if(!is.na(val)) {
        if (param %in% above_par) x[['above']][[coh, param]] <- val
        if (param %in% base_par) x[['paramsBase']][[coh, param]] <- val
        if (param %in% transp_par) x[['paramsTransp']][[coh, param]] <- val
        if (param %in% anatomy_par) x[['paramsAnatomy']][[coh, param]] <- val
        if (param %in% waterstorage_par) x[['paramsWaterStorage']][[coh, param]] <- val
      }
    }
  }

  # return the modified input
  return(x)
}