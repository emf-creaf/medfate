defaultManagementFunction<-function(x, args, verbose = FALSE) {
  
  # Next year arguments
  nextArgs = list(type = args$type, 
                  thinning = args$thinning, 
                  thinningBA = args$thinningBA, 
                  thinningHB = args$thinningHB, 
                  thinningPerc = args$thinningPerc,
                  finalMeanDBH = args$finalMeanDBH, 
                  finalPerc = args$finalPerc,
                  plantingSpecies = args$plantingSpecies, 
                  plantingDBH = args$plantingDBH, 
                  plantingHeight = args$plantingHeight, 
                  plantingDensity = args$plantingDensity,
                  finalPreviousStage = args$finalPreviousStage,
                  finalYearsToCut = args$finalYearsToCut)
  ntree = nrow(x$treeData)
  nshrub = nrow(x$shrubData)
  # Initialize output data
  N_tree_cut = rep(0, ntree)
  Cover_shrub_cut = rep(0, nshrub)
  planted_forest = emptyforest()
  
  BAtotal = stand_basalArea(x)
  meanDBH = sum(x$treeData$N * x$treeData$DBH, na.rm=TRUE)/sum(x$treeData$N, na.rm=TRUE)
  
  if(type == "irregular") action = "thinning"
  else if(type == "regular") {
    if(verbose) cat(paste0("  mean DBH: ", round(meanDBH,1), " threshold ", args$finalMeanDBH))
    if((meanDBH>args$finalMeanDBH && args$finalPreviousStage==0) || ((args$finalPreviousStage>0) && (args$finalYearsToCut==0))) {
      action = "finalcut"
    } else {
      action = "thinning"
    }
  }
  planting = FALSE
  
  if(verbose) cat(paste0("  action: ", action))
  
  if(action=="thinning") {
    thin = FALSE
    if(!is.na(args$thinningBA)) {
      if(verbose) cat(paste0("  BA: ", round(BAtotal,1), " threshold ", round(args$thinningBA)))
      if(BAtotal > args$thinningBA) thin = TRUE
    }
    else if(!is.na(thinningHB)) {
    }
    
  }
  # Return
  return(list(N_tree_cut = N_tree_cut,
              Cover_shrub_cut = Cover_shrub_cut,
              planted_forest = planted_forest, 
              management_args = nextArgs))
}