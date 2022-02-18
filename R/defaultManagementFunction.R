defaultManagementFunction<-function(x, args, verbose = FALSE) {
  
  ntree = nrow(x$treeData)
  nshrub = nrow(x$shrubData)
  # Initialize output data
  N_tree_cut = rep(0, ntree)
  Cover_shrub_cut = rep(0, nshrub)
  planted_forest = emptyforest()
  
  Ntotal = sum(x$treeData$N, na.rm=TRUE)
  BAtotal = stand_basalArea(x)
  meanDBH = sum(x$treeData$N * x$treeData$DBH, na.rm=TRUE)/sum(x$treeData$N, na.rm=TRUE)
  
  action = "none"

  if(args$type == "irregular") action = "thinning"
  else if(args$type == "regular") {
    if(verbose) cat(paste0("  mean DBH: ", round(meanDBH,1), " threshold ", args$finalMeanDBH))
    if(meanDBH>args$finalMeanDBH && args$finalPreviousStage==0) { # If meanDBH exceeds threshold start final cuts
      action = "finalcut"
    } else if(args$finalPreviousStage>0) {
      if(args$finalYearsToCut==0) { # If we are in final cuts and this year needs cutting...
        action = "finalcut"
      } else { # If we are in final cuts but we have to wait, decrease year by one
        args$finalYearsToCut = args$finalYearsToCut - 1
      }
    } else {
      action = "thinning"
    }
  }
  planting = FALSE
  
  if(verbose) cat(paste0("  action: ", action))
  
  if(action=="thinning") {
    thin = FALSE
    
    enoughYearsForThinning = (args$yearsSinceThinning >= args$minThinningInterval)
    if(is.na(enoughYearsForThinning)) enoughYearsForThinning = TRUE
    
    if(args$thinningMetric=="BA") {
      if(verbose) cat(paste0("  Basal area: ", round(BAtotal,1), " threshold ", round(args$thinningThreshold)))
      if((BAtotal > args$thinningThreshold) && enoughYearsForThinning) thin = TRUE
    }
    else if(args$thinningMetric=="N") {
      if(verbose) cat(paste0("  Density: ", round(Ntotal,1), " threshold ", round(args$thinningThreshold)))
      if((Ntotal > args$thinningThreshold) && enoughYearsForThinning) thin = TRUE
    }
    else if(args$thinningMetric=="HB") {
      HB = stand_hartBeckingIndex(x)
      if(is.na(HB)) stop("NA Hart-becking index")
      if(verbose) cat(paste0("  Hart-Becking: ", round(HB,1), " threshold ", round(args$thinningThreshold)))
      if(HB < args$thinningThreshold && args$yearsSinceThinning >= args$minThinningInterval) thin = TRUE
    }
    else {
      stop(paste0("Non-recognized thinning metric '", args$thinningMetric,"'.\n"))
    }
    if(thin) {
      BA2remove = BAtotal*(args$thinningPerc/100)
      BAremoved = 0
      args$yearsSinceThinning = 1
      
      if(verbose) cat(paste0(", type: ",args$thinning,", BA to extract: ", round(BA2remove,1)))
      
      if(args$thinning %in% c("below", "above")) {
        o = order(x$treeData$DBH, decreasing = (args$thinning=="above"))
        cohort = 1
        while(BA2remove > 0) {
          r = (x$treeData$DBH[o[cohort]]/200)
          n = x$treeData$N[o[cohort]]
          BAcohort = pi*(r^2)*n
          if(BAcohort > BA2remove) {
            N_tree_cut[o[cohort]] = BA2remove/(pi*(r^2))
            BA2remove = 0
          } else {
            N_tree_cut[o[cohort]] = n
            BA2remove = BA2remove - BAcohort
          }
          cohort = cohort + 1
        }
      }
      else if (args$thinning%in% c("below-systematic", "above-systematic")) { #Cut half of target BA from below/above and the other half using trees from all diameters (keep current distribution)
        # Remove half as systematic
        sec = pi*(x$treeData$DBH/200)^2
        propN  = x$treeData$N/sum(x$treeData$N)
        HalfBA2remove = BA2remove/2
        while(HalfBA2remove > 0) {
          if(sum(sec*propN)> HalfBA2remove) { # Correct to avoid cutting above existences
            propN = propN*(HalfBA2remove/sum(sec*propN))
          }
          Nremoved = Nremoved + propN
          HalfBA2remove = HalfBA2remove - sum(sec*propN)
        }
        BA2remove = BA2remove/2
        #Remove remaining as below/above
        o = order(x$treeData$DBH, decreasing = (args$thinning=="above-systematic"))
        cohort = 1
        while(BA2remove > 0) {
          r = (x$treeData$DBH[o[cohort]]/200)
          n = x$treeData$N[o[cohort]] - N_tree_cut[o[cohort]]
          BAcohort = pi*(r^2)*n
          if(BAcohort > BA2remove) {
            N_tree_cut[o[cohort]] =  N_tree_cut[o[cohort]] + BA2remove/(pi*(r^2))
            BA2remove = 0
          } else {
            N_tree_cut[o[cohort]] = N_tree_cut[o[cohort]]+n
            BA2remove = BA2remove - BAcohort
          }
          cohort = cohort + 1
        }
      }
      else if (args$thinning=="systematic") { #Cut trees from all diameters (keep current distribution)
        sec = pi*(x$treeData$DBH/200)^2
        propN  = x$treeData$N/sum(x$treeData$N)
        while(BA2remove > 0) {
          if(sum(sec*propN)> BA2remove) { # Correct to avoid cutting above existences
            propN = propN*(BA2remove/sum(sec*propN))
          }
          N_tree_cut = N_tree_cut + propN
          BA2remove = BA2remove - sum(sec*propN)
        }
      }
      else {
        s = strsplit(args$thinning,"/")[[1]]
        # print(s)
        ncl = length(s)
        breaks = rep(0, ncl)
        props = rep(0, ncl)
        for(i in 1:ncl) {
          s2  = strsplit(s[i],"-")[[1]]
          breaks[i] = as.numeric(s2[1])
          props[i] = as.numeric(s2[2])
        }
        breaks[ncl] = max(breaks[ncl], max(x$treeData$DBH)+1) # To include largest trees
        breaks = c(0, breaks)
        # print(breaks)
        props = props/sum(props)
        cl = as.numeric(cut(x$treeData$DBH, breaks))
        propNmat = matrix(0, ntree, ncol = ncl)
        sec = pi*(x$treeData$DBH/200)^2
        BAcl = rep(0, ncl)
        for(i in 1:ncl) {
          # print(sum(cl==i))
          propNmat[cl==i,i] = x$treeData$N[cl==i]/sum(x$treeData$N[cl==i])
          BAcl[i] = sum(sec[cl==i]*x$treeData$N[cl==i])
        }
        # print(propNmat)
        BA2removecl = BA2remove*props
        # print(cbind(BA2removecl, BAcl))
        BA2removecl = pmin(BA2removecl, BAcl) # Do not allow to cut more than is available for a given class
        for(i in 1:ncl) {
          propN = propNmat[,i]
          while(BA2removecl[i] > 0) {
            if(sum(sec*propN)> BA2removecl[i]) { # Correct to avoid cutting above existences
              propN = propN*(BA2removecl[i]/sum(sec*propN))
            }
            N_tree_cut = N_tree_cut + propN
            BA2removecl[i] = BA2removecl[i] - sum(sec*propN)
          }
        }
      }
    }
    else {
      action = "none"
      args$yearsSinceThinning = args$yearsSinceThinning + 1
      
    }
  }
  else if(action=="finalcut") {
    #Proportions in different stages
    if(verbose) cat(paste0(", previous stage: ",args$finalPreviousStage))
    args$finalPreviousStage = args$finalPreviousStage + 1
    args$finalYearsToCut = args$finalYearsBetweenCuts
    props = as.numeric(strsplit(args$finalPerc,"-")[[1]])/100
    BA2remove = BAtotal*props[args$finalPreviousStage]
    if(verbose) cat(paste0(", % to extract: ", round(props[args$finalPreviousStage]*100),", BA to extract: ", round(BA2remove,1)))
    BAremoved = 0
    #Cut from above in final cuts
    o = order(x$treeData$DBH, decreasing = TRUE)
    cohort = 1
    while((BA2remove > 0) && (cohort<=length(o))) {
      # print(BA2remove)
      d = x$treeData$DBH[o[cohort]]
      r = (d/200)
      n = x$treeData$N[o[cohort]]
      if((d > 12.5) || (!is.na(args$plantingSpecies))) { # Do not cut regeneration if there is no planting
        BAcohort = pi*(r^2)*n
        if(BAcohort > BA2remove) {
          N_tree_cut[o[cohort]] = BA2remove/(pi*(r^2))
          BA2remove = 0
        } else {
          N_tree_cut[o[cohort]] = n
          BA2remove = BA2remove - BAcohort
        }
      }
      cohort = cohort + 1
    }
    if(args$finalPreviousStage==length(props)) {
      args$finalPreviousStage = 0
      if(!is.na(args$plantingSpecies)) planting = TRUE
    }
    if(verbose) cat(paste0(", final stage: ",args$finalPreviousStage))
  }
  # If tree felling occurred, apply understory clearing
  if(sum(N_tree_cut)>0) {
    stand_shrub_cover = sum(x$shrubData$Cover, na.rm=TRUE)
    stand_to_rem = pmax(0,stand_shrub_cover - args$understoryMaximumCover)
    Cover_shrub_cut = x$shrubData$Cover*(stand_to_rem/stand_shrub_cover)
  }
  # Apply planting if needed (regular management)
  if(planting) {
    action = paste0(action, " + planting")
    if(verbose) cat(paste0(", planting"))
    planted_forest = emptyforest(ntree=1)
    planted_forest$treeData$Species[1] = args$plantingSpecies
    planted_forest$treeData$DBH[1] = args$plantingDBH
    planted_forest$treeData$Height[1] = args$plantingHeight
    planted_forest$treeData$N[1] = args$plantingDensity
  }
  # Return
  return(list(action = action,
              N_tree_cut = N_tree_cut,
              Cover_shrub_cut = Cover_shrub_cut,
              planted_forest = planted_forest, 
              management_args = args))
}
defaultManagementArguments<-function(){
  return(list(
    type = "irregular",
    thinning = "below", 
    thinningMetric = "BA", 
    thinningThreshold = 20,
    thinningPerc = 30,
    minThinningInterval = 10,
    yearsSinceThinning = NA,
    finalMeanDBH = 20, 
    finalPerc = "40-60-100",
    finalPreviousStage = 0,
    finalYearsBetweenCuts = 10,
    finalYearsToCut = NA,
    plantingSpecies = NA,  
    plantingDBH = 5, 
    plantingHeight = 5, 
    plantingDensity = 200,
    understoryMaximumCover = 3
  ))
}