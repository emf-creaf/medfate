#' Default forest management actions
#'
#' Function \code{defaultManagementFunction} implements actions for 'regular' and 'irregular' management models of monospecific or mixed stands, 
#' whereas function \code{defaultManagementArguments} returns a list with default values for the parameters regulating management. 
#' Both functions are meant to be used in simulations with \code{\link{fordyn}}.
#' 
#' @param x An object of class \code{\link{forest}}
#' @param args A list of arguments regulating management actions, e.g. the list returned by \code{defaultManagementArguments}
#' @param verbose A logical flag enabling console printing
#' 
#' @details This function implements silvicultural actions following either 'regular' or 'irregular' management models. 
#' Irregular models are implemented by executing thinning operations only, whereas regular models include both thinning and a set of final cuts. 
#' Thinning occurs anytime a stand-level metric (e.g. basal area) crosses a given threshold, and different kinds of thinning operations are allowed. 
#' Unrealistic high frequency thinning can be avoided by setting a minimum number of years to happen between thinning operations. 
#' Final cuts start whenever mean DBH exceeds a given threshold, and may include different cuts separated a number of years. 
#' The function can be applied to target management of specific taxa (instead of assuming a monospecific stand), but the thresholds that determine
#' thinning operations apply to stand-level metrics. Mean DBH will be calculated for the target species only.
#' Planting is only allowed under regular management models, and is applied after the last final cut. 
#' Understory clearings are assumed to occur anytime there is an intervention on trees, an only a residual shrub cover is left.
#' 
#' \emph{Thinning types}:
#'   \itemize{
#'     \item{\code{above}: Extract largest trees (according to DBH) until thinning objective is met.}
#'     \item{\code{below}: Extract smallest trees (according to DBH) until thinning objective is met.}
#'     \item{\code{systematic}: Extract equally from all size classes until thinning objective is met.}
#'     \item{\code{above-systematic}: Extract half the objective as systematic thinning and the other hald as above thinning.}
#'     \item{\code{below-systematic}: Extract half the objective as systematic thinning and the other hald as below thinning.}
#'     \item{\code{free string}: A string specifying the proportion of tree cuts from size classes, with size classes separated by "/" and each one composed of a number specifying the upper limit and a number indicating its proportion, separated by "-" (e.g. "10-50/40-30/60-20").}
#'   }
#' 
#' @return 
#' Function \code{defaultManagementFunction} returns a list with the following items:
#'   \itemize{
#'     \item{\code{"action"}: A string identifying the action performed (e.g. "thinning").}
#'     \item{\code{"N_tree_cut"}: A vector with the density of trees removed.}
#'     \item{\code{"Cover_shrub_cut"}: A vector with the cover of shrubs removed.} 
#'     \item{\code{"planted_forest"}: An object of class \code{\link{forest}} with the new plant cohorts resulting from tree/shrub planting.}
#'     \item{\code{"management_args"}: A list of management arguments to be used in the next call to the management function.}
#'   }
#' 
#' Function \code{defaultManagementArguments} returns a list with default arguments:
#'  \itemize{
#'    \item{\code{"type"}: Management model, either 'regular' or 'irregular'.}
#'    \item{\code{"targetTreeSpecies"}: Either \code{"all"} for unspecific cuttings or a numeric vector of target tree species to be selected for cutting operations.}
#'    \item{\code{"thinning"}: Kind of thinning to be applied in irregular models or in regular models before the final cuts. Options are 'below', 'above', 'systematic', 'below-systematic', 'above-systematic' or a string with the proportion of cuts to be applied to different diameter sizes (see details).}
#'    \item{\code{"thinningMetric"}: The stand-level metric used to decide whether thinning is applied, either 'BA' (basal area), 'N' (density) or 'HB' (Hart-Becking index).} 
#'    \item{\code{"thinningThreshold"}: The threshold value of the stand-level metric causing the thinning decision.}
#'    \item{\code{"thinningPerc"}: Percentage of stand's basal area to be removed in thinning operations.}
#'    \item{\code{"minThinningInterval"}: Minimum number of years between thinning operations.}
#'    \item{\code{"yearsSinceThinning"}: State variable to count the years since the last thinning ocurred.}
#'    \item{\code{"finalMeanDBH"}: Mean DBH threshold to start final cuts.}
#'    \item{\code{"finalPerc"}: String with percentages of basal area to be removed in final cuts, separated by '-' (e.g. "40-60-100").}
#'    \item{\code{"finalPreviousStage"}: Integer state variable to store the stage of final cuts ('0' before starting final cuts).}
#'    \item{\code{"finalYearsBetweenCuts"}: Number of years separating final cuts.}
#'    \item{\code{"finalYearsToCut"}: State variable to count the years to be passed before new final cut is applied.}
#'    \item{\code{"plantingSpecies"}: Species code to be planted. If missing, planting does not occur and only natural regeneration is allowed.}
#'    \item{\code{"plantingDBH"}: Initial DBH (cm) of planted species.}
#'    \item{\code{"plantingHeight"}: Initial height (cm) of planted species.}
#'    \item{\code{"plantingDensity"}: Initial density (ind./ha) of the planted species.}
#'    \item{\code{"understoryMaximumCover"}: Percentage of overall shrub cover to be left after any silvicultural intervention. 
#'    If missing, shrub cover will not be left unmodified.}
#'  }
#' 
#' @author 
#' Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' Aitor \enc{Améztegui}{Caceres}, UdL
#' 
#' Jose-Ramon Gonzalez Olabarria, CTFC
#' 
#' @seealso \code{\link{fordyn}}
#' 
#' @examples
#' # Load example forest object
#' data(exampleforestMED)
#'   
#' # Define arguments
#' args = defaultManagementArguments()
#'   
#' # Call management function
#' f = defaultManagementFunction(exampleforestMED, args)
#'   
#' #list names
#' names(f)
#'   
#' # Action performed
#' f$action
#'   
#' # Number of trees cut for each cohort
#' f$N_tree_cut
#'   
#' # Percent cover of shrubs removed
#' f$Cover_shrub_cut
#'   
#' @name defaultManagementFunction
defaultManagementFunction<-function(x, args, verbose = FALSE) {
  
  if(is.null(args)) stop("Please supply a list of management arguments")
  
  # Set missing final previous stage to 0
  if(is.na(args$finalPreviousStage)) args$finalPreviousStage = 0
  # If missing years since thinning, allow thinning by setting it equal to the minimum years for thinning
  if(is.na(args$yearsSinceThinning)) args$yearsSinceThinning = args$minThinningInterval
  
  ntree = nrow(x$treeData)
  nshrub = nrow(x$shrubData)
  # Initialize output data
  N_tree_cut = rep(0, ntree)
  Cover_shrub_cut = rep(0, nshrub)
  planted_forest = emptyforest()
  
  Ntotal = sum(x$treeData$N, na.rm=TRUE)
  BAtotal = stand_basalArea(x)
  
  if(length(args$targetTreeSpecies)==0) args$targetTreeSpecies = "all"
    
  if(args$targetTreeSpecies[1] == "all") {
    isTarget = rep(TRUE, ntree)
  } else {
    isTarget = (x$treeData$Species %in% args$targetTreeSpecies)
  }
  
  # Calculate mean DBH (of target species if specified)
  meanDBH = sum(x$treeData$N[isTarget] * x$treeData$DBH[isTarget], na.rm=TRUE)/sum(x$treeData$N[isTarget], na.rm=TRUE)
  
  action <- "none"

  if(args$type == "irregular") {
    action <- "thinning"
  } else if(args$type == "regular") {
    if(!is.na(meanDBH)) {
      if(is.na(args$finalMeanDBH)) stop("Argument 'finalMeanDBH' is missing in a regular management model.")
      if(verbose) cat(paste0("  mean DBH: ", round(meanDBH,1), " threshold ", args$finalMeanDBH))
      if(meanDBH > args$finalMeanDBH && args$finalPreviousStage==0) { # If meanDBH exceeds threshold start final cuts
        action <- "finalcut"
      } else if(args$finalPreviousStage>0) {
        if(args$finalYearsToCut==0) { # If we are in final cuts and this year needs cutting...
          action <- "finalcut"
        } else { # If we are in final cuts but we have to wait, decrease year by one
          args$finalYearsToCut <- args$finalYearsToCut - 1
        }
      } else {
        action <- "thinning"
      }
    } else {
      warning("Mean DBH for the target species is missing. Cannot apply management.")
    }
  }
  planting = FALSE
  
  if(verbose) cat(paste0("  action: ", action))
  
  if(action=="thinning") {
    thin = FALSE
    
    enoughYearsForThinning = (args$yearsSinceThinning >= args$minThinningInterval)

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
      if(!is.na(HB)) { # HB can be missing if not enough trees > minDBH
        if(verbose) cat(paste0("  Hart-Becking: ", round(HB,1), " threshold ", round(args$thinningThreshold)))
        if((HB < args$thinningThreshold) && enoughYearsForThinning) thin = TRUE
      }
    }
    else {
      stop(paste0("Non-recognized thinning metric '", args$thinningMetric,"'.\n"))
    }
    if(thin && (sum(isTarget)>0)) {
      
      #DBH, N and N_cut of target trees
      dbhTarget = x$treeData$DBH[isTarget]
      nTarget = x$treeData$N[isTarget]
      nCutTarget = N_tree_cut[isTarget]
      
      BAtarget = sum(pi*(dbhTarget/200)^2*nTarget)
      
      # Cannot remove more than BA of target species
      BA2remove = min(BAtarget, BAtotal*(args$thinningPerc/100))
      BAremoved = 0
      args$yearsSinceThinning = 1
      
      if(verbose) cat(paste0(", type: ",args$thinning,", BA to extract: ", round(BA2remove,1), ", target cohorts: ", sum(isTarget)))
      
      if(args$thinning %in% c("below", "above")) {
        o = order(dbhTarget, decreasing = (args$thinning=="above"))
        cohort = 1
        while(BA2remove > 0) {
          r = (dbhTarget[o[cohort]]/200)
          n = nTarget[o[cohort]]
          BAcohort = pi*(r^2)*n
          if(BAcohort > BA2remove) {
            nCutTarget[o[cohort]] = BA2remove/(pi*(r^2))
            BA2remove = 0
          } else {
            nCutTarget[o[cohort]] = n
            BA2remove = BA2remove - BAcohort
          }
          cohort = cohort + 1
        }
      }
      else if (args$thinning%in% c("below-systematic", "above-systematic")) { #Cut half of target BA from below/above and the other half using trees from all diameters (keep current distribution)
        # Remove half as systematic
        sec = pi*(dbhTarget/200)^2
        propN  = nTarget/sum(nTarget)
        HalfBA2remove = BA2remove/2
        while(HalfBA2remove > 0) {
          if(sum(sec*propN)> HalfBA2remove) { # Correct to avoid cutting above existences
            propN = propN*(HalfBA2remove/sum(sec*propN))
          }
          nCutTarget = nCutTarget + propN
          HalfBA2remove = HalfBA2remove - sum(sec*propN)
        }
        BA2remove = BA2remove/2
        #Remove remaining as below/above
        o = order(dbhTarget, decreasing = (args$thinning=="above-systematic"))
        cohort = 1
        while(BA2remove > 0) {
          r = (dbhTarget[o[cohort]]/200)
          n = nTarget[o[cohort]] - nCutTarget[o[cohort]]
          BAcohort = pi*(r^2)*n
          if(BAcohort > BA2remove) {
            nCutTarget[o[cohort]] =  nCutTarget[o[cohort]] + BA2remove/(pi*(r^2))
            BA2remove = 0
          } else {
            nCutTarget[o[cohort]] = nCutTarget[o[cohort]]+n
            BA2remove = BA2remove - BAcohort
          }
          cohort = cohort + 1
        }
        #Copy cut vector of target trees to overall cut vector
        N_tree_cut[isTarget] = nCutTarget
      }
      else if (args$thinning=="systematic") { #Cut trees from all diameters (keep current distribution)
        sec = pi*(dbhTarget/200)^2
        propN  = nTarget/sum(nTarget)
        while(BA2remove > 0) {
          if(sum(sec*propN)> BA2remove) { # Correct to avoid cutting above existences
            propN = propN*(BA2remove/sum(sec*propN))
          }
          nCutTarget = nCutTarget + propN
          BA2remove = BA2remove - sum(sec*propN)
        }
        #Copy cut vector of target trees to overall cut vector
        N_tree_cut[isTarget] = nCutTarget
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
        breaks[ncl] = max(breaks[ncl], max(dbhTarget)+1) # To include largest trees
        breaks = c(0, breaks)
        # print(breaks)
        props = props/sum(props)
        cl = as.numeric(cut(dbhTarget, breaks))
        propNmat = matrix(0, sum(isTarget), ncol = ncl)
        sec = pi*(dbhTarget/200)^2
        BAcl = rep(0, ncl)
        for(i in 1:ncl) {
          # print(sum(cl==i))
          propNmat[cl==i,i] = nTarget[cl==i]/sum(nTarget[cl==i])
          BAcl[i] = sum(sec[cl==i]*nTarget[cl==i])
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
            nCutTarget = nCutTarget + propN
            BA2removecl[i] = BA2removecl[i] - sum(sec*propN)
          }
        }
      }
      
      #Copy cut vector of target trees to overall cut vector
      N_tree_cut[isTarget] = nCutTarget
    }
    else {
      action = "none"
      args$yearsSinceThinning = args$yearsSinceThinning + 1
      
    }
  }
  else if(action=="finalcut") {
    
    #DBH, N and N_cut of target trees
    dbhTarget = x$treeData$DBH[isTarget]
    nTarget = x$treeData$N[isTarget]
    nCutTarget = N_tree_cut[isTarget]
    
    BAtarget = sum(pi*(dbhTarget/200)^2*nTarget)
    
    #Proportions in different stages
    if(verbose) cat(paste0(", previous stage: ",args$finalPreviousStage))
    args$finalPreviousStage = args$finalPreviousStage + 1
    args$finalYearsToCut = args$finalYearsBetweenCuts
    props = as.numeric(strsplit(args$finalPerc,"-")[[1]])/100
    
    #BA to extract
    BA2remove = BAtarget*props[args$finalPreviousStage]
    
    if(verbose) cat(paste0(", % to extract: ", round(props[args$finalPreviousStage]*100),", BA to extract: ", round(BA2remove,1)))
    BAremoved = 0
    #Cut from above in final cuts
    o = order(dbhTarget, decreasing = TRUE)
    cohort = 1
    while((BA2remove > 0) && (cohort<=length(o))) {
      # print(BA2remove)
      d = dbhTarget[o[cohort]]
      r = (d/200)
      n = nTarget[o[cohort]]
      if((d > 12.5) || (!is.na(args$plantingSpecies))) { # Do not cut regeneration if there is no planting
        BAcohort = pi*(r^2)*n
        if(BAcohort > BA2remove) {
          nCutTarget[o[cohort]] = BA2remove/(pi*(r^2))
          BA2remove = 0
        } else {
          nCutTarget[o[cohort]] = n
          BA2remove = BA2remove - BAcohort
        }
      }
      cohort = cohort + 1
    }
    
    #Copy cut vector of target trees to overall cut vector
    N_tree_cut[isTarget] = nCutTarget
    
    if(args$finalPreviousStage==length(props)) {
      args$finalPreviousStage = 0
      if(!is.na(args$plantingSpecies)) planting = TRUE
    }
    if(verbose) cat(paste0(", final stage: ",args$finalPreviousStage))
  }
  
  # If tree felling occurred, apply understory clearing
  if(sum(N_tree_cut)>0) {
    understory_max_cover <- args$understoryMaximumCover
    if(!is.na(understory_max_cover)) {
      stand_shrub_cover = sum(x$shrubData$Cover, na.rm=TRUE)
      stand_to_rem = pmax(0,stand_shrub_cover - understory_max_cover)
      Cover_shrub_cut = x$shrubData$Cover*(stand_to_rem/stand_shrub_cover)
    }
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

#' @rdname defaultManagementFunction
defaultManagementArguments<-function(){
  return(list(
    type = "irregular",
    targetTreeSpecies = "all",
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