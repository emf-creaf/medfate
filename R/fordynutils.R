.checkForestColumns <-function(forest) {
  ntree <- nrow(forest$treeData)
  if(!("Z100" %in% names(forest$treeData))) forest$treeData$Z100 <- as.numeric(rep(NA, ntree))
  if(!("Age" %in% names(forest$treeData))) forest$treeData$Age <- as.numeric(rep(NA, ntree))
  if(!("ObsID" %in% names(forest$treeData))) forest$treeData$ObsID <- as.character(rep(NA, ntree))
  forest$treeData <- forest$treeData[,c("Species", "N", "DBH", "Height", "Z50", "Z95","Z100", "Age", "ObsID")]
  nshrub <- nrow(forest$shrubData)
  if(!("Z100" %in% names(forest$shrubData))) forest$shrubData$Z100 <- as.numeric(rep(NA, nshrub))
  if(!("Age" %in% names(forest$shrubData))) forest$shrubData$Age <- as.numeric(rep(NA, nshrub))
  if(!("ObsID" %in% names(forest$shrubData))) forest$shrubData$ObsID <- as.character(rep(NA, nshrub))
  forest$shrubData <- forest$shrubData[,c("Species", "Cover", "Height", "Z50", "Z95","Z100", "Age", "ObsID")]
  return(forest)
}
.litterUpdate <- function(initialLitterData, xo) {
  numCohorts <- nrow(xo$cohorts)
  litterLeaves <- data.frame(Species = xo$cohorts$Name,
                             Type = rep("leaves", numCohorts),
                             Necromass = xo$internalStructuralLitter$leaves)
  litterTwigs <- data.frame(Species = xo$cohorts$Name,
                            Type = rep("twigs", numCohorts),
                            Necromass = xo$internalLitter$twigs)
  litterSmallBranches <- data.frame(Species = xo$cohorts$Name,
                                    Type = rep("smallbranches", numCohorts),
                                    Necromass = xo$internalLitter$smallbranches)
  litterData <- data.frame(Species = character(0),
                           Type = character(0),
                           Necromass = numeric(0))
  for(i in 1:numCohorts) {
    r <- which((litterData$Species == litterLeaves$Species[i]) & (litterData$Type == litterLeaves$Type[i]))
    if(length(r)==1) {
      litterData$Necromass[r] <-litterData$Necromass[r] + litterLeaves$Necromass[i]
    } else {
      litterData <- rbind(litterData, litterLeaves[i, , drop = FALSE])
    }
    r <- which((litterData$Species == litterTwigs$Species[i]) & (litterData$Type == litterTwigs$Type[i]))
    if(length(r)==1) {
      litterData$Necromass[r] <-litterData$Necromass[r] + litterTwigs$Necromass[i]
    } else {
      litterData <- rbind(litterData, litterTwigs[i, , drop = FALSE])
    }
    r <- which((litterData$Species == litterSmallBranches$Species[i]) & (litterData$Type == litterSmallBranches$Type[i]))
    if(length(r)==1) {
      litterData$Necromass[r] <-litterData$Necromass[r] + litterSmallBranches$Necromass[i]
    } else {
      litterData <- rbind(litterData, litterSmallBranches[i, , drop = FALSE])
    }
  }
  # Add existing litter
  if(!is.null(initialLitterData)) {
    for(i in 1:nrow(initialLitterData)) {
      r <- which((litterData$Species == initialLitterData$Species[i]) & (litterData$Type == initialLitterData$Type[i]))
      if(length(r)==1) {
        litterData$Necromass[r] <-litterData$Necromass[r] + initialLitterData$Necromass[i]
      } else {
        litterData <- rbind(litterData, initialLitterData[i, , drop = FALSE])
      }
    }
  }
  return(litterData)
}
.nextYearForest<-function(forest, xo, SpParams, control,
                          planted_forest = emptyforest(addcolumns = c("Z100", "Age", "ObsID")),
                          recr_forest = emptyforest(addcolumns = c("Z100", "Age", "ObsID")),
                          resp_forest = emptyforest(addcolumns = c("Z100", "Age", "ObsID"))) {
  
  #Check completeness of columns in forest objects
  forest <- .checkForestColumns(forest)
  planted_forest <- .checkForestColumns(planted_forest)
  recr_forest <- .checkForestColumns(recr_forest)
  resp_forest <- .checkForestColumns(resp_forest)
  
  #Determine offsets
  treeOffset <- 0
  shrubOffset <- 0
  coh_names <- row.names(xo$cohorts)
  ctype <- .cohortType(coh_names)
  isTree <- ctype=="tree"
  isShrub <- ctype=="shrub"
  if(length(coh_names)>0) {
    coh_numbers <- sapply(strsplit(coh_names, "_"), function(x) {
      as.numeric(substr(x[[1]], 2, nchar(x[[1]])))
    })
    if(sum(isTree)>0) treeOffset <- max(coh_numbers[isTree])
    if(sum(isShrub)>0) shrubOffset <- max(coh_numbers[isShrub])
  } 
  
  # 1. Remove empty cohorts if required
  emptyTrees <- rep(FALSE, nrow(forest$treeData))
  emptyShrubs <- rep(FALSE, nrow(forest$shrubData))
  # print(forest)
  if(control$removeEmptyCohorts) {
    emptyTrees <- (forest$treeData$N < control$minimumTreeCohortDensity)
    if(control$keepCohortsWithObsID) { # Exclude from removal if cohort has non-missing ObsID
      if("ObsID" %in% names(forest$treeData)) {
        emptyTrees[!is.na(forest$treeData)] <- FALSE
      }
    }
    if(control$shrubDynamics) {
      emptyShrubs <- (forest$shrubData$Cover < control$minimumShrubCohortCover)
      if(control$keepCohortsWithObsID) { # Exclude from removal if cohort has non-missing ObsID
        if("ObsID" %in% names(forest$shrubData)) {
          emptyShrubs[!is.na(forest$shrubData)] <- FALSE
        }
      }
    }
  }
  emptyCohorts <- c(emptyTrees, emptyShrubs)
  
  # print(emptyTrees)
  # print(emptyShrubs)
  if(sum(emptyCohorts)>0) {
    forest$treeData <- forest$treeData[!emptyTrees,, drop=FALSE] 
    forest$shrubData <- forest$shrubData[!emptyShrubs,, drop=FALSE] 
    # Remove from growth input object
    xo$cohorts <- xo$cohorts[!emptyCohorts, , drop=FALSE] 
    xo$above <- xo$above[!emptyCohorts, , drop=FALSE] 
    xo$below <- xo$below[!emptyCohorts, , drop=FALSE] 
    xo$belowLayers$V <- xo$belowLayers$V[!emptyCohorts, , drop=FALSE] 
    xo$belowLayers$L <- xo$belowLayers$L[!emptyCohorts, , drop=FALSE] 
    if(control$transpirationMode!="Granier") {
      xo$belowLayers$VGrhizo_kmax <- xo$belowLayers$VGrhizo_kmax[!emptyCohorts, , drop=FALSE]
      xo$belowLayers$VCroot_kmax <- xo$belowLayers$VCroot_kmax[!emptyCohorts, , drop=FALSE]
      xo$belowLayers$RhizoPsi <- xo$belowLayers$RhizoPsi[!emptyCohorts, , drop=FALSE]
    }
    xo$belowLayers$Wpool <- xo$belowLayers$Wpool[!emptyCohorts, , drop=FALSE]
    xo$paramsPhenology <- xo$paramsPhenology[!emptyCohorts,, drop=FALSE]
    xo$paramsAnatomy <- xo$paramsAnatomy[!emptyCohorts,, drop=FALSE]
    xo$paramsInterception <- xo$paramsInterception[!emptyCohorts,, drop=FALSE]
    xo$paramsTranspiration <- xo$paramsTranspiration[!emptyCohorts,, drop=FALSE]
    xo$paramsWaterStorage <- xo$paramsWaterStorage[!emptyCohorts,, drop=FALSE]
    xo$paramsGrowth <- xo$paramsGrowth[!emptyCohorts,, drop=FALSE]
    xo$paramsAllometries <- xo$paramsAllometries[!emptyCohorts,, drop=FALSE]
    xo$internalPhenology <- xo$internalPhenology[!emptyCohorts,, drop=FALSE]
    xo$internalWater <- xo$internalWater[!emptyCohorts,, drop=FALSE]
    xo$internalCarbon <- xo$internalCarbon[!emptyCohorts,, drop=FALSE]
    xo$internalAllocation <- xo$internalAllocation[!emptyCohorts,, drop=FALSE]
    xo$internalMortality <- xo$internalMortality[!emptyCohorts,, drop=FALSE]
  }
  
  
  # 2.1 Generate above-ground data
  #planted
  planted_above <- forest2aboveground(planted_forest, SpParams, NA, TRUE)
  row.names(planted_above) <- plant_ID(planted_forest, SpParams, treeOffset, shrubOffset)
  treeOffset <- treeOffset + nrow(planted_forest$treeData)
  shrubOffset <- shrubOffset + nrow(planted_forest$shrubData)
  #recruitment
  recr_above <- forest2aboveground(recr_forest, SpParams, NA, TRUE)
  row.names(recr_above) <- plant_ID(recr_forest, SpParams, treeOffset, shrubOffset)
  treeOffset <- treeOffset + nrow(recr_forest$treeData)
  shrubOffset <- shrubOffset + nrow(recr_forest$shrubData)
  #resprouting
  resp_above <- forest2aboveground(resp_forest, SpParams, NA, TRUE)
  row.names(resp_above) <- plant_ID(resp_forest, SpParams, treeOffset, shrubOffset)
  treeOffset <- treeOffset + nrow(resp_forest$treeData)
  shrubOffset <- shrubOffset + nrow(resp_forest$shrubData)
  #surviving
  forest_above <- forest2aboveground(forest, SpParams, NA, TRUE)
  row.names(forest_above) <- row.names(xo$cohorts)
  forest_above$LAI_live[!is.na(forest_above$DBH)] <- xo$above$LAI_live[!is.na(forest_above$DBH)]
  forest_above$LAI_expanded[!is.na(forest_above$DBH)] <- xo$above$LAI_expanded[!is.na(forest_above$DBH)]
  forest_above$LAI_nocomp[!is.na(forest_above$DBH)] <- xo$above$LAI_nocomp[!is.na(forest_above$DBH)]
  if(control$shrubDynamics) {
    forest_above$LAI_live[is.na(forest_above$DBH)] <- xo$above$LAI_live[is.na(forest_above$DBH)]
    forest_above$LAI_expanded[is.na(forest_above$DBH)] <- xo$above$LAI_expanded[is.na(forest_above$DBH)]
    forest_above$LAI_nocomp[is.na(forest_above$DBH)] <- xo$above$LAI_nocomp[is.na(forest_above$DBH)]
  }

  # 2.2 Merge aboveground data (first trees)
  above_all <- rbind(forest_above[!is.na(forest_above$DBH),, drop = FALSE], 
                     planted_above[!is.na(planted_above$DBH),, drop = FALSE],
                     recr_above[!is.na(recr_above$DBH),, drop = FALSE],
                     resp_above[!is.na(resp_above$DBH),, drop = FALSE],
                     forest_above[is.na(forest_above$DBH),, drop = FALSE],
                     planted_above[is.na(planted_above$DBH),, drop = FALSE],
                     recr_above[is.na(recr_above$DBH),, drop = FALSE],
                     resp_above[is.na(resp_above$DBH),, drop = FALSE])
  
  # 3. Logical vector for replacement
  repl_vec <- c(rep(TRUE, nrow(forest$treeData)),
                rep(FALSE, nrow(planted_forest$treeData)),
                rep(FALSE, nrow(recr_forest$treeData)),
                rep(FALSE, nrow(resp_forest$treeData)),
                rep(control$shrubDynamics, nrow(forest$shrubData)),
                rep(FALSE, nrow(planted_forest$shrubData)),
                rep(FALSE, nrow(recr_forest$shrubData)),
                rep(FALSE, nrow(resp_forest$shrubData)))
  sel_vec <- c(rep(TRUE, nrow(forest$treeData)),
               rep(control$shrubDynamics, nrow(forest$shrubData)))
  
  
  # 4.1 Increase forest age for existing cohorts
  forest$treeData$Age <- forest$treeData$Age + 1
  forest$shrubData$Age <- forest$shrubData$Age + 1
  
  # 4.2 Merge cohorts in forest object
  forest$treeData <- rbind(forest$treeData, planted_forest$treeData, recr_forest$treeData, resp_forest$treeData)
  row.names(forest$treeData) <- NULL
  forest$shrubData <- rbind(forest$shrubData, planted_forest$shrubData, recr_forest$shrubData, resp_forest$shrubData)
  row.names(forest$shrubData) <- NULL

  # 4.3 Copy seedling bank from recr_forest
  if(!is.null(recr_forest$seedlingBank)) {
    forest$seedlingBank <- recr_forest$seedlingBank
  }
  
  # 4.4 Litter dynamics
  if(control$allowLitterDynamics) {
    forest$litterData <- .litterUpdate(forest$litterData, xo)
  }

  # 5.1 Prepare growth input for next year
  FCCSprops = fuel_FCCS(forest, SpParams);
  treeZ50 <- forest$treeData$Z50
  treeZ95 <- forest$treeData$Z95
  treeZ100 <- rep(NA, nrow(forest$treeData))
  if("Z100" %in% names(forest$treeData)) treeZ100 <- forest$treeData$Z100
  shrubZ50 <- forest$shrubData$Z50
  shrubZ95 <- forest$shrubData$Z95
  shrubZ100 <- rep(NA, nrow(forest$shrubData))
  if("Z100" %in% names(forest$shrubData)) shrubZ100 <- forest$shrubData$Z100
  herbZ50 <- numeric(0)
  herbZ95 <- numeric(0)
  herbZ100 <- numeric(0)
  if("herbData" %in% names(forest)) {
    herbZ50 <- forest$herbData$Z50
    herbZ95 <- forest$herbData$Z95
    herbZ100 <- rep(NA, nrow(forest$herbData))
    if("Z100" %in% names(forest$herbData)) herbZ100 <- forest$herbData$Z100
  }
  xi <- .growthInput(above = above_all, 
                     Z50 = c(treeZ50, shrubZ50, herbZ50),
                     Z95 = c(treeZ95, shrubZ95, herbZ95),
                     Z100 = c(treeZ100, shrubZ100, herbZ100),
                     xo$soil, FCCSprops, SpParams, control)
  if("herbLAI" %in% names(xo)) xi$herbLAI <- xo$herbLAI
  if("herbLAImax" %in% names(xo)) xi$herbLAImax <- xo$herbLAImax
  
  # 5.2 Replace previous state for surviving cohorts (except age)
  xi$cohorts[repl_vec,] <- xo$cohorts[sel_vec,, drop=FALSE]
  xi$above[repl_vec, (names(xi$above)!="Age")] <- xo$above[sel_vec,(names(xo$above)!="Age"), drop=FALSE]
  xi$below[repl_vec,] <- xo$below[sel_vec,, drop=FALSE]
  xi$belowLayers$V[repl_vec,] <- xo$belowLayers$V[sel_vec,, drop=FALSE]
  xi$belowLayers$L[repl_vec,] <- xo$belowLayers$L[sel_vec,, drop=FALSE]
  if(control$transpirationMode!="Granier") {
    xi$belowLayers$VGrhizo_kmax[repl_vec,] <- xo$belowLayers$VGrhizo_kmax[sel_vec,, drop=FALSE]
    xi$belowLayers$VCroot_kmax[repl_vec,] <- xo$belowLayers$VCroot_kmax[sel_vec,, drop=FALSE]
    xi$belowLayers$RhizoPsi[repl_vec,] <- xo$belowLayers$RhizoPsi[sel_vec,, drop=FALSE]
  }
  xi$belowLayers$Wpool[repl_vec,] <- xo$belowLayers$Wpool[sel_vec,, drop=FALSE]
  xi$paramsPhenology[repl_vec,] <- xo$paramsPhenology[sel_vec,, drop=FALSE]
  xi$paramsAnatomy[repl_vec,] <- xo$paramsAnatomy[sel_vec,, drop=FALSE]
  xi$paramsInterception[repl_vec,] <- xo$paramsInterception[sel_vec,, drop=FALSE]
  xi$paramsTranspiration[repl_vec,] <- xo$paramsTranspiration[sel_vec,, drop=FALSE]
  xi$paramsWaterStorage[repl_vec,] <- xo$paramsWaterStorage[sel_vec,, drop=FALSE]
  xi$paramsGrowth[repl_vec,] <- xo$paramsGrowth[sel_vec,, drop=FALSE]
  xi$paramsAllometries[repl_vec,] <- xo$paramsAllometries[sel_vec,, drop=FALSE]
  xi$internalPhenology[repl_vec,] <- xo$internalPhenology[sel_vec,, drop=FALSE]
  xi$internalWater[repl_vec,] <- xo$internalWater[sel_vec,, drop=FALSE]
  xi$internalCarbon[repl_vec,] <- xo$internalCarbon[sel_vec,, drop=FALSE]
  xi$internalAllocation[repl_vec,] <- xo$internalAllocation[sel_vec,, drop=FALSE]
  
  # This causes loss of cohort identity and reinitizalization of vegetation state variables
  # except for cohorts with ObsID when control variable 'keepCohortsWithObsID = TRUE'
  if(control$dynamicallyMergeCohorts) {
    merged_forest <- forest_mergeTrees(forest, byDBHclass = TRUE, keepCohortsWithObsID = control$keepCohortsWithObsID)
    merged_forest <- forest_mergeShrubs(merged_forest, byHeightclass = TRUE, keepCohortsWithObsID = control$keepCohortsWithObsID)
    # Only replace growth input if merging caused a reduction in woody cohorts
    if((nrow(merged_forest$treeData) < nrow(forest$treeData)) || (nrow(merged_forest$shrubData) < nrow(forest$shrubData))) {
      xi_merged <- growthInput(merged_forest, xo$soil, SpParams, control)
      if(control$keepCohortsWithObsID && ("ObsID" %in% names(xi$above))) { # Try to copy values from spared cohorts into new growth input
        # Determine sel_vec for spared cohorts
        sel_vec <- rep(FALSE, nrow(xi$above))
        sel_vec[!is.na(xi$above$ObsID)] <- TRUE
        repl_vec <- rep(FALSE, nrow(xi_merged$above))
        repl_vec[!is.na(xi_merged$above$ObsID)] <- TRUE
        
        # Replace previous state for spared cohorts
        xi_merged$cohorts[repl_vec,] <- xi$cohorts[sel_vec,, drop=FALSE]
        xi_merged$above[repl_vec,] <- xi$above[sel_vec,, drop=FALSE]
        xi_merged$below[repl_vec,] <- xi$below[sel_vec,, drop=FALSE]
        xi_merged$belowLayers$V[repl_vec,] <- xi$belowLayers$V[sel_vec,, drop=FALSE]
        xi_merged$belowLayers$L[repl_vec,] <- xi$belowLayers$L[sel_vec,, drop=FALSE]
        if(control$transpirationMode!="Granier") {
          xi_merged$belowLayers$VGrhizo_kmax[repl_vec,] <- xi$belowLayers$VGrhizo_kmax[sel_vec,, drop=FALSE]
          xi_merged$belowLayers$VCroot_kmax[repl_vec,] <- xi$belowLayers$VCroot_kmax[sel_vec,, drop=FALSE]
          xi_merged$belowLayers$RhizoPsi[repl_vec,] <- xi$belowLayers$RhizoPsi[sel_vec,, drop=FALSE]
        }
        xi_merged$belowLayers$Wpool[repl_vec,] <- xi$belowLayers$Wpool[sel_vec,, drop=FALSE]
        xi_merged$paramsPhenology[repl_vec,] <- xi$paramsPhenology[sel_vec,, drop=FALSE]
        xi_merged$paramsAnatomy[repl_vec,] <- xi$paramsAnatomy[sel_vec,, drop=FALSE]
        xi_merged$paramsInterception[repl_vec,] <- xi$paramsInterception[sel_vec,, drop=FALSE]
        xi_merged$paramsTranspiration[repl_vec,] <- xi$paramsTranspiration[sel_vec,, drop=FALSE]
        xi_merged$paramsWaterStorage[repl_vec,] <- xi$paramsWaterStorage[sel_vec,, drop=FALSE]
        xi_merged$paramsGrowth[repl_vec,] <- xi$paramsGrowth[sel_vec,, drop=FALSE]
        xi_merged$paramsAllometries[repl_vec,] <- xi$paramsAllometries[sel_vec,, drop=FALSE]
        xi_merged$internalPhenology[repl_vec,] <- xi$internalPhenology[sel_vec,, drop=FALSE]
        xi_merged$internalWater[repl_vec,] <- xi$internalWater[sel_vec,, drop=FALSE]
        xi_merged$internalCarbon[repl_vec,] <- xi$internalCarbon[sel_vec,, drop=FALSE]
        xi_merged$internalAllocation[repl_vec,] <- xi$internalAllocation[sel_vec,, drop=FALSE]
      }
      # Replace growthInput and forest by the one after merging
      forest <- merged_forest
      xi <- xi_merged
    }
  }
  return(list(forest = forest, xi = xi))
}

.summarizeCohorts<-function(step, 
                            treeTableYear, shrubTableYear,
                            deadTreeTableYear, deadShrubTableYear,
                            cutTreeTableYear, cutShrubTableYear) {
  treecohnames = unique(c(treeTableYear$Cohort,deadTreeTableYear$Cohort))
  if(!is.null(cutTreeTableYear)) treecohnames = unique(c(treecohnames, cutTreeTableYear$Cohort))
  shrubcohnames = unique(c(shrubTableYear$Cohort,deadShrubTableYear$Cohort))
  if(!is.null(cutShrubTableYear)) shrubcohnames = unique(c(shrubcohnames, cutShrubTableYear$Cohort))
  cohnames = c(treecohnames, shrubcohnames)
  isTree = c(rep(TRUE, length(treecohnames)), rep(FALSE, length(shrubcohnames)))
  cohortSummary = data.frame("Step" = rep(step, length(cohnames)),
                             "Species" = rep(NA, length(cohnames)),
                             "Cohort" = cohnames,
                             "TreeDensityLive" = rep(0, length(cohnames)),
                             "TreeBasalAreaLive"= rep(0, length(cohnames)),
                             "ShrubCoverLive"= rep(0, length(cohnames)),
                             "BasalAreaDead" = rep(0, length(cohnames)),
                             "ShrubCoverDead" = rep(0, length(cohnames)),
                             "BasalAreaCut" = rep(0, length(cohnames)),
                             "ShrubCoverCut" = rep(0, length(cohnames)))
  cohortSummary$ShrubCoverLive[isTree] = NA
  cohortSummary$ShrubCoverDead[isTree] = NA
  cohortSummary$ShrubCoverCut[isTree] = NA
  cohortSummary$TreeDensityLive[!isTree] = NA
  cohortSummary$TreeBasalAreaLive[!isTree] = NA
  cohortSummary$BasalAreaDead[!isTree] = NA
  cohortSummary$BasalAreaCut[!isTree] = NA
  
  ba_live = .treeBasalArea(treeTableYear$N, treeTableYear$DBH)
  for(i in 1:nrow(treeTableYear)) {
    icoh = which(cohortSummary$Cohort==treeTableYear$Cohort[i])
    cohortSummary$Cohort[icoh] = treeTableYear$Cohort[i]
    cohortSummary$Species[icoh] = treeTableYear$Species[i]
    cohortSummary$TreeDensityLive[icoh] = treeTableYear$N[i]
    cohortSummary$TreeBasalAreaLive[icoh] = ba_live[i]
  }
  for(i in 1:nrow(shrubTableYear)) {
    icoh = which(cohortSummary$Cohort==shrubTableYear$Cohort[i])
    cohortSummary$Cohort[icoh] = shrubTableYear$Cohort[i]
    cohortSummary$Species[icoh] = shrubTableYear$Species[i]
    cohortSummary$ShrubCoverLive[icoh] = shrubTableYear$Cover[i]
  }
  ba_dead = .treeBasalArea(deadTreeTableYear$N, deadTreeTableYear$DBH)
  for(i in 1:nrow(deadTreeTableYear)) {
    icoh = which(cohortSummary$Cohort==deadTreeTableYear$Cohort[i])
    cohortSummary$Cohort[icoh] = deadTreeTableYear$Cohort[i]
    cohortSummary$Species[icoh] = deadTreeTableYear$Species[i]
    cohortSummary$BasalAreaDead[icoh] = ba_dead[i]
  }
  for(i in 1:nrow(deadShrubTableYear)) {
    icoh = which(cohortSummary$Cohort==deadShrubTableYear$Cohort[i])
    cohortSummary$Cohort[icoh] = deadShrubTableYear$Cohort[i]
    cohortSummary$Species[icoh] = deadShrubTableYear$Species[i]
    cohortSummary$ShrubCoverDead[icoh] = deadShrubTableYear$Cover[i]
  }
  if(!is.null(cutTreeTableYear)) {
    ba_cut = .treeBasalArea(cutTreeTableYear$N, cutTreeTableYear$DBH)
    for(i in 1:nrow(cutTreeTableYear)) {
      icoh = which(cohortSummary$Cohort==cutTreeTableYear$Cohort[i])
      cohortSummary$Cohort[icoh] = cutTreeTableYear$Cohort[i]
      cohortSummary$Species[icoh] = cutTreeTableYear$Species[i]
      cohortSummary$BasalAreaCut[icoh] = ba_cut[i]
    }
  }
  if(!is.null(cutShrubTableYear)) {
    for(i in 1:nrow(cutShrubTableYear)) {
      icoh = which(cohortSummary$Cohort==cutShrubTableYear$Cohort[i])
      cohortSummary$Cohort[icoh] = cutShrubTableYear$Cohort[i]
      cohortSummary$Species[icoh] = cutShrubTableYear$Species[i]
      cohortSummary$ShrubCoverCut[icoh] = cutShrubTableYear$Cover[i]
    }
  }
  return(cohortSummary)
}
.summarizeSpecies<-function(step, cohSum, x, SpParams) {
  nl_sp = tapply(cohSum$TreeDensityLive, cohSum$Species, sum, na.rm=FALSE)
  bal_sp = tapply(cohSum$TreeBasalAreaLive, cohSum$Species, sum, na.rm=FALSE)
  bad_sp = tapply(cohSum$BasalAreaDead, cohSum$Species, sum, na.rm=FALSE)
  bac_sp = tapply(cohSum$BasalAreaCut, cohSum$Species, sum, na.rm=FALSE)
  shl_sp = tapply(cohSum$ShrubCoverLive, cohSum$Species, sum, na.rm=FALSE)
  shd_sp = tapply(cohSum$ShrubCoverDead, cohSum$Species, sum, na.rm=FALSE)
  shc_sp = tapply(cohSum$ShrubCoverCut, cohSum$Species, sum, na.rm=FALSE)
  spSumYear <-data.frame("Step" = rep(step, length(bal_sp)),
                         "Species" = names(bal_sp),
                         "NumCohorts" = as.numeric(table(cohSum$Species)),
                         "TreeDensityLive"= as.numeric(nl_sp),
                         "TreeBasalAreaLive"= as.numeric(bal_sp),
                         "ShrubCoverLive"= as.numeric(shl_sp),
                         "BasalAreaDead" = as.numeric(bad_sp),
                         "ShrubCoverDead" = as.numeric(shd_sp),
                         "BasalAreaCut" = as.numeric(bac_sp),
                         "ShrubCoverCut" = as.numeric(shc_sp))
}
.summarizeStand<-function(step, cohSum, x) {
  isTree = !is.na(x$above$DBH)
  HB = .hartBeckingIndex(x$above$N[isTree], x$above$H[isTree], x$above$DBH[isTree])
  domH = .dominantTreeHeight(x$above$N[isTree], x$above$H[isTree], x$above$DBH[isTree])
  domDBH = .dominantTreeDiameter(x$above$N[isTree], x$above$DBH[isTree])
  qmDBH = .quadraticMeanTreeDiameter(x$above$N[isTree], x$above$DBH[isTree])
  standSumYear = data.frame("Step" = step,
                            "NumTreeSpecies" = length(unique(x$above$SP[isTree])),
                            "NumTreeCohorts" = sum(isTree),
                            "NumShrubSpecies" = length(unique(x$above$SP[!isTree])),
                            "NumShrubCohorts" = sum(!isTree),
                            "TreeDensityLive"= sum(cohSum$TreeDensityLive, na.rm = TRUE),
                            "TreeBasalAreaLive"= sum(cohSum$TreeBasalAreaLive, na.rm = TRUE),
                            "DominantTreeHeight"= domH,
                            "DominantTreeDiameter"= domDBH,
                            "QuadraticMeanTreeDiameter"= qmDBH,
                            "HartBeckingIndex"= HB,
                            "ShrubCoverLive"= sum(cohSum$ShrubCoverLive, na.rm=TRUE),
                            "BasalAreaDead" = sum(cohSum$BasalAreaDead, na.rm=TRUE),
                            "ShrubCoverDead" = sum(cohSum$ShrubCoverDead, na.rm=TRUE),
                            "BasalAreaCut" = sum(cohSum$BasalAreaCut, na.rm=TRUE),
                            "ShrubCoverCut" = sum(cohSum$ShrubCoverCut, na.rm=TRUE))
}

.createTreeTable<-function(step, year, x) {
  range <- numeric(0)
  ctype  <- .cohortType(row.names(x$cohorts))
  isTree <- ctype=="tree"
  if(sum(isTree)>0) range <- 1:sum(isTree)

  tt<-data.frame(Step = rep(step, length(range)), 
                 Year = rep(year, length(range)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$cohorts$Name[range],
                 DBH = x$above$DBH[range],
                 Height = x$above$H[range],
                 N = x$above$N[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range],
                 Z100 = x$below$Z100[range])
  if("Age" %in% names(x$above)) {
    tt$Age <- x$above$Age[range]
  } else {
    tt$Age <- as.numeric(rep(NA, sum(isTree)))
  }
  if("ObsID" %in% names(x$above)) {
    tt$ObsID <- x$above$ObsID[range]
  } else {
    tt$ObsID <- as.character(rep(NA, sum(isTree)))
  }
  tt = tt[tt$N>0,, drop=FALSE]
  return(tt)
}
.createDeadTreeTable<-function(step, year, x) {
  range <- numeric(0)
  ctype  <- .cohortType(row.names(x$cohorts))
  isTree <- ctype=="tree"
  if(sum(isTree)>0) range <- 1:sum(isTree)
  
  dtt<-data.frame(Step = rep(step, length(range)), 
                  Year = rep(year, length(range)),
                  Cohort = row.names(x$cohorts)[range],
                  Species = x$cohorts$Name[range],
                  DBH = x$above$DBH[range],
                  Height = x$above$H[range],
                  N = x$internalMortality$N_dead[range],
                  N_starvation = x$internalMortality$N_starvation[range],
                  N_dessication = x$internalMortality$N_dessication[range],
                  N_burnt = x$internalMortality$N_burnt[range],
                  Z50 = x$below$Z50[range],
                  Z95 = x$below$Z95[range],
                  Z100 = x$below$Z100[range])
  if("Age" %in% names(x$above)) {
    dtt$Age <- x$above$Age[range]
  } else {
    dtt$Age <- as.numeric(rep(NA, sum(isTree)))
  }
  if("ObsID" %in% names(x$above)) {
    dtt$ObsID <- x$above$ObsID[range]
  } else {
    dtt$ObsID <- as.character(rep(NA, sum(isTree)))
  }
  dtt = dtt[dtt$N>0,, drop = FALSE]
  return(dtt)
}
.createCutTreeTable<-function(step, year, x, N_cut) {
  ctype  <- .cohortType(row.names(x$cohorts))
  isTree <- ctype=="tree"
  range <- numeric(0)
  if(length(N_cut)>0) range = 1:length(N_cut)
  ctt<-data.frame(Step = rep(step, length(N_cut)), 
                 Year = rep(year, length(N_cut)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$cohorts$Name[range],
                 DBH = x$above$DBH[range],
                 Height = x$above$H[range],
                 N = N_cut,
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range],
                 Z100 = x$below$Z100[range])
  if("Age" %in% names(x$above)) {
    ctt$Age <- x$above$Age[range]
  } else {
    ctt$Age <- as.numeric(rep(NA, sum(isTree)))
  }
  if("ObsID" %in% names(x$above)) {
    ctt$ObsID <- x$above$ObsID[range]
  } else {
    ctt$ObsID <- as.character(rep(NA, sum(isTree)))
  }
  ctt = ctt[ctt$N>0,, drop=FALSE]
  return(ctt)
}
.createShrubTable<-function(step, year, x) {
  ctype  <- .cohortType(row.names(x$cohorts))
  nt <- sum(ctype=="tree")
  ns <- sum(ctype=="shrub")
  range <- numeric(0)
  if(ns>0) range = (nt+1):(nt+ns)
  st<-data.frame(Step = rep(step, ns), 
                 Year = rep(year, ns),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$cohorts$Name[range],
                 Cover = x$above$Cover[range],
                 Height = x$above$H[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range],
                 Z100 = x$below$Z100[range])
  if("Age" %in% names(x$above)) {
    st$Age <- x$above$Age[range]
  } else {
    st$Age <- as.numeric(rep(NA, ns))
  }
  if("ObsID" %in% names(x$above)) {
    st$ObsID <- x$above$ObsID[range]
  } else {
    st$ObsID <- as.character(rep(NA, ns))
  }
  st = st[st$Cover>0,, drop = FALSE]
  return(st)
}
.createDeadShrubTable<-function(step, year, x) {
  ctype  <- .cohortType(row.names(x$cohorts))
  nt = sum(ctype=="tree")
  ns <- sum( ctype=="shrub")
  range <- numeric(0)
  if(ns>0) range = (nt+1):(nt+ns)
  dst<-data.frame(Step = rep(step, ns), 
                  Year = rep(year, ns),
                  Cohort = row.names(x$cohorts)[range],
                  Species = x$cohorts$Name[range],
                  Cover = x$internalMortality$Cover_dead[range],
                  Cover_starvation = x$internalMortality$Cover_starvation[range],
                  Cover_dessication = x$internalMortality$Cover_dessication[range],
                  Cover_burnt = x$internalMortality$Cover_burnt[range],
                  Height = x$above$H[range],
                  Z50 = x$below$Z50[range],
                  Z95 = x$below$Z95[range],
                  Z100 = x$below$Z100[range])
  if("Age" %in% names(x$above)) {
    dst$Age <- x$above$Age[range]
  } else {
    dst$Age <- as.numeric(rep(NA, ns))
  }
  if("ObsID" %in% names(x$above)) {
    dst$ObsID <- x$above$ObsID[range]
  } else {
    dst$ObsID <- as.character(rep(NA, ns))
  }
  dst = dst[dst$Cover>0,,drop=FALSE]
  return(dst)
}
.createCutShrubTable<-function(step, year, x, Cover_cut) {
  ctype  <- .cohortType(row.names(x$cohorts))
  nt = sum(ctype=="tree")
  ns <- sum( ctype=="shrub")
  range <- numeric(0)
  if(ns>0) range = (nt+1):(nt+ns)
  cst<-data.frame(Step = rep(step, length(range)), 
                 Year = rep(year, length(range)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$cohorts$Name[range],
                 Cover = Cover_cut,
                 Height = x$above$H[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range],
                 Z100 = x$below$Z100[range])
  if("Age" %in% names(x$above)) {
    cst$Age <- x$above$Age[range]
  } else {
    cst$Age <- as.numeric(rep(NA, ns))
  }
  if("ObsID" %in% names(x$above)) {
    cst$ObsID <- x$above$ObsID[range]
  } else {
    cst$ObsID <- as.character(rep(NA, ns))
  }
  cst = cst[cst$Cover>0,, drop = FALSE]
  return(cst)
}