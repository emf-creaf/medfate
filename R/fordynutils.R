.nextYearForest<-function(forest, xo, SpParams, control,
                          planted_forest = emptyforest(),
                          recr_forest = emptyforest(),
                          resp_forest = emptyforest()) {
  
  #Determine offsets
  treeOffset <- 0
  shrubOffset <- 0
  coh_names <- row.names(xo$cohorts)
  isTree <- substr(coh_names,1,1)=="T"
  if(length(coh_names)>0) {
    coh_numbers <- sapply(strsplit(coh_names, "_"), function(x) {
      as.numeric(substr(x[[1]], 2, nchar(x[[1]])))
    })
    if(sum(isTree)>0) treeOffset <- max(coh_numbers[isTree])
    if(sum(!isTree)>0) shrubOffset <- max(coh_numbers[!isTree])
  } 
  
  # 1. Remove empty cohorts if required
  emptyTrees <- rep(FALSE, nrow(forest$treeData))
  emptyShrubs <- rep(FALSE, nrow(forest$shrubData))
  if(control$removeEmptyCohorts) {
    emptyTrees <- (forest$treeData$N < control$minimumTreeCohortDensity)
    if(control$shrubDynamics) emptyShrubs <- (forest$shrubData$Cover < control$minimumShrubCohortCover)
  }
  emptyCohorts <- c(emptyTrees, emptyShrubs)
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
  if(control$shrubDynamics) {
    forest_above$LAI_live[is.na(forest_above$DBH)] <- xo$above$LAI_live[is.na(forest_above$DBH)]
    forest_above$LAI_expanded[is.na(forest_above$DBH)] <- xo$above$LAI_expanded[is.na(forest_above$DBH)]
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
  
  # 4. Merge cohorts in forest object
  forest$treeData <- rbind(forest$treeData, planted_forest$treeData, recr_forest$treeData, resp_forest$treeData)
  forest$shrubData <- rbind(forest$shrubData, planted_forest$shrubData, recr_forest$shrubData, resp_forest$shrubData)
  
  
  # 5.1 Prepare growth input for next year
  FCCSprops = fuel_FCCS(forest, SpParams);
  xi <- .growthInput(above = above_all, 
                     Z50 = c(forest$treeData$Z50, forest$shrubData$Z50),
                     Z95 = c(forest$treeData$Z95, forest$shrubData$Z95),
                     xo$soil, FCCSprops, SpParams, control)
  xi$herbLAI <- xo$herbLAI
  xi$herbLAImax <- xo$herbLAImax
  
  # 5.2 Replace previous state for surviving cohorts
  xi$cohorts[repl_vec,] <- xo$cohorts[sel_vec,, drop=FALSE]
  xi$above[repl_vec,] <- xo$above[sel_vec,, drop=FALSE]
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
  if(control$dynamicallyMergeCohorts) {
    merged_forest <- forest_mergeTrees(forest, byDBHclass = TRUE)
    merged_forest <- forest_mergeShrubs(merged_forest, byHeightclass = TRUE)
    # Only replace if merging caused a reduction in woody cohorts
    if((nrow(merged_forest$treeData) < nrow(forest$treeData)) || (nrow(merged_forest$shrubData) < nrow(forest$shrubData))) {
      forest <- merged_forest
      xi <- forest2growthInput(forest, xo$soil, SpParams, control)
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
  isTree <- !is.na(x$above$DBH)
  if(sum(isTree)>0) range <- 1:sum(isTree)

  tt<-data.frame(Step = rep(step, length(range)), 
                 Year = rep(year, length(range)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$cohorts$Name[range],
                 DBH = x$above$DBH[range],
                 Height = x$above$H[range],
                 N = x$above$N[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range])
  tt = tt[tt$N>0,, drop=FALSE]
  return(tt)
}
.createDeadTreeTable<-function(step, year, x) {
  range <- numeric(0)
  isTree <- !is.na(x$above$DBH)
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
                  Z95 = x$below$Z95[range])
  dtt = dtt[dtt$N>0,, drop = FALSE]
  return(dtt)
}
.createCutTreeTable<-function(step, year, x, N_cut) {
  range = numeric(0)
  if(length(N_cut)>0) range = 1:length(N_cut)
  tt<-data.frame(Step = rep(step, length(N_cut)), 
                 Year = rep(year, length(N_cut)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$cohorts$Name[range],
                 DBH = x$above$DBH[range],
                 Height = x$above$H[range],
                 N = N_cut,
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range])
  tt = tt[tt$N>0,, drop=FALSE]
  return(tt)
}
.createShrubTable<-function(step, year, x) {
  isShrub = !is.na(x$above$Cover)
  nt = sum(!isShrub)
  numCohorts = length(isShrub)
  range = numeric(0)
  if(numCohorts>nt) range = (nt+1):numCohorts
  st<-data.frame(Step = rep(step, length(range)), 
                 Year = rep(year, length(range)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$cohorts$Name[range],
                 Cover = x$above$Cover[range],
                 Height = x$above$H[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range])
  st = st[st$Cover>0,, drop = FALSE]
  return(st)
}
.createDeadShrubTable<-function(step, year, x) {
  isShrub = !is.na(x$above$Cover)
  nt = sum(!isShrub)
  numCohorts = length(isShrub)
  range = numeric(0)
  if(numCohorts>nt) range = (nt+1):numCohorts
  dst<-data.frame(Step = rep(step, length(range)), 
                  Year = rep(year, length(range)),
                  Cohort = row.names(x$cohorts)[range],
                  Species = x$cohorts$Name[range],
                  Cover = x$internalMortality$Cover_dead[range],
                  Cover_starvation = x$internalMortality$Cover_starvation[range],
                  Cover_dessication = x$internalMortality$Cover_dessication[range],
                  Cover_burnt = x$internalMortality$Cover_burnt[range],
                  Height = x$above$H[range],
                  Z50 = x$below$Z50[range],
                  Z95 = x$below$Z95[range])
  dst = dst[dst$Cover>0,,drop=FALSE]
  return(dst)
}
.createCutShrubTable<-function(step, year, x, Cover_cut) {
  isShrub = !is.na(x$above$Cover)
  nt = sum(!isShrub)
  numCohorts = length(isShrub)
  range = numeric(0)
  if(numCohorts>nt) range = (nt+1):numCohorts
  st<-data.frame(Step = rep(step, length(range)), 
                 Year = rep(year, length(range)),
                 Cohort = row.names(x$cohorts)[range],
                 Species = x$cohorts$Name[range],
                 Cover = Cover_cut,
                 Height = x$above$H[range],
                 Z50 = x$below$Z50[range],
                 Z95 = x$below$Z95[range])
  st = st[st$Cover>0,, drop = FALSE]
  return(st)
}