dynamics<-function(forest, soil, SpParams,
                   meteo, control,
                   latitude , elevation = NA, slope = NA, aspect = NA) {
  
  control$modifyInput = FALSE
  control$verbose = FALSE
  
  dates = as.Date(row.names(meteo))
  years = as.numeric(format(dates, "%Y"))
  yearsUnique = unique(years)
  nYears = length(yearsUnique)
  growthResults = vector("list", nYears)
  names(growthResults) = paste0("Step_", 1:nYears)
  forestStructures = vector("list", nYears+1)
  names(forestStructures) = c("Initial", paste0("Step_", 1:nYears))
  forestStructures[[1]] = forest


  summarizeCohorts<-function(step, x) {
    cohortSummary = data.frame("Step" = rep(step, nrow(x$cohorts)),
                               "Species" = x$cohorts$SP,
                               "Name" = x$cohorts$Name,
                               "Cohort" = row.names(x$cohorts),
                               "LeafAreaIndex" = x$above$LAI_expanded,
                               "DensityLive" = x$above$N,
                               "TreeBasalAreaLive"= x$above$N*pi*(x$above$DBH/200)^2,
                               "ShrubCoverLive"= x$above$Cover)
    return(cohortSummary)
  }
  summarizeSpecies<-function(step, cohSum, x) {
    lai_sp = tapply(cohSum$LeafAreaIndex, cohSum$Species, sum, na.rm=TRUE)
    nl_sp = tapply(cohSum$DensityLive, cohSum$Species, sum, na.rm=TRUE)
    bal_sp = tapply(cohSum$TreeBasalAreaLive, cohSum$Species, sum, na.rm=TRUE)
    shl_sp = tapply(cohSum$ShrubCoverLive, cohSum$Species, sum, na.rm=TRUE)
    mh_sp = tapply(x$above$H, x$above$SP, max, na.rm=TRUE)
    spSumYear <-data.frame("Step" = rep(step, length(lai_sp)),
                           "Species" = as.numeric(names(lai_sp)),
                           "LeafAreaIndex" = as.numeric(lai_sp),
                           "DensityLive"= as.numeric(nl_sp),
                           "TreeBasalAreaLive"= as.numeric(bal_sp),
                           "ShrubCoverLive"= as.numeric(shl_sp),
                           "MaxHeight"= as.numeric(mh_sp))
  }
  summarizeStand<-function(step, cohSum, x) {
    standSumYear = data.frame("Step" = step,
                              "LeafAreaIndex" = sum(cohSum$LeafAreaIndex, na.rm = TRUE),
                              "DensityLive"= sum(cohSum$DensityLive, na.rm = TRUE),
                              "TreeBasalAreaLive"= sum(cohSum$TreeBasalAreaLive, na.rm = TRUE),
                              "ShrubCoverLive"= sum(cohSum$ShrubCoverLive, na.rm=TRUE),
                              "MaxHeight"= max(x$above$H, na.rm=TRUE))
  }
  createTreeTable<-function(step, year, x) {
    isTree = !is.na(x$above$DBH)
    nt = sum(isTree)
    tt<-data.frame(Step = step, Year = year,
                   Cohort = row.names(x$cohorts)[1:nt],
                   Species = x$above$SP[1:nt],
                   Name = x$cohorts$Name[1:nt],
                   N = x$above$N[1:nt],
                   DBH = x$above$DBH[1:nt],
                   Height = x$above$H[1:nt])
    if(control$removeDeadCohorts) tt = tt[tt$N>control$minimumDensity,]
    return(tt)
  }
  createShrubTable<-function(step, year, x) {
    isShrub = !is.na(x$above$Cover)
    nt = sum(!isShrub)
    numCohorts = length(isShrub)
    range = numeric(0)
    if(numCohorts>nt) range = (nt+1):numCohorts
    st<-data.frame(Step = step, Year = year,
                   Cohort = row.names(x$cohorts)[range],
                   Species = x$above$SP[range],
                   Name = x$cohorts$Name[range],
                   Cover = x$above$Cover[range],
                   Height = x$above$H[range])
    if(control$removeDeadCohorts) st = st[st$N>control$minimumDensity,]
    return(st)
  }
  
  #Initialization
  treeOffset = nrow(forest$treeData)
  shrubOffset = nrow(forest$shrubData)
  xi = forest2growthInput(forest, soil, SpParams, control)
  
  #initial summaries
  cohortSummary <-summarizeCohorts(0, xi)
  speciesSummary<-summarizeSpecies(0,cohortSummary, xi)
  standSummary<-summarizeStand(0,cohortSummary, xi)

  #initial tree/shrub tables
  treeTable = createTreeTable(0, NA, xi)
  shrubTable = createShrubTable(0, NA, xi)
  
  #Simulations
  for(iYear in 1:nYears) {
    year = yearsUnique[iYear]
    cat(paste0("Simulating forest dynamics for year ", year, " (", iYear,"/", nYears,")\n"))
    meteoYear = meteo[years==year,]
    cat(paste0("   (a) Growth/mortality\n"))
    Gi = growth(xi, meteoYear, latitude = latitude, elevation = elevation, slope = slope, aspect = aspect)
    
    # Store growth results
    growthResults[[iYear]] = Gi
    
    # Retrieve modified growth output
    xo = Gi$growthInput
    
    # Update forest structural variables
    isTree = is.na(xo$above$Cover)
    forest$treeData$N  = xo$above$N[isTree]
    forest$treeData$DBH  = xo$above$DBH[isTree]
    forest$treeData$Height  = xo$above$H[isTree]
    if(control$shrubDynamics) {
      forest$shrubData$Cover  = xo$above$Cover[!isTree]
      forest$shrubData$Height  = xo$above$H[!isTree]
    }
    
    # Remove dead cohorts if required
    deadTrees = rep(FALSE, nrow(forest$treeData))
    deadShrubs = rep(FALSE, nrow(forest$shrubData))
    if(control$removeDeadCohorts) {
      deadTrees = (forest$treeData$N < control$minimumDensity)
      if(control$shrubDynamics) deadShrubs = (forest$shrubData$Cover < control$minimumDensity)
    }
    deadCohorts = c(deadTrees, deadShrubs)
    if(sum(deadCohorts)>0) {
      cat(paste0("   (-) Removing dead cohorts: ", paste(row.names(xo$above)[deadCohorts], collapse=","),"\n"))
      forest$treeData = forest$treeData[!deadTrees,, drop=FALSE] 
      forest$shrubData = forest$shrubData[!deadShrubs,, drop=FALSE] 
      # Remove from growth input object
      xo$cohorts = xo$cohorts[!deadCohorts, , drop=FALSE] 
      xo$above = xo$above[!deadCohorts, , drop=FALSE] 
      xo$below = xo$below[!deadCohorts, , drop=FALSE] 
      xo$belowLayers$V = xo$belowLayers$V[!deadCohorts, , drop=FALSE] 
      xo$belowLayers$L = xo$belowLayers$L[!deadCohorts, , drop=FALSE] 
      if(control$transpirationMode=="Sperry") {
        xo$belowLayers$VGrhizo_kmax = xo$belowLayers$VGrhizo_kmax[!deadCohorts, , drop=FALSE]
        xo$belowLayers$VCroot_kmax = xo$belowLayers$VCroot_kmax[!deadCohorts, , drop=FALSE]
        xo$belowLayers$RhizoPsi = xo$belowLayers$RhizoPsi[!deadCohorts, , drop=FALSE]
      }
      xo$belowLayers$Wpool = xo$belowLayers$Wpool[!deadCohorts, , drop=FALSE]
      xo$paramsPhenology = xo$paramsPhenology[!deadCohorts,, drop=FALSE]
      xo$paramsAnatomy = xo$paramsAnatomy[!deadCohorts,, drop=FALSE]
      xo$paramsInterception = xo$paramsInterception[!deadCohorts,, drop=FALSE]
      xo$paramsTranspiration = xo$paramsTranspiration[!deadCohorts,, drop=FALSE]
      xo$paramsWaterStorage = xo$paramsWaterStorage[!deadCohorts,, drop=FALSE]
      xo$paramsGrowth = xo$paramsGrowth[!deadCohorts,, drop=FALSE]
      xo$paramsAllometries = xo$paramsAllometries[!deadCohorts,, drop=FALSE]
      xo$internalPhenology = xo$internalPhenology[!deadCohorts,, drop=FALSE]
      xo$internalWater = xo$internalWater[!deadCohorts,, drop=FALSE]
      xo$internalCarbon = xo$internalCarbon[!deadCohorts,, drop=FALSE]
      xo$internalAllocation = xo$internalAllocation[!deadCohorts,, drop=FALSE]
      xo$internalMortality = xo$internalMortality[!deadCohorts,, drop=FALSE]
      xo$internalRings = xo$internalRings[!deadCohorts]
    }

    
    # Simulate species recruitment
    cat(paste0("   (b) Recruitment\n"))
    treeSpp = unique(forest$treeData$Species)
    shrubSpp = unique(forest$shrubData$Species)
    recr_forest = emptyforest(ntree = length(treeSpp), nshrub=length(shrubSpp))
    if(length(treeSpp)>0) {
      recr_forest$treeData$Species = treeSpp
      recr_forest$treeData$N = 100
      recr_forest$treeData$DBH = 1
      recr_forest$treeData$Height = 100
      for(i in 1:length(treeSpp)) {
        j = which(forest$treeData$Species==treeSpp[i])[1]
        recr_forest$treeData$Z50[i] = forest$treeData$Z50[j]
        recr_forest$treeData$Z95[i] = forest$treeData$Z95[j]
      }
    }
    if((length(shrubSpp)>0)) {
      recr_forest$shrubData$Species = shrubSpp
      recr_forest$shrubData$Cover = 1
      recr_forest$shrubData$Height = 10
      for(i in 1:length(shrubSpp)) {
        j = which(forest$shrubData$Species==shrubSpp[i])[1]
        recr_forest$shrubData$Z50[i] = forest$shrubData$Z50[j]
        recr_forest$shrubData$Z95[i] = forest$shrubData$Z95[j]
      }
    }
    if(!control$shrubDynamics) recr_forest$shrubData = recr_forest$shrubData[numeric(0),, drop=FALSE]
    
    # Generate above-ground data
    recr_above = forest2aboveground(recr_forest, SpParams, NA, "MED")
    row.names(recr_above) = plant_ID(recr_forest, treeOffset, shrubOffset)
    treeOffset = treeOffset + nrow(recr_forest$treeData)
    shrubOffset = shrubOffset + nrow(recr_forest$shrubData)
    forest_above = forest2aboveground(forest, SpParams, NA, "MED")
    row.names(forest_above) = row.names(xo$cohorts)
    forest_above$LAI_live = xo$above$LAI_live
    forest_above$LAI_expanded = xo$above$LAI_expanded

    # Merge above-ground data (first trees)
    above_all = rbind(forest_above[!is.na(forest_above$DBH),, drop = FALSE], 
                      recr_above[!is.na(recr_above$DBH),, drop = FALSE],
                      forest_above[is.na(forest_above$DBH),, drop = FALSE], 
                      recr_above[is.na(recr_above$DBH),, drop = FALSE])
    
    # Logical vector for replacement
    repl_vec <- c(rep(TRUE, nrow(forest$treeData)),
                  rep(FALSE, nrow(recr_forest$treeData)),
                  rep(TRUE, nrow(forest$shrubData)),
                  rep(FALSE, nrow(recr_forest$shrubData)))
    
    # Merge cohorts in forest object
    forest$treeData = rbind(forest$treeData, recr_forest$treeData)
    forest$shrubData = rbind(forest$shrubData, recr_forest$shrubData)
    
    
    # Prepare growth input for next year
    xi = growthInput(above = above_all,
                     Z50 = c(forest$treeData$Z50, forest$shrubData$Z50),
                     Z95 = c(forest$treeData$Z95, forest$shrubData$Z95),
                     xo$soil, SpParams, control)
    
    # Replace previous state for surviving cohorts
    xi$cohorts[repl_vec,] <- xo$cohorts
    xi$above[repl_vec,] <- xo$above
    xi$below[repl_vec,] <- xo$below
    xi$belowLayers$V[repl_vec,] <- xo$belowLayers$V
    xi$belowLayers$L[repl_vec,] <- xo$belowLayers$L
    if(control$transpirationMode=="Sperry") {
      xi$belowLayers$VGrhizo_kmax[repl_vec,] <- xo$belowLayers$VGrhizo_kmax
      xi$belowLayers$VCroot_kmax[repl_vec,] <- xo$belowLayers$VCroot_kmax
      xi$belowLayers$RhizoPsi[repl_vec,] <- xo$belowLayers$RhizoPsi
    }
    xi$belowLayers$Wpool[repl_vec,] <- xo$belowLayers$Wpool
    xi$paramsPhenology[repl_vec,] <- xo$paramsPhenology
    xi$paramsAnatomy[repl_vec,] <- xo$paramsAnatomy
    xi$paramsInterception[repl_vec,] <- xo$paramsInterception
    xi$paramsTranspiration[repl_vec,] <- xo$paramsTranspiration
    xi$paramsWaterStorage[repl_vec,] <- xo$paramsWaterStorage
    xi$paramsGrowth[repl_vec,] <- xo$paramsGrowth
    xi$paramsAllometries[repl_vec,] <- xo$paramsAllometries
    xi$internalPhenology[repl_vec,] <- xo$internalPhenology
    xi$internalWater[repl_vec,] <- xo$internalWater
    xi$internalCarbon[repl_vec,] <- xo$internalCarbon
    xi$internalAllocation[repl_vec,] <- xo$internalAllocation
    xi$internalMortality[repl_vec,] <- xo$internalMortality
    xi$internalRings[repl_vec] <- xo$internalRings

    cat(paste0("   (c) Summaries\n"))
    
    # Store current forest state (after recruitment)
    forestStructures[[iYear+1]] = forest
    
    # Process summaries (after recruitment)
    cohSumYear = summarizeCohorts(iYear, xi)
    speciesSummary = rbind(speciesSummary,  summarizeSpecies(iYear,cohSumYear, xi))
    standSummary = rbind(standSummary,  summarizeStand(iYear,cohSumYear, xi))
    cohortSummary = rbind(cohortSummary, cohSumYear)
    
    # Update tree/shrub tables (after recruitment)
    treeTable = rbind(treeTable, createTreeTable(iYear, year, xi))
    shrubTable = rbind(shrubTable, createShrubTable(iYear, year, xi))
    
  }
  res = list(
    "StandSummary" = standSummary,
    "SpeciesSummary" = speciesSummary,
    "CohortSummary" = cohortSummary,
    "TreeTable" = treeTable,
    "ShrubTable" = shrubTable,
    "ForestStructures" = forestStructures,
    "GrowthResults" = growthResults)
  class(res)<-c("forestdynamics", "list")
  return(res)
}