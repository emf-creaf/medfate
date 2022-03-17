recruitment<-function(forest, SpParams, control,
                      minMonthTemp, moistureIndex) {
  if((nrow(forest$treeData)>0) || (nrow(forest$shrubData)>0)) {
    PARperc = light_PARground(forest, SpParams)
  } else {
    PARperc = 100
  } 
  
  treeSpp = numeric(0)
  shrubSpp = numeric(0)
  if(is.null(control$seedRain)) {
    treeSpp = forest$treeData$Species
    if(length(treeSpp)>0) {
      sph = species_parameter(treeSpp, SpParams, "SeedProductionHeight")
      sph[is.na(sph)] = control$seedProductionTreeHeight
      treeSpp = treeSpp[forest$treeData$Height > sph]
      treeSpp = unique(treeSpp)
    }
    if(control$shrubDynamics) {
      shrubSpp = forest$shrubData$Species
      if(length(shrubSpp)>0) {
        sph = species_parameter(shrubSpp, SpParams, "SeedProductionHeight")
        sph[is.na(sph)] = control$seedProductionShrubHeight
        shrubSpp = shrubSpp[forest$shrubData$Height > sph]
        shrubSpp = unique(shrubSpp)
      }
    } 
  } else {
    isTree = !(species_characterParameter(control$seedRain, SpParams, "GrowthForm")=="Shrub")
    treeSpp = control$seedRain[isTree]
    if(control$shrubDynamics) shrubSpp = control$seedRain[!isTree]
  }
  recr_forest = emptyforest(ntree = length(treeSpp), nshrub=length(shrubSpp))
  
  ## Determine if species can recruit 
  tree_recr_selection = logical(0)
  tree_minFPAR = numeric(0)
  if(length(treeSpp)>0) {
    recr_forest$treeData$Species = treeSpp
    recr_forest$treeData$N = species_parameter(treeSpp, SpParams, "RecrTreeDensity")
    recr_forest$treeData$N[is.na(recr_forest$treeData$N)] = control$recrTreeDensity
    recr_forest$treeData$DBH = species_parameter(treeSpp, SpParams, "RecrTreeDBH")
    recr_forest$treeData$DBH[is.na(recr_forest$treeData$DBH)] = control$recrTreeDBH
    recr_forest$treeData$Height = species_parameter(treeSpp, SpParams, "RecrTreeHeight")
    recr_forest$treeData$Height[is.na(recr_forest$treeData$Height)] = control$recrTreeHeight
    recr_forest$treeData$Z50 = species_parameter(treeSpp, SpParams, "RecrZ50")
    recr_forest$treeData$Z50[is.na(recr_forest$treeData$Z50)] = control$recrTreeZ50
    recr_forest$treeData$Z95 = species_parameter(treeSpp, SpParams, "RecrZ95")
    recr_forest$treeData$Z95[is.na(recr_forest$treeData$Z95)] = control$recrTreeZ95
    minTemp = species_parameter(treeSpp, SpParams, "MinTempRecr")
    minMoisture = species_parameter(treeSpp, SpParams, "MinMoistureRecr")
    minTemp[is.na(minTemp)] = control$minTempRecr
    minMoisture[is.na(minMoisture)] = control$minMoistureRecr
    tree_minFPAR = species_parameter(treeSpp, SpParams, "MinFPARRecr")
    tree_minFPAR[is.na(tree_minFPAR)] = control$minFPARRecr
    tree_recr_selection = (minMonthTemp > minTemp) & (moistureIndex > minMoisture) & (PARperc > tree_minFPAR)
  }
  shrub_recr_selection = logical(0)
  shrub_minFPAR = numeric(0)
  if((length(shrubSpp)>0)) {
    recr_forest$shrubData$Species = shrubSpp
    recr_forest$shrubData$Cover = species_parameter(shrubSpp, SpParams, "RecrShrubCover")
    recr_forest$shrubData$Cover[is.na(recr_forest$shrubData$Cover)] = control$recrShrubCover
    recr_forest$shrubData$Height = species_parameter(shrubSpp, SpParams, "RecrShrubHeight")
    recr_forest$shrubData$Height[is.na(recr_forest$shrubData$Height)] = control$recrShrubHeight
    recr_forest$shrubData$Z50 = species_parameter(shrubSpp, SpParams, "RecrZ50")
    recr_forest$shrubData$Z50[is.na(recr_forest$shrubData$Z50)] = control$recrShrubZ50
    recr_forest$shrubData$Z95 = species_parameter(shrubSpp, SpParams, "RecrZ95")
    recr_forest$shrubData$Z95[is.na(recr_forest$shrubData$Z95)] = control$recrShrubZ95
    minTemp = species_parameter(shrubSpp, SpParams, "MinTempRecr")
    minMoisture = species_parameter(shrubSpp, SpParams, "MinMoistureRecr")
    minTemp[is.na(minTemp)] = control$minTempRecr
    minMoisture[is.na(minMoisture)] = control$minMoistureRecr
    shrub_minFPAR = species_parameter(shrubSpp, SpParams, "MinFPARRecr")
    shrub_minFPAR[is.na(shrub_minFPAR)] = control$minFPARRecr
    shrub_recr_selection = (minMonthTemp > minTemp) & (moistureIndex > minMoisture) & (PARperc > shrub_minFPAR)
    if(!control$shrubDynamics) shrub_recr_selection = rep(FALSE, nrow(recr_forest$shrubData))
  }

  tree_recr_prob = species_parameter(treeSpp, SpParams, "ProbRecr")
  tree_recr_prob[is.na(tree_recr_prob)] = control$probRecr
  tree_recr_prob[!tree_recr_selection] = 0
  shrub_recr_prob = species_parameter(shrubSpp, SpParams, "ProbRecr")
  shrub_recr_prob[is.na(shrub_recr_prob)] = control$probRecr
  shrub_recr_prob[!shrub_recr_selection] = 0
  
  if(control$recruitmentMode=="stochastic") {
    recr_forest$treeData$N = rbinom(length(treeSpp), size = 1, prob = tree_recr_prob)*rpois(length(treeSpp), recr_forest$treeData$N)
    recr_forest$shrubData$Cover = rbinom(length(shrubSpp), size = 1, prob = shrub_recr_prob)*rpois(length(shrubSpp), recr_forest$shrubData$Cover)
  } else {
    recr_forest$treeData$N = recr_forest$treeData$N*tree_recr_prob
    recr_forest$shrubData$Cover = recr_forest$shrubData$Cover*shrub_recr_prob
  }
  recr_forest$treeData = recr_forest$treeData[recr_forest$treeData$N>0, ,drop=FALSE]
  recr_forest$shrubData = recr_forest$shrubData[recr_forest$shrubData$Cover>0, ,drop=FALSE]
  return(recr_forest)
}

fordyn<-function(forest, soil, SpParams,
                 meteo, control,
                 latitude , elevation = NA, slope = NA, aspect = NA,
                 management_function = NULL, management_args = NULL) {
  
  # Modify control parameters
  verboseDyn = control$verbose
  control$verbose = FALSE
  control$subdailyResults = FALSE
  
  dates = as.Date(row.names(meteo))
  years = as.numeric(format(dates, "%Y"))
  months = as.numeric(format(dates, "%m"))
  yearsUnique = unique(years)
  nYears = length(yearsUnique)
  growthResults = vector("list", nYears)
  names(growthResults) = paste0("Step_", 1:nYears)
  forestStructures = vector("list", nYears+1)
  names(forestStructures) = c("Initial", paste0("Step_", 1:nYears))
  
  
  #Initialization
  if(inherits(forest, "fordyn")) {
    if(verboseDyn) cat(paste0("Initialisation from previous run\n"))
    xi = forest$NextInputObject
    forest = forest$NextForestObject
  } else {
    #Subset columns relevant for fordyn (in case there are other)
    if(control$allowRecruitment) {
      forest$treeData = forest$treeData[,c("Species","DBH", "Height","N","Z50","Z95")]
      forest$shrubData = forest$shrubData[,c("Species","Height","Cover", "Z50","Z95")]
    }
    xi = forest2growthInput(forest, soil, SpParams, control)
  }
  forestStructures[[1]] = forest
  treeOffset = nrow(forest$treeData)
  shrubOffset = nrow(forest$shrubData)

  #initial tree/shrub tables
  treeTable = .createTreeTable(0, NA, xi)
  shrubTable = .createShrubTable(0, NA, xi)
  deadTreeTable = .createDeadTreeTable(0, NA, xi)
  deadShrubTable = .createDeadShrubTable(0, NA, xi)
  cutTreeTable = treeTable[numeric(),,drop = FALSE]
  cutShrubTable = shrubTable[numeric(),,drop = FALSE]
  
  #initial summaries
  cohortSummary <-.summarizeCohorts(0, 
                                    treeTable, shrubTable,
                                    deadTreeTable, deadShrubTable,
                                    cutTreeTable, cutShrubTable)
  speciesSummary<-.summarizeSpecies(0,cohortSummary, xi, SpParams)
  standSummary<-.summarizeStand(0,cohortSummary, xi)

  
  #Simulations
  for(iYear in 1:nYears) {
    year = yearsUnique[iYear]
    if(verboseDyn) cat(paste0("Simulating year ", year, " (", iYear,"/", nYears,"): "))
    meteoYear = meteo[years==year,]
    monthsYear = months[years==year]
    # 1.1 Calls growth model
    if(verboseDyn) cat(paste0(" (a) Growth/mortality"))
    Gi = growth(xi, meteoYear, latitude = latitude, elevation = elevation, slope = slope, aspect = aspect)
    
    # 1.2 Store growth results
    growthResults[[iYear]] = Gi
    
    # 1.3 Retrieve modified growth output
    xo = Gi$growthOutput
    
    # 2.2 Update dead tree/shrub tables
    deadTreeTableYear = .createDeadTreeTable(iYear, year, xo)
    deadShrubTableYear = .createDeadShrubTable(iYear, year, xo)
    
    # 2.2 Update forest structural variables
    isTree = is.na(xo$above$Cover)
    forest$treeData$N  = xo$above$N[isTree]
    forest$treeData$DBH  = xo$above$DBH[isTree]
    forest$treeData$Height  = xo$above$H[isTree]
    if(control$shrubDynamics) {
      forest$shrubData$Cover  = xo$above$Cover[!isTree]
      forest$shrubData$Height  = xo$above$H[!isTree]
    }
    
    # 2.3 Call management function if required
    cutTreeTableYear = NULL
    cutShrubTableYear = NULL
    if(!is.null(management_function)) {
      res = do.call(management_function, list(x = forest, args= management_args, verbose = FALSE))
      if(verboseDyn) cat(paste0(" & management [", res$action,"]"))
      # Update forest and xo objects
      forest$treeData$N <- pmax(0,forest$treeData$N - res$N_tree_cut)
      xo$above$N[isTree] <- forest$treeData$N
      forest$shrubData$Cover <- pmax(0,forest$shrubData$Cover - res$Cover_shrub_cut)
      xo$above$Cover[!isTree] <- forest$shrubData$Cover
      # Update cut tables
      cutTreeTableYear = .createCutTreeTable(iYear, year, xo, res$N_tree_cut)
      cutShrubTableYear = .createCutShrubTable(iYear, year, xo, res$Cover_shrub_cut)
      # Retrieve plantation information
      planted_forest <- res$planted_forest
      if(nrow(planted_forest$treeData)>0) {
        for(i in 1:nrow(planted_forest$treeData)) {
          planted_forest$treeData$Z50[i] = species_parameter(planted_forest$treeData$Species[i], SpParams,"RecrZ50")
          planted_forest$treeData$Z95[i] = species_parameter(planted_forest$treeData$Species[i], SpParams,"RecrZ95")
          if(is.na(planted_forest$treeData$Z50[i])) planted_forest$treeData$Z50[i] = 250
          if(is.na(planted_forest$treeData$Z95[i])) planted_forest$treeData$Z95[i] = 500
        }
      }
      if(nrow(planted_forest$shrubData)>0) {
        for(i in 1:nrow(planted_forest$shrubData)) {
          planted_forest$shrubData$Z50[i] = species_parameter(planted_forest$shrubData$Species[i], SpParams,"RecrZ50")
          planted_forest$shrubData$Z95[i] = species_parameter(planted_forest$shrubData$Species[i], SpParams,"RecrZ95")
          if(is.na(planted_forest$shrubData$Z50[i])) planted_forest$shrubData$Z50[i] = 100
          if(is.na(planted_forest$shrubData$Z95[i])) planted_forest$shrubData$Z95[i] = 300
        }
      }
      
      # Store new management arguments (may have changed)
      management_args <- res$management_args
    } else {
      planted_forest = emptyforest()
    }
    # 2.4 Remove empty cohorts if required
    emptyTrees = rep(FALSE, nrow(forest$treeData))
    emptyShrubs = rep(FALSE, nrow(forest$shrubData))
    if(control$removeEmptyCohorts) {
      emptyTrees = (forest$treeData$N < control$minimumCohortDensity)
      if(control$shrubDynamics) emptyShrubs = (forest$shrubData$Cover < control$minimumCohortDensity)
    }
    emptyCohorts = c(emptyTrees, emptyShrubs)
    if(sum(emptyCohorts)>0) {
      forest$treeData = forest$treeData[!emptyTrees,, drop=FALSE] 
      forest$shrubData = forest$shrubData[!emptyShrubs,, drop=FALSE] 
      # Remove from growth input object
      xo$cohorts = xo$cohorts[!emptyCohorts, , drop=FALSE] 
      xo$above = xo$above[!emptyCohorts, , drop=FALSE] 
      xo$below = xo$below[!emptyCohorts, , drop=FALSE] 
      xo$belowLayers$V = xo$belowLayers$V[!emptyCohorts, , drop=FALSE] 
      xo$belowLayers$L = xo$belowLayers$L[!emptyCohorts, , drop=FALSE] 
      if(control$transpirationMode=="Sperry") {
        xo$belowLayers$VGrhizo_kmax = xo$belowLayers$VGrhizo_kmax[!emptyCohorts, , drop=FALSE]
        xo$belowLayers$VCroot_kmax = xo$belowLayers$VCroot_kmax[!emptyCohorts, , drop=FALSE]
        xo$belowLayers$RhizoPsi = xo$belowLayers$RhizoPsi[!emptyCohorts, , drop=FALSE]
      }
      xo$belowLayers$Wpool = xo$belowLayers$Wpool[!emptyCohorts, , drop=FALSE]
      xo$paramsPhenology = xo$paramsPhenology[!emptyCohorts,, drop=FALSE]
      xo$paramsAnatomy = xo$paramsAnatomy[!emptyCohorts,, drop=FALSE]
      xo$paramsInterception = xo$paramsInterception[!emptyCohorts,, drop=FALSE]
      xo$paramsTranspiration = xo$paramsTranspiration[!emptyCohorts,, drop=FALSE]
      xo$paramsWaterStorage = xo$paramsWaterStorage[!emptyCohorts,, drop=FALSE]
      xo$paramsGrowth = xo$paramsGrowth[!emptyCohorts,, drop=FALSE]
      xo$paramsAllometries = xo$paramsAllometries[!emptyCohorts,, drop=FALSE]
      xo$internalPhenology = xo$internalPhenology[!emptyCohorts,, drop=FALSE]
      xo$internalWater = xo$internalWater[!emptyCohorts,, drop=FALSE]
      xo$internalCarbon = xo$internalCarbon[!emptyCohorts,, drop=FALSE]
      xo$internalAllocation = xo$internalAllocation[!emptyCohorts,, drop=FALSE]
      xo$internalMortality = xo$internalMortality[!emptyCohorts,, drop=FALSE]
      xo$internalRings = xo$internalRings[!emptyCohorts]
    }

    
    # 3. Simulate species recruitment
    if(control$allowRecruitment) {
      if(verboseDyn) cat(paste0(", (b) Recruitment"))
      monthlyTemp = tapply(meteoYear$MeanTemperature, monthsYear, FUN="mean", na.rm=TRUE)
      minMonthTemp = min(monthlyTemp, na.rm=TRUE)
      moistureIndex = sum(meteoYear$Precipitation, na.rm=TRUE)/sum(meteoYear$PET, na.rm=TRUE)
      # if(verboseDyn) cat(paste0("       Coldest month mean temp. (Celsius): ", round(minMonthTemp,2), "   Moisture index: ", round(moistureIndex,2), "   FPAR (%): ", round(PARperc,1), "\n"))
      recr_forest = recruitment(forest, SpParams, control, minMonthTemp, moistureIndex)
    } else {
      recr_forest = emptyforest()
    }

    
    # 4.1 Generate above-ground data
    planted_above = forest2aboveground(planted_forest, SpParams, NA, "MED")
    row.names(planted_above) = plant_ID(planted_forest, treeOffset, shrubOffset)
    treeOffset = treeOffset + nrow(planted_forest$treeData)
    shrubOffset = shrubOffset + nrow(planted_forest$shrubData)
    recr_above = forest2aboveground(recr_forest, SpParams, NA, "MED")
    row.names(recr_above) = plant_ID(recr_forest, treeOffset, shrubOffset)
    treeOffset = treeOffset + nrow(recr_forest$treeData)
    shrubOffset = shrubOffset + nrow(recr_forest$shrubData)
    forest_above = forest2aboveground(forest, SpParams, NA, "MED")
    row.names(forest_above) = row.names(xo$cohorts)
    forest_above$LAI_live[!is.na(forest_above$DBH)] = xo$above$LAI_live[!is.na(forest_above$DBH)]
    forest_above$LAI_expanded[!is.na(forest_above$DBH)] = xo$above$LAI_expanded[!is.na(forest_above$DBH)]
    if(control$shrubDynamics) {
      forest_above$LAI_live[is.na(forest_above$DBH)] = xo$above$LAI_live[is.na(forest_above$DBH)]
      forest_above$LAI_expanded[is.na(forest_above$DBH)] = xo$above$LAI_expanded[is.na(forest_above$DBH)]
    }

    # 4.2 Merge above-ground data (first trees)
    above_all = rbind(forest_above[!is.na(forest_above$DBH),, drop = FALSE], 
                      planted_above[!is.na(planted_above$DBH),, drop = FALSE],
                      recr_above[!is.na(recr_above$DBH),, drop = FALSE],
                      forest_above[is.na(forest_above$DBH),, drop = FALSE],
                      planted_above[is.na(planted_above$DBH),, drop = FALSE],
                      recr_above[is.na(recr_above$DBH),, drop = FALSE])
    
    # 4.3 Logical vector for replacement
    repl_vec <- c(rep(TRUE, nrow(forest$treeData)),
                  rep(FALSE, nrow(planted_forest$treeData)),
                  rep(FALSE, nrow(recr_forest$treeData)),
                  rep(control$shrubDynamics, nrow(forest$shrubData)),
                  rep(FALSE, nrow(planted_forest$shrubData)),
                  rep(FALSE, nrow(recr_forest$shrubData)))
    sel_vec = c(rep(TRUE, nrow(forest$treeData)),
                rep(control$shrubDynamics, nrow(forest$shrubData)))
    
    # 4.4 Merge cohorts in forest object
    forest$treeData = rbind(forest$treeData, planted_forest$treeData, recr_forest$treeData)
    forest$shrubData = rbind(forest$shrubData, planted_forest$shrubData, recr_forest$shrubData)
    
    
    # 4.5 Prepare growth input for next year
    xi = growthInput(above = above_all,
                     Z50 = c(forest$treeData$Z50, forest$shrubData$Z50),
                     Z95 = c(forest$treeData$Z95, forest$shrubData$Z95),
                     xo$soil, SpParams, control)
    
    # 4.6 Replace previous state for surviving cohorts
    xi$cohorts[repl_vec,] <- xo$cohorts[sel_vec,, drop=FALSE]
    xi$above[repl_vec,] <- xo$above[sel_vec,, drop=FALSE]
    xi$below[repl_vec,] <- xo$below[sel_vec,, drop=FALSE]
    xi$belowLayers$V[repl_vec,] <- xo$belowLayers$V[sel_vec,, drop=FALSE]
    xi$belowLayers$L[repl_vec,] <- xo$belowLayers$L[sel_vec,, drop=FALSE]
    if(control$transpirationMode=="Sperry") {
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
    xi$internalRings[repl_vec] <- xo$internalRings[sel_vec]

    
    # 5.1 Store current forest state (after recruitment)
    forestStructures[[iYear+1]] = forest
    
    # 5.2 Process summaries (after recruitment)
    treeTableYear = .createTreeTable(iYear, year, xi)
    shrubTableYear = .createShrubTable(iYear, year, xi)
    cohSumYear = .summarizeCohorts(iYear,
                                   treeTableYear, shrubTableYear,
                                   deadTreeTableYear, deadShrubTableYear,
                                   cutTreeTableYear, cutShrubTableYear)
    speciesSummary = rbind(speciesSummary,  .summarizeSpecies(iYear,cohSumYear, xi, SpParams))
    standSummary = rbind(standSummary,  .summarizeStand(iYear,cohSumYear, xi))
    cohortSummary = rbind(cohortSummary, cohSumYear)
    
    # 5.3 Update tree/shrub tables (after recruitment)
    treeTable = rbind(treeTable, treeTableYear)
    shrubTable = rbind(shrubTable, shrubTableYear)
    deadTreeTable = rbind(deadTreeTable, deadTreeTableYear)
    deadShrubTable = rbind(deadShrubTable, deadShrubTableYear)
    if(!is.null(cutTreeTableYear)) cutTreeTable = rbind(cutTreeTable, cutTreeTableYear)
    if(!is.null(cutShrubTableYear)) cutShrubTable = rbind(cutShrubTable, cutShrubTableYear)
    
    if(verboseDyn) cat(paste0("\n"))
  }
  res = list(
    "StandSummary" = standSummary,
    "SpeciesSummary" = speciesSummary,
    "CohortSummary" = cohortSummary,
    "TreeTable" = treeTable,
    "DeadTreeTable" = deadTreeTable,
    "CutTreeTable" = cutTreeTable,
    "ShrubTable" = shrubTable,
    "DeadShrubTable" = deadShrubTable,
    "CutShrubTable" = cutShrubTable,
    "ForestStructures" = forestStructures,
    "GrowthResults" = growthResults,
    "ManagementArgs" = management_args,
    "NextInputObject" = xi,
    "NextForestObject" = forest)
  class(res)<-c("fordyn", "list")
  return(res)
}