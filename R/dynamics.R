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
                               "LeafAreaIndex" = x$above$LAI_live,
                               "DensityLive" = x$above$N,
                               "DensityDead" = x$internalMortality$N_dead,
                               "TreeBasalAreaLive"= x$above$N*pi*(x$above$DBH/200)^2,
                               "TreeBasalAreaDead"= x$internalMortality$N_dead*pi*(x$above$DBH/200)^2,
                               "ShrubCoverLive"= x$above$Cover,
                               "ShrubCoverDead"= x$above$Cover*(x$internalMortality$N_dead/x$above$N))
    return(cohortSummary)
  }
  summarizeSpecies<-function(step, cohSum, x) {
    lai_sp = tapply(cohSum$LeafAreaIndex, cohSum$Species, sum, na.rm=TRUE)
    nl_sp = tapply(cohSum$DensityLive, cohSum$Species, sum, na.rm=TRUE)
    nd_sp = tapply(cohSum$DensityDead, cohSum$Species, sum, na.rm=TRUE)
    bal_sp = tapply(cohSum$TreeBasalAreaLive, cohSum$Species, sum, na.rm=TRUE)
    bad_sp = tapply(cohSum$TreeBasalAreaDead, cohSum$Species, sum, na.rm=TRUE)
    shl_sp = tapply(cohSum$ShrubCoverLive, cohSum$Species, sum, na.rm=TRUE)
    shd_sp = tapply(cohSum$ShrubCoverDead, cohSum$Species, sum, na.rm=TRUE)
    mh_sp = tapply(x$above$H, x$above$SP, max, na.rm=TRUE)
    spSumYear <-data.frame("Step" = rep(step, length(lai_sp)),
                           "Species" = as.numeric(names(lai_sp)),
                           "LeafAreaIndex" = as.numeric(lai_sp),
                           "DensityLive"= as.numeric(nl_sp),
                           "DensityDead"= as.numeric(nd_sp),
                           "TreeBasalAreaLive"= as.numeric(bal_sp),
                           "TreeBasalAreaDead"= as.numeric(bad_sp),
                           "ShrubCoverLive"= as.numeric(shl_sp),
                           "ShrubCoverDead"= as.numeric(shd_sp),
                           "MaxHeight"= as.numeric(mh_sp))
  }
  summarizeStand<-function(step, cohSum, x) {
    standSumYear = data.frame("Step" = step,
                              "LeafAreaIndex" = sum(cohSum$LeafAreaIndex, na.rm = TRUE),
                              "DensityLive"= sum(cohSum$DensityLive, na.rm = TRUE),
                              "DensityDead"= sum(cohSum$DensityDead, na.rm = TRUE),
                              "TreeBasalAreaLive"= sum(cohSum$ShrubCoverLive, na.rm = TRUE),
                              "TreeBasalAreaDead"= sum(cohSum$ShrubCoverDead, na.rm = TRUE),
                              "ShrubCoverLive"= sum(cohSum$TreeBasalAreaLive, na.rm=TRUE),
                              "ShrubCoverDead"= sum(cohSum$TreeBasalAreaDead, na.rm=TRUE),
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
    return(tt)
  }
  createShrubTable<-function(step, year, x) {
    isShrub = !is.na(x$above$Cover)
    nt = sum(!isShrub)
    numCohorts = length(isShrub)
    range = numeric(0)
    if(numCohorts>nt) range = (nt+1):numCohorts
    tt<-data.frame(Step = step, Year = year,
                   Cohort = row.names(x$cohorts)[range],
                   Species = x$above$SP[range],
                   Name = x$cohorts$Name[range],
                   Cover = x$above$Cover[range],
                   Height = x$above$H[range])
    return(tt)
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
    cat(paste0("Growth simulation for year ", year, " (", iYear,"/", nYears,")\n"))
    meteoYear = meteo[years==year,]
    Gi = growth(xi, meteoYear, latitude = latitude, elevation = elevation, slope = slope, aspect = aspect)
    # Modified growth output
    xo = Gi$growthInput
    
    # Store growth results
    growthResults[[iYear]] = Gi
    # Process summaries
    cohSumYear = summarizeCohorts(iYear, xo)
    speciesSummary = rbind(speciesSummary,  summarizeSpecies(iYear,cohSumYear, xo))
    standSummary = rbind(standSummary,  summarizeStand(iYear,cohSumYear, xo))
    cohortSummary = rbind(cohortSummary, cohSumYear)
    # Update tree/shrub tables
    treeTable = rbind(treeTable, createTreeTable(iYear, year, xo))
    shrubTable = rbind(shrubTable, createShrubTable(iYear, year, xo))
    
    # Update forest state for next year
    isTree = is.na(xo$above$Cover)
    forest$treeData$N  = xo$above$N[isTree]
    forest$treeData$DBH  = xo$above$DBH[isTree]
    forest$treeData$Height  = xo$above$H[isTree]
    forest$shrubData$Cover  = xo$above$Cover[!isTree]
    forest$shrubData$Height  = xo$above$H[!isTree]
    
    # Store forest state
    forestStructures[[iYear+1]] = forest
    # Simulate recruitment
    recr_forest = forest
    recr_forest$treeData$N = 100
    recr_forest$treeData$DBH = 1
    recr_forest$treeData$Height = 100
    recr_forest$shrubData$Cover = 1
    recr_forest$shrubData$Height = 100
    recr_above = forest2aboveground(recr_forest, SpParams, NA, "MED")
    row.names(recr_above) = plant_ID(recr_forest, treeOffset, shrubOffset)
    treeOffset = treeOffset + nrow(recr_forest$treeData)
    shrubOffset = shrubOffset + nrow(recr_forest$shrubData)
    forest_above = forest2aboveground(forest, SpParams, NA, "MED")
    # merge in forest
    forest$treeData = rbind(forest$treeData, recr_forest$treeData)
    forest$shrubData = rbind(forest$shrubData, recr_forest$shrubData)
    # Prepare input for next year
    xi = growthInput(above = rbind(forest_above, recr_above),
                     Z50 = c(forest$treeData$Z50, forest$shrubData$Z50),
                     Z95 = c(forest$treeData$Z95, forest$shrubData$Z95),
                     soil, SpParams, control)
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