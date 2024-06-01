#' Plant regeneration
#'
#' Functions to simulate annual plant regeneration from seed recruitment or from resprouting
#' 
#' @param forest An object of class \code{\link{forest}}.
#' @param seedBank A data frame with columns 'Species' and 'Percent', describing a seed bank.
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsDefinition}}).
#' @param refillSpecies A string vector of species names corresponding to seed rain to refill seed bank.
#' @param refillPercent A numeric vector of indicating the percentage of seed bank refilling (if missing then seed bank is set to 100\%).
#' @param minPercent A minimum percent of seed bank to retain entry in \code{seedBank} element of \code{forest}.
#' @param control A list with default control parameters (see \code{\link{defaultControl}}).
#' @param minMonthTemp Minimum month temperature.
#' @param moistureIndex Moisture index (annual precipitation over annual potential evapotranspiration).
#' @param verbose Boolean flag to indicate console output during calculations.
#' 
#' @details 
#' \itemize{
#'   \item{\code{regeneration_seedproduction} evaluates if reproductive individuals (i.e. sufficiently tall individuals) are present.} 
#'   \item{\code{regeneration_seedrefill} fills seed bank of input \code{forest} object with seed rain.}
#'   \item{\code{regeneration_seedmortality} updates seed bank of input \code{forest} object according to annual seed mortality.}
#'   \item{\code{regeneration_recruitment} evaluates recruitment from the seed bank (or local seed production if seed bank is missing). Minimum month temperature and moisture index values are used to determine if recruitment was successful. 
#' Species also require a minimum amount of light at the ground level.}
#'   \item{\code{regeneration_resprouting} evaluates resprouting occurs after “mortality” from die-back (including drought- or pathogen-induced dessication), 
#' cutting or burning of the aerial part in a species with resprouting ability, 
#' but not after carbon starvation or baseline mortality (unspecific mortality causes).}
#' }
#' 
#' @return 
#' \itemize{
#'   \item{\code{regeneration_seedproduction} returns a list of species names}
#'   \item{\code{regeneration_seedrefill} and \code{regeneration_seedmortality} return a copy of the input \code{data.frame} object with an update seed bank. }
#'   \item{\code{regeneration_resprouting} and \code{regeneration_recruitment} return a new object of class \code{\link{forest}} with the new plant cohorts.}
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{fordyn}}
#' 
#' @examples 
#' #Load example plot plant data
#' data(exampleforest)
#' 
#' #Default species parameterization
#' data(SpParamsMED)
#' 
#' #Initialize control parameters
#' control <- defaultControl("Granier")
#' control$recruitmentMode = "deterministic" 
#' 
#' #Recruitment limits
#' plant_parameter(exampleforest, SpParamsMED, "MinTempRecr")
#' plant_parameter(exampleforest, SpParamsMED, "MinMoistureRecr")
#' 
#' #Compare seed recruitment outcomes
#' regeneration_recruitment(exampleforest, SpParamsMED, control, 0, 0.25)
#' regeneration_recruitment(exampleforest, SpParamsMED, control, 3, 0.25)
#' 
#' @name regeneration
#' @keywords internal
regeneration_seedproduction<-function(forest, SpParams, control) {
  treeSpp <- forest$treeData$Species
  if(length(treeSpp)>0) {
    sph <- species_parameter(treeSpp, SpParams, "SeedProductionHeight")
    sph[is.na(sph)] <- control$seedProductionTreeHeight
    treeSpp <- treeSpp[forest$treeData$Height > sph]
    treeSpp <- unique(treeSpp)
  }
  if(control$shrubDynamics) {
    shrubSpp <- forest$shrubData$Species
    if(length(shrubSpp)>0) {
      sph <- species_parameter(shrubSpp, SpParams, "SeedProductionHeight")
      sph[is.na(sph)] <- control$seedProductionShrubHeight
      shrubSpp <- shrubSpp[forest$shrubData$Height > sph]
      shrubSpp <- unique(shrubSpp)
    }
  } else {
    shrubSpp <- character(0)
  } 
  return(unique(c(treeSpp, shrubSpp)))
}

#' @rdname regeneration
#' @keywords internal
regeneration_seedrefill <-function(seedBank, refillSpecies, refillPercent = NULL) {
  # Initialize seed bank if not present
  if(is.null(seedBank)) {
    seedBank <- data.frame(Species = character(0), Percent = numeric(0))
  }
  
  #If not supplied, all species refill to 100%
  if(is.null(refillPercent)) refillPercent <- rep(100, length(refillSpecies))
  
  #If refillPercent is missing for some species, it will be refilled to 100%
  refillPercent[is.na(refillPercent)] <- 100
  
  if(length(refillSpecies) > 0) {
    for(isp in 1:length(refillSpecies)) {
      sp <- refillSpecies[isp]
      if(sp %in% seedBank$Species) {
        seedBank$Percent[seedBank$Species==sp] <- min(100.0, seedBank$Percent[seedBank$Species==sp] + refillPercent[isp])
      } else {
        seedBank <- rbind(seedBank, data.frame("Species" = sp, "Percent" = refillPercent[isp]))
      }
    }
  }
  return(seedBank)
}

#' @rdname regeneration
#' @keywords internal
regeneration_seedmortality <- function(seedBank, SpParams, minPercent = 1) {
  if(is.null(seedBank)) {
    seedBank <- data.frame(Species = character(0), Percent = numeric(0))
    return(seedBank)
  }
  seed_longevity <- species_parameter(seedBank$Species, SpParams, "SeedLongevity")
  seed_longevity[is.na(seed_longevity)] <- 2.0 # Default 2 years
  seedBank$Percent <- seedBank$Percent*exp(-1.0/seed_longevity)
  seedBank <- seedBank[seedBank$Percent >= minPercent, , drop = FALSE]
  return(seedBank)
}

#' @rdname regeneration
#' @keywords internal
regeneration_recruitment<-function(forest, SpParams, control,
                      minMonthTemp, moistureIndex, verbose = FALSE) {
  if((nrow(forest$treeData)>0) || (nrow(forest$shrubData)>0)) {
    PARperc <- light_PARground(forest, SpParams)
  } else {
    PARperc <- 100
  } 
  if(verbose) cat(paste0(" [Cold (C): ", round(minMonthTemp,2), " / Moist: ", round(moistureIndex,2), " / FPAR (%): ", round(PARperc,1), "]"))
  
  # If seed bank is defined take seeds from it. Otherwise, use local seed production
  if("seedBank" %in% names(forest)) {
    seedSpp <- forest$seedBank$Species 
    seedPercent <- forest$seedBank$Percent
  } else {
    seedSpp <- regeneration_seedproduction(forest, SpParams, control)
    seedPercent <- rep(100.0, length(seedSpp))
  }
  isTree <- !(species_characterParameter(seedSpp, SpParams, "GrowthForm")=="Shrub")
  treeSpp <- seedSpp[isTree]
  treePercent <- seedPercent[isTree]
  if(control$shrubDynamics) {
    shrubSpp <- seedSpp[!isTree]
    shrubPercent <- seedPercent[!isTree]
  } else {
    shrubSpp <- character(0)
    shrubPercent <- numeric(0)
  }
  recr_forest <- emptyforest(ntree = length(treeSpp), nshrub=length(shrubSpp))
  ## Determine if species can recruit 
  tree_recr_selection <- logical(0)
  tree_minFPAR <- numeric(0)
  if(length(treeSpp)>0) {
    recr_forest$treeData$Species <- treeSpp
    recr_forest$treeData$N <- species_parameter(treeSpp, SpParams, "RecrTreeDensity")
    recr_forest$treeData$N[is.na(recr_forest$treeData$N)] <- control$recrTreeDensity
    recr_forest$treeData$N <- recr_forest$treeData$N*(treePercent/100.0) #Apply reduction due to seed bank size
    recr_forest$treeData$DBH <- species_parameter(treeSpp, SpParams, "RecrTreeDBH")
    recr_forest$treeData$DBH[is.na(recr_forest$treeData$DBH)] <- control$recrTreeDBH
    recr_forest$treeData$Height <- species_parameter(treeSpp, SpParams, "RecrTreeHeight")
    recr_forest$treeData$Height[is.na(recr_forest$treeData$Height)] <- control$recrTreeHeight
    recr_forest$treeData$Z50 <- species_parameter(treeSpp, SpParams, "RecrZ50")
    recr_forest$treeData$Z50[is.na(recr_forest$treeData$Z50)] <- control$recrTreeZ50
    recr_forest$treeData$Z95 <- species_parameter(treeSpp, SpParams, "RecrZ95")
    recr_forest$treeData$Z95[is.na(recr_forest$treeData$Z95)] <- control$recrTreeZ95
    minTemp <- species_parameter(treeSpp, SpParams, "MinTempRecr")
    minMoisture <- species_parameter(treeSpp, SpParams, "MinMoistureRecr")
    minTemp[is.na(minTemp)] <- control$minTempRecr
    minMoisture[is.na(minMoisture)] <- control$minMoistureRecr
    tree_minFPAR <- species_parameter(treeSpp, SpParams, "MinFPARRecr")
    tree_minFPAR[is.na(tree_minFPAR)] <- control$minFPARRecr
    tree_recr_selection <- (minMonthTemp > minTemp) & (moistureIndex > minMoisture) & (PARperc > tree_minFPAR)
    # if(verbose) print(cbind(treeSpp, tree_recr_selection))
  }
  shrub_recr_selection <- logical(0)
  shrub_minFPAR <- numeric(0)
  if((length(shrubSpp)>0)) {
    recr_forest$shrubData$Species <- shrubSpp
    recr_forest$shrubData$Cover <- species_parameter(shrubSpp, SpParams, "RecrShrubCover")
    recr_forest$shrubData$Cover[is.na(recr_forest$shrubData$Cover)] <- control$recrShrubCover
    recr_forest$shrubData$Cover <- recr_forest$shrubData$Cover*(shrubPercent/100.0) #Apply reduction due to seed bank size
    recr_forest$shrubData$Height <- species_parameter(shrubSpp, SpParams, "RecrShrubHeight")
    recr_forest$shrubData$Height[is.na(recr_forest$shrubData$Height)] <- control$recrShrubHeight
    recr_forest$shrubData$Z50 <- species_parameter(shrubSpp, SpParams, "RecrZ50")
    recr_forest$shrubData$Z50[is.na(recr_forest$shrubData$Z50)] <- control$recrShrubZ50
    recr_forest$shrubData$Z95 <- species_parameter(shrubSpp, SpParams, "RecrZ95")
    recr_forest$shrubData$Z95[is.na(recr_forest$shrubData$Z95)] <- control$recrShrubZ95
    minTemp <- species_parameter(shrubSpp, SpParams, "MinTempRecr")
    minMoisture <- species_parameter(shrubSpp, SpParams, "MinMoistureRecr")
    minTemp[is.na(minTemp)] <- control$minTempRecr
    minMoisture[is.na(minMoisture)] <- control$minMoistureRecr
    shrub_minFPAR <- species_parameter(shrubSpp, SpParams, "MinFPARRecr")
    shrub_minFPAR[is.na(shrub_minFPAR)] <- control$minFPARRecr
    shrub_recr_selection <- (minMonthTemp > minTemp) & (moistureIndex > minMoisture) & (PARperc > shrub_minFPAR)
    if(!control$shrubDynamics) shrub_recr_selection <- rep(FALSE, nrow(recr_forest$shrubData))
  }
  
  tree_recr_prob <- species_parameter(treeSpp, SpParams, "ProbRecr")
  tree_recr_prob[is.na(tree_recr_prob)] <- control$probRecr
  tree_recr_prob[!tree_recr_selection] <- 0
  shrub_recr_prob <- species_parameter(shrubSpp, SpParams, "ProbRecr")
  shrub_recr_prob[is.na(shrub_recr_prob)] <- control$probRecr
  shrub_recr_prob[!shrub_recr_selection] <- 0
  if(control$recruitmentMode=="stochastic") {
    recr_forest$treeData$N <- rbinom(length(treeSpp), size = 1, prob = tree_recr_prob)*rpois(length(treeSpp), recr_forest$treeData$N)
    recr_forest$shrubData$Cover <- rbinom(length(shrubSpp), size = 1, prob = shrub_recr_prob)*rpois(length(shrubSpp), recr_forest$shrubData$Cover)
  } else {
    recr_forest$treeData$N <- recr_forest$treeData$N*tree_recr_prob
    recr_forest$shrubData$Cover <- recr_forest$shrubData$Cover*shrub_recr_prob
  }
  recr_forest$treeData <- recr_forest$treeData[recr_forest$treeData$N>0, ,drop=FALSE]
  recr_forest$shrubData <- recr_forest$shrubData[recr_forest$shrubData$Cover>0, ,drop=FALSE]
  
  if("herbCover" %in% names(recr_forest)) recr_forest$herbCover <- NULL
  if("herbHeight" %in% names(recr_forest)) recr_forest$herbHeight <- NULL
  if("seedBank" %in% names(recr_forest)) recr_forest$seedBank <- NULL
  
  return(recr_forest)
}


#' @rdname regeneration
#' 
#' @param internalMortality A data frame with mortality occurred in the last year of simulation. 
#' @param management_results The result of calling a management function (see \code{\link{defaultManagementFunction}}).
#' @keywords internal
regeneration_resprouting <- function(forest, internalMortality, SpParams, control, 
                        management_results = NULL) {
  n_trees <- nrow(forest$treeData)
  n_shrubs <- nrow(forest$shrubData)
  
  ## Resprouting survivorship
  resp_dist_trees <- species_parameter(forest$treeData$Species, SpParams, "RespDist")
  resp_dist_trees[is.na(resp_dist_trees)] <- 0
  resp_dist_shrubs <- species_parameter(forest$shrubData$Species, SpParams, "RespDist")
  resp_dist_shrubs[is.na(resp_dist_shrubs)] <- 0
  resp_fire_trees <- species_parameter(forest$treeData$Species, SpParams, "RespFire")
  resp_fire_trees[is.na(resp_fire_trees)] <- 0
  resp_fire_shrubs <- species_parameter(forest$shrubData$Species, SpParams, "RespFire")
  resp_fire_shrubs[is.na(resp_fire_shrubs)] <- 0
  resp_clip_trees <- species_parameter(forest$treeData$Species, SpParams, "RespClip")
  resp_clip_trees[is.na(resp_clip_trees)] <- 0
  resp_clip_shrubs <- species_parameter(forest$shrubData$Species, SpParams, "RespClip")
  resp_clip_shrubs[is.na(resp_clip_shrubs)] <- 0
  
  N_resprouting_dessication <- internalMortality$N_dessication
  N_resprouting_burnt <- internalMortality$N_burnt
  if(n_trees>0) {
    # Resprouting survivorship (fire + dessication)
    N_surv_dessication <- N_resprouting_dessication[1:n_trees]*resp_dist_trees
    N_surv_fire <- N_resprouting_burnt[1:n_trees]*resp_fire_trees
    N_surv <- N_surv_dessication + N_surv_fire
    # Resprouting vigor depends on cm2 of stump area
    N_resprouting <- N_surv*pi*(forest$treeData$DBH/2)^2*1.82*10^(-0.053*5) 
  } else { 
    N_resprouting <- numeric(0)
  }
  Cover_resprouting_dessication <- internalMortality$Cover_dessication
  Cover_resprouting_burnt <- internalMortality$Cover_burnt
  if(n_shrubs>0) {
    Cover_resprouting <- Cover_resprouting_dessication[(n_trees+1):(n_trees+n_shrubs)]*resp_dist_shrubs
    Cover_resprouting <- Cover_resprouting + Cover_resprouting_burnt[(n_trees+1):(n_trees+n_shrubs)]*resp_fire_shrubs
  } else {
    Cover_resprouting <- numeric(0)
  }
  if(!is.null(management_results)) { ## Add tree/shrub cuts
    N_surv_cut <- management_results$N_tree_cut*resp_clip_trees
    N_resprouting_cut <- N_surv_cut*pi*(forest$treeData$DBH/2)^2*1.82*10^(-0.053*5) 
    N_resprouting <- N_resprouting + N_resprouting_cut
    Cover_resprouting <- Cover_resprouting + management_results$Cover_shrub_cut*resp_clip_shrubs
  }
  
  # Copy forest to inherit species and belowground information
  resp_forest <- forest 
  resp_forest$treeData$N <- pmin(control$recrTreeDensity, N_resprouting) #Truncates maximum density to control variable
  
  resp_forest$treeData$DBH <- species_parameter(resp_forest$treeData$Species, SpParams, "RecrTreeDBH")
  resp_forest$treeData$DBH[is.na(resp_forest$treeData$DBH)] <- control$recrTreeDBH
  resp_forest$treeData$Height <- species_parameter(resp_forest$treeData$Species, SpParams, "RecrTreeHeight")
  resp_forest$treeData$Height[is.na(resp_forest$treeData$Height)] <- control$recrTreeHeight
  
  resp_forest$shrubData$Cover <- Cover_resprouting
  resp_forest$shrubData$Height <- species_parameter(resp_forest$shrubData$Species, SpParams, "RecrShrubHeight")
  resp_forest$shrubData$Height[is.na(resp_forest$shrubData$Height)] <- control$recrShrubHeight
  
  # Trim species with no resprouting
  resp_forest$treeData <- resp_forest$treeData[resp_forest$treeData$N > 0, , drop = FALSE]
  resp_forest$shrubData <- resp_forest$shrubData[resp_forest$shrubData$Cover > 0, , drop = FALSE]
  
  if("herbCover" %in% names(resp_forest)) resp_forest$herbCover <- NULL
  if("herbHeight" %in% names(resp_forest)) resp_forest$herbHeight <- NULL
  if("seedBank" %in% names(resp_forest)) resp_forest$seedBank <- NULL

  return(resp_forest)
}
