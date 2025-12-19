#' Plant regeneration
#'
#' Functions to simulate annual plant regeneration from seed recruitment or from resprouting
#' 
#' @param forest An object of class \code{\link{forest}}.
#' @param seedBank A data frame with columns 'Species' and 'Percent', describing a seed bank.
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsDefinition}}).
#' @param refillSpecies A string vector of species names corresponding to seed rain to refill seed bank.
#' @param refillPercent A numeric vector of indicating the percentage of seed bank refilling (if missing then seed bank is set to 100%).
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
#' control$recruitmentMode = "annual/deterministic" 
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
    # Whether the tree reaches seed production diameter or height
    # We could use the Dmat ~ Dmax relationship in https://onlinelibrary.wiley.com/doi/10.1111/ele.14500
    spd <- species_parameter(treeSpp, SpParams, "SeedProductionDiameter", fillMissing = FALSE)
    selDiam  <- (forest$treeData$DBH > spd)
    sph <- species_parameter(treeSpp, SpParams, "SeedProductionHeight", fillMissing = FALSE)
    sph[is.na(sph)] <- control$seedProductionTreeHeight
    selHeight <- (forest$treeData$Height > sph)
    selDiam[is.na(selDiam)] <- selHeight[is.na(selDiam)] # Replace missing selection values with height criterion
    treeSpp <- unique(treeSpp[selDiam])
  }
  if(control$shrubDynamics) {
    shrubSpp <- forest$shrubData$Species
    if(length(shrubSpp)>0) {
      sph <- species_parameter(shrubSpp, SpParams, "SeedProductionHeight", fillMissing = FALSE)
      sph[is.na(sph)] <- control$seedProductionShrubHeight
      selHeight <- forest$shrubData$Height > sph
      shrubSpp <- unique(shrubSpp[selHeight])
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
regeneration_germination <- function(forest, SpParams, control) {
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
  tree_recr_prob <- species_parameter(treeSpp, SpParams, "ProbRecr")
  tree_recr_prob[is.na(tree_recr_prob)] <- control$probRecr
  shrub_recr_prob <- species_parameter(shrubSpp, SpParams, "ProbRecr")
  shrub_recr_prob[is.na(shrub_recr_prob)] <- control$probRecr
  if(control$recruitmentMode %in% c("annual/stochastic", "daily/stochastic", "stochastic")) {
    treePercent <- rbinom(length(treeSpp), size = 1, prob = tree_recr_prob)*treePercent
    shrubPercent <- rbinom(length(shrubSpp), size = 1, prob = shrub_recr_prob)*shrubPercent
  } else {
    treePercent <- treePercent*tree_recr_prob
    shrubPercent <- shrubPercent*shrub_recr_prob
  }
  nseedling <- length(seedSpp)
  # Create new seedling data
  seedlingBank <- data.frame(Species = c(treeSpp, shrubSpp),
                             Percent = c(treePercent, shrubPercent),
                             Age = as.numeric(rep(0, nseedling)),
                             Z50 = as.numeric(rep(0, nseedling)),
                             Z95 = as.numeric(rep(0, nseedling)),
                             Z100 = as.numeric(rep(NA, nseedling)))
  return(seedlingBank)
}

#' @rdname regeneration
#' @keywords internal
regeneration_seedlings_daily <- function(forest, SpParams, control,
                                         growthResult, verbose = FALSE) {

  
  newSeedlingBank <- regeneration_germination(forest, SpParams, control)
  
  # If exists, merge with existing forest seedling bank 
  if(!is.null(forest$seedlingBank)) {
    seedlingBank <- rbind(newSeedlingBank, forest$seedlingBank)
  } else {
    seedlingBank <- newSeedlingBank
  }
  
  if(nrow(seedlingBank)>0) {
    # Determine thresholds
    minFPAR <- species_parameter(seedlingBank$Species, SpParams, "MinFPARRecr")
    minTemp <- species_parameter(seedlingBank$Species, SpParams, "MinTempRecr")
    minFPAR[is.na(minFPAR)] <- control$minFPARRecr
    minTemp[is.na(minTemp)] <- control$minTempRecr
    
    VCstem_P12 <- species_parameter(seedlingBank$Species, SpParams, "VCstem_P12", TRUE, TRUE)
    VCstem_P50 <- species_parameter(seedlingBank$Species, SpParams, "VCstem_P50", TRUE, TRUE)
    VCstem_P88 <- species_parameter(seedlingBank$Species, SpParams, "VCstem_P88", TRUE, TRUE)
    VCstem_slope = (88.0 - 12.0)/(abs(VCstem_P88) - abs(VCstem_P12))
    Psi_Extract <- species_parameter(seedlingBank$Species, SpParams, "Psi_Extract", TRUE, TRUE)
    Exp_Extract <- species_parameter(seedlingBank$Species, SpParams, "Exp_Extract", TRUE, TRUE)
    Z95adult <- species_parameter(seedlingBank$Species, SpParams, "Z95")
    if("RecrAge" %in% names(SpParams)) RecrAge <- species_parameter(seedlingBank$Species, SpParams, "RecrTreeDensity")
    if("recrAge" %in% names(control)) {
      RecrAge[is.na(RecrAge)] <- control$recrAge
    } else {
      RecrAge[is.na(RecrAge)] <- 5
    }
    
    # Daily increment in root depth
    incrZ95 <- Z95adult/(RecrAge*365.26)

    # Retrieve environmental_cues
    soil_widths <- growthResult$growthInput$soil$widths
    soilPsi <- growthResult$Soil$Psi[, 1:length(soil_widths)]
    LgroundPAR <- growthResult$Stand$LgroundPAR
    months <- as.numeric(format(as.Date(growthResult$weather$dates), "%m"))
    monthlyMinTemp <- tapply(growthResult$weather$MinTemperature, months, FUN="mean", na.rm=TRUE)
    monthlyMaxTemp <- tapply(growthResult$weather$MaxTemperature, months, FUN="mean", na.rm=TRUE)
    monthlyTemp <- 0.606*monthlyMaxTemp + 0.394*monthlyMinTemp
    minMonthTemp <- min(monthlyTemp, na.rm=TRUE)
    
    # Daily root growth and seedling mortality
    for(i in 1:length(LgroundPAR)) {
      #root growth
      seedlingBank$Z95 <- seedlingBank$Z95 + incrZ95
      seedlingBank$Z50 <- exp(log(seedlingBank$Z95)/1.4)
      if(control$truncateRootDistribution) {
        seedlingBank$Z100 <- exp(log(seedlingBank$Z95)/0.95);
      }
      seedling_v <- root_ldrDistribution(seedlingBank$Z50, seedlingBank$Z95, seedlingBank$Z100, 
                                         soil_widths)
      #mortality
      for(j in 1:nrow(seedlingBank)) {
        seedling_psi <- hydraulics_averagePsi(soilPsi[i, ], seedling_v[j, ], exp_extract = Exp_Extract[j], psi_extract = Psi_Extract[j])
        stem_plc <- 1.0 - hydraulics_xylemConductanceSigmoid(seedling_psi, 1.0, VCstem_P50[j], VCstem_slope[j]);
        m_light <- mortality_dailyProbability(LgroundPAR[i]/100, minFPAR[j]/100)
        m_dry <- mortality_dailyProbability(1.0 - stem_plc, 0.4) 
        seedlingBank$Percent[j] <- seedlingBank$Percent[j]*(1 - max(m_light, m_dry)) 
      }
    }
    # Minimum temperature filter 
    recr_selection <- (minMonthTemp > minTemp)
    
    # Remove rows
    recr_selection[seedlingBank$Percent == 0] <- FALSE
    seedlingBank <- seedlingBank[recr_selection, , drop = FALSE]
    row.names(seedlingBank) <- NULL
    
    # Increment seedling age
    seedlingBank$Age <- seedlingBank$Age + 1 
  }
  
  return(seedlingBank)
}


#' @rdname regeneration
#' @keywords internal
regeneration_seedlings <- function(forest, SpParams, control,
                                   minMonthTemp, moistureIndex, verbose = FALSE) {
  
  newSeedlingBank <- regeneration_germination(forest, SpParams, control)
  
  # If exists, merge with existing forest seedling bank 
  if(!is.null(forest$seedlingBank)) {
    seedlingBank <- rbind(newSeedlingBank, forest$seedlingBank)
  } else {
    seedlingBank <- newSeedlingBank
  }
  
  # Calculate light in the ground
  if((nrow(forest$treeData)>0) || (nrow(forest$shrubData)>0)) {
    PARperc <- light_PARground(forest, SpParams)
  } else {
    PARperc <- 100
  } 
  if(verbose) cat(paste0(" [Cold (C): ", round(minMonthTemp,2), " / Moist: ", round(moistureIndex,2), " / FPAR (%): ", round(PARperc,1), "]"))
  
  
  # Determine thresholds
  minFPAR <- species_parameter(seedlingBank$Species, SpParams, "MinFPARRecr")
  minTemp <- species_parameter(seedlingBank$Species, SpParams, "MinTempRecr")
  minMoisture <- species_parameter(seedlingBank$Species, SpParams, "MinMoistureRecr")
  minFPAR[is.na(minFPAR)] <- control$minFPARRecr
  minTemp[is.na(minTemp)] <- control$minTempRecr
  minMoisture[is.na(minMoisture)] <- control$minMoistureRecr
  # Seedling mortality
  recr_selection <- (minMonthTemp > minTemp) & (moistureIndex > minMoisture) & (PARperc > minFPAR)
  recr_selection[seedlingBank$Percent == 0] <- FALSE
  seedlingBank <- seedlingBank[recr_selection, , drop = FALSE]
  row.names(seedlingBank) <- NULL
  # Increment seedling age and update root distribution
  seedlingBank$Age <- seedlingBank$Age + 1 
  Z95adult <- species_parameter(seedlingBank$Species, SpParams, "Z95")
  if("RecrAge" %in% names(SpParams)) RecrAge <- species_parameter(seedlingBank$Species, SpParams, "RecrTreeDensity")
  if("recrAge" %in% names(control)) {
    RecrAge[is.na(RecrAge)] <- control$recrAge
  } else {
    RecrAge[is.na(RecrAge)] <- 5
  }
  seedlingBank$Z95 <- Z95adult*(seedlingBank$Age/RecrAge)
  seedlingBank$Z50 <- exp(log(seedlingBank$Z95)/1.4)
  if(control$truncateRootDistribution) {
    seedlingBank$Z100 <- exp(log(seedlingBank$Z95)/0.95);
  }
  return(seedlingBank)
}

.seedlings2recruits <- function(seedlingBank, SpParams, control) {
  
  # Get recruitment age (with back-compatibility)
  if("RecrAge" %in% names(SpParams)) RecrAge <- species_parameter(seedlingBank$Species, SpParams, "RecrTreeDensity")
  if("recrAge" %in% names(control)) {
    RecrAge[is.na(RecrAge)] <- control$recrAge
  } else {
    RecrAge[is.na(RecrAge)] <- 5
  }
  
  ## Determine if species can recruit 
  recruitment <- seedlingBank$Age >= RecrAge
  recrBank <- seedlingBank[recruitment, , drop = FALSE]
  if(length(recrBank$Species)>0) {
    isTree <- !(species_characterParameter(recrBank$Species, SpParams, "GrowthForm")=="Shrub")
  } else {
    isTree <- logical(0)
  }
  recrBankTree <- recrBank[isTree, , drop = FALSE]  
  recrBankShrub <- recrBank[!isTree, , drop = FALSE]  
  if(!control$shrubDynamics) {
    recrBankShrub <- recrBank[numeric(0), , drop = FALSE]
  }
  
  recr_forest <- emptyforest(ntree = nrow(recrBankTree), nshrub=nrow(recrBankShrub),
                             addcolumns = c("Z100", "Age", "ObsID"))
  
  if(nrow(recrBankTree)>0) {
    recr_forest$treeData$Species <- recrBankTree$Species
    recr_forest$treeData$N <- species_parameter(recrBankTree$Species, SpParams, "RecrTreeDensity")
    recr_forest$treeData$N[is.na(recr_forest$treeData$N)] <- control$recrTreeDensity
    recr_forest$treeData$N <- recr_forest$treeData$N*(recrBankTree$Percent/100.0) #Apply reduction due to seed bank size
    recr_forest$treeData$DBH <- species_parameter(recrBankTree$Species, SpParams, "RecrTreeDBH")
    recr_forest$treeData$DBH[is.na(recr_forest$treeData$DBH)] <- control$recrTreeDBH
    recr_forest$treeData$Height <- species_parameter(recrBankTree$Species, SpParams, "RecrTreeHeight")
    recr_forest$treeData$Height[is.na(recr_forest$treeData$Height)] <- control$recrTreeHeight
    recr_forest$treeData$Z50 <- recrBankTree$Z50
    recr_forest$treeData$Z95 <- recrBankTree$Z95
    recr_forest$treeData$Z100 <- recrBankTree$Z100
    recr_forest$treeData$Age <- recrBankTree$Age 
  }
  if(nrow(recrBankShrub)>0) {
    recr_forest$shrubData$Species <- recrBankTree$Species
    recr_forest$shrubData$Cover <- species_parameter(recrBankShrub$Species, SpParams, "RecrShrubCover")
    recr_forest$shrubData$Cover[is.na(recr_forest$shrubData$Cover)] <- control$recrShrubCover
    recr_forest$shrubData$Cover <- recr_forest$shrubData$Cover*(recrBankShrub$Percent/100.0) #Apply reduction due to seed bank size
    recr_forest$shrubData$Height <- species_parameter(recrBankShrub$Species, SpParams, "RecrShrubHeight")
    recr_forest$shrubData$Height[is.na(recr_forest$shrubData$Height)] <- control$recrShrubHeight
    recr_forest$shrubData$Z50 <- recrBankShrub$Z50
    recr_forest$shrubData$Z95 <- recrBankShrub$Z95
    recr_forest$shrubData$Z100 <- recrBankShrub$Z100
    recr_forest$shrubData$Age <- recrBankShrub$Age 
  }
  
  recr_forest$treeData <- recr_forest$treeData[recr_forest$treeData$N>0, ,drop=FALSE]
  recr_forest$shrubData <- recr_forest$shrubData[recr_forest$shrubData$Cover>0, ,drop=FALSE]
  
  
  if("herbData" %in% names(recr_forest)) recr_forest$herbData <- NULL
  if("herbCover" %in% names(recr_forest)) recr_forest$herbCover <- NULL
  if("herbHeight" %in% names(recr_forest)) recr_forest$herbHeight <- NULL
  if("seedBank" %in% names(recr_forest)) recr_forest$seedBank <- NULL
  if("litterData" %in% names(recr_forest)) recr_forest$litterData <- NULL
  
  #Remaining seedling bank
  recr_forest$seedlingBank <- seedlingBank[!recruitment, ,drop  = FALSE]
  
  return(recr_forest)
}

#' @rdname regeneration
#' @keywords internal
regeneration_recruitment<-function(forest, SpParams, control,
                      minMonthTemp, moistureIndex, verbose = FALSE) {
  
  seedlingBank <- regeneration_seedlings(forest, SpParams, control,
                                         minMonthTemp, moistureIndex, verbose)

  return(.seedlings2recruits(seedlingBank, SpParams, control))
}

#' @rdname regeneration
#' @param growthResult An object of class 'growth'.
#' @keywords internal
regeneration_recruitment_daily<-function(forest, SpParams, control,
                                         growthResult, verbose = FALSE) {
  
  seedlingBank <- regeneration_seedlings_daily(forest, SpParams, control,
                                               growthResult, verbose)
  
  return(.seedlings2recruits(seedlingBank, SpParams, control))
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
  
  if("herbData" %in% names(resp_forest)) resp_forest$herbData <- NULL
  if("herbCover" %in% names(resp_forest)) resp_forest$herbCover <- NULL
  if("herbHeight" %in% names(resp_forest)) resp_forest$herbHeight <- NULL
  if("seedBank" %in% names(resp_forest)) resp_forest$seedBank <- NULL
  if("litterData" %in% names(resp_forest)) resp_forest$litterData <- NULL
  
  return(resp_forest)
}
