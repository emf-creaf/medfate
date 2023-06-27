#' Forest dynamics
#' 
#' Function \code{fordyn} implements a forest dynamics model that simulates 
#' growth, mortality, recruitment and (optionally) management actions in a given forest stand 
#' during a period specified in the input climatic data.
#' 
#' @param forest An object of class \code{\link{forest}}. Alternatively, the output of a previous run, if continuing a previous simulation.
#' @param soil An object of class \code{\link{soil}}.
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsDefinition}}).
#' @param meteo A data frame with daily weather data series (see \code{\link{spwb}}).
#' @param control A list with default control parameters (see \code{\link{defaultControl}}).
#' @param latitude Latitude (in degrees).
#' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
#' @param CO2ByYear A named numeric vector with years as names and atmospheric CO2 concentration (in ppm) as values. Used to specify annual changes in CO2 concentration along the simulation (as an alternative to specifying daily values in \code{meteo}).
#' @param management_function A function that implements forest management actions (see details).
#' @param management_args A list of additional arguments to be passed to the \code{management_function}.
#' 
#' @details Function \code{fordyn} simulates forest dynamics for annual time steps, building on other simulation functions. For each simulated year, the function performs the following steps:
#' \enumerate{
#'   \item{Calls function \code{\link{growth}} to simulate daily water/carbon balance, growth and mortality processes and update the forest object.}
#'   \item{If required, calls function \code{management_function}, using as parameters the forest object and \code{management_args}, which may result in a density reduction for existing plant cohorts and/or a set of new planted cohorts.}
#'   \item{Simulate natural recruitment (for species present in the stand or given in a seed rain input).}
#'   \item{Prepares the input of function \code{\link{growth}} for the next annual time step.}
#'   \item{Store forest status, management arguments, and summaries.}
#' }
#' 
#' To enable forest management, the user needs to provide a function that implements it, which is passed to \code{fordyn} via its argument \code{management_function}. Such function should have  the following arguments:
#'   \itemize{
#'     \item{\code{"x"}: the \code{\link{forest}} object representing the stand to be managed.}
#'     \item{\code{"args"}: a list of parameters regulating the behavior of the management function.} 
#'     \item{\code{"verbose"}: a logical flag to enable console output during the execution of the management function.}
#'   }
#' and return a list with the following elements:
#'   \itemize{
#'     \item{\code{"action"}: A string identifying the action performed (e.g. "thinning").}
#'     \item{\code{"N_tree_cut"}: A vector with the density of trees removed.}
#'     \item{\code{"Cover_shrub_cut"}: A vector with the cover of shrubs removed.} 
#'     \item{\code{"planted_forest"}: An object of class \code{\link{forest}} with the new plant cohorts resulting from tree/shrub planting.}
#'     \item{\code{"management_args"}: A list of management arguments to be used in the next call to the management function.}
#'   }
#' 
#' An example of management function is provided in \code{\link{defaultManagementFunction}}.
#' 
#' @return A list of class 'fordyn' with the following elements:
#' \itemize{
#'   \item{\code{"StandSummary"}: A data frame with stand-level summaries (tree basal area, tree density, shrub cover, etc.) at the beginning of the simulation and after each simulated year.}
#'   \item{\code{"SpeciesSummary"}: A data frame with species-level summaries (tree basal area, tree density, shrub cover, etc.) at the beginning of the simulation and after each simulated year.}
#'   \item{\code{"CohortSummary"}: A data frame with cohort-level summaries (tree basal area, tree density, shrub cover, etc.) at the beginning of the simulation and after each simulated year.}
#'   \item{\code{"TreeTable"}: A data frame with tree-cohort data (species, density, diameter, height, etc.) at the beginning of the simulation (if any) and after each simulated year.}
#'   \item{\code{"DeadTreeTable"}: A data frame with dead tree-cohort data (species, density, diameter, height, etc.) at the beginning of the simulation and after each simulated year.}
#'   \item{\code{"CutTreeTable"}: A data frame with cut tree data (species, density, diameter, height, etc.) after each simulated year.}
#'   \item{\code{"ShrubTable"}: A data frame with shrub-cohort data (species, density, cover, height, etc.) at the beginning of the simulation and after each simulated year.}
#'   \item{\code{"DeadShrubTable"}: A data frame with dead shrub-cohort data (species, density, cover, height, etc.) at the beginning of the simulation (if any) and after each simulated year.}
#'   \item{\code{"CutShrubTable"}: A data frame with cut shrub data (species, density, cover, height, etc.) after each simulated year.}
#'   \item{\code{"ForestStructures"}: A list with the \code{\link{forest}} object of the stand at the beginning of the simulation and after each simulated year.}
#'   \item{\code{"GrowthResults"}: A list with the results of calling function \code{\link{growth}} for each simulated year.}
#'   \item{\code{"ManagementArgs"}: A list of management arguments to be used in another call to \code{fordyn}.}
#'   \item{\code{"NextInputObject"}: An object of class \code{growthInput} to be used in a subsequent simulation.}
#'   \item{\code{"NextForestObject"}: An object of class \code{forest} to be used in a subsequent simulation.}
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{growth}}, \code{\link{regeneration}}, \code{\link{plot.growth}}, \code{\link{defaultManagementFunction}}
#' 
#' @references 
#' De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini M, García-Valdés R, Nadal-Sala D, Sabaté S, 
#' Martin-StPaul N, Morin X, D'Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A trait-enabled model to simulate 
#' Mediterranean forest function and dynamics at regional scales. 
#' Geoscientific Model Development 16: 3165-3201 (https://doi.org/10.5194/gmd-16-3165-2023).
#' 
#' @examples 
#' \donttest{
#' #Load example daily meteorological data
#' data(examplemeteo)
#' #Prepare a two-year meteorological data with half precipitation during 
#' #the second year
#' meteo2001 <- examplemeteo
#' meteo2002 <- examplemeteo
#' meteo2002$Precipitation <- meteo2002$Precipitation/2
#' row.names(meteo2002) <- seq(as.Date("2002-01-01"), 
#'                            as.Date("2002-12-31"), by="day")
#' meteo_01_02 <- rbind(meteo2001, meteo2002)
#' 
#' #Load example plot plant data
#' data(exampleforestMED)
#' 
#' #Default species parameterization
#' data(SpParamsMED)
#' 
#' #Initialize control parameters
#' control <- defaultControl("Granier")
#' 
#' #Initialize soil with default soil params (4 layers)
#' examplesoil <- soil(defaultSoilParams(4))
#' 
#' #Call simulation function
#' fd<-fordyn(exampleforestMED, examplesoil, 
#'            SpParamsMED, meteo_01_02, control,
#'            latitude = 41.82592, elevation = 100)
#' 
#' #Stand-level summaries
#' fd$StandSummary
#' 
#' #Tree table by annual steps
#' fd$TreeTable
#' 
#' #Dead tree table by annual steps
#' fd$DeadTreeTable
#' }
#' 
fordyn<-function(forest, soil, SpParams,
                 meteo, control,
                 latitude , elevation = NA, slope = NA, aspect = NA,
                 CO2ByYear = numeric(0),
                 management_function = NULL, management_args = NULL) {
  
  # Modify control parameters
  verboseDyn <- control$verbose
  control$verbose <- FALSE
  control$subdailyResults <- FALSE
  
  if("dates" %in% names(meteo)) {
    dates <- as.Date(meteo$dates)
  } else {
    dates <- as.Date(row.names(meteo))
  }
  years <- as.numeric(format(dates, "%Y"))
  months <- as.numeric(format(dates, "%m"))
  yearsUnique <- unique(years)
  nYears <- length(yearsUnique)
  growthResults <- vector("list", nYears)
  names(growthResults) <- paste0("Step_", 1:nYears)
  forestStructures <- vector("list", nYears+1)
  names(forestStructures) <- c("Initial", paste0("Step_", 1:nYears))
  
  
  #Initialization
  if(inherits(forest, "fordyn")) {
    if(verboseDyn) cat(paste0("Initialisation from previous run\n"))
    xi <- forest$NextInputObject
    forest <- forest$NextForestObject
    #Check number of cohorts of both objects
    if(nrow(xi$above) != (nrow(forest$treeData) + nrow(forest$shrubData))) stop("growthInput and forest objects do not match!") 
  } else {
    if(is.numeric(forest$treeData$Species)) {
      forest$treeData$Species <- .speciesCharacterParameterFromSpIndex(forest$treeData$Species, SpParams, "Name")
    }
    if(is.numeric(forest$shrubData$Species)) {
      forest$shrubData$Species <- .speciesCharacterParameterFromSpIndex(forest$shrubData$Species, SpParams, "Name")
    }
    #Subset columns relevant for fordyn (in case there are other)
    if(control$allowRecruitment) {
      forest$treeData <- forest$treeData[,c("Species","DBH", "Height","N","Z50","Z95")]
      forest$shrubData <- forest$shrubData[,c("Species","Height","Cover", "Z50","Z95")]
    }
    #Fill missing root params
    if(control$fillMissingRootParams) {
      treeSPZ95 <- species_parameter(forest$treeData$Species, SpParams, "Z95")
      treeSPZ50 <- species_parameter(forest$treeData$Species, SpParams, "Z50")
      treeSPZ50[is.na(treeSPZ50)] <- exp(log(treeSPZ95[is.na(treeSPZ50)])/1.4)
      shrubSPZ50 <- species_parameter(forest$shrubData$Species, SpParams, "Z50")
      shrubSPZ95 <- species_parameter(forest$shrubData$Species, SpParams, "Z95")
      shrubSPZ50[is.na(shrubSPZ50)] <- exp(log(shrubSPZ95[is.na(shrubSPZ50)])/1.4)
      forest$treeData$Z50[is.na(forest$treeData$Z50)] <- treeSPZ50[is.na(forest$treeData$Z50)]
      forest$treeData$Z95[is.na(forest$treeData$Z95)] <- treeSPZ95[is.na(forest$treeData$Z95)]
      forest$shrubData$Z50[is.na(forest$shrubData$Z50)] <- shrubSPZ50[is.na(forest$shrubData$Z50)]
      forest$shrubData$Z95[is.na(forest$shrubData$Z95)] <- shrubSPZ95[is.na(forest$shrubData$Z95)]
    }
    xi <- forest2growthInput(forest, soil, SpParams, control)
  }
  forestStructures[[1]] <- forest

  #initial tree/shrub tables
  treeTable <- .createTreeTable(0, NA, xi)
  shrubTable <- .createShrubTable(0, NA, xi)
  deadTreeTable <- .createDeadTreeTable(0, NA, xi)
  deadShrubTable <- .createDeadShrubTable(0, NA, xi)
  cutTreeTable <- treeTable[numeric(),,drop = FALSE]
  cutShrubTable <- shrubTable[numeric(),,drop = FALSE]
  
  #initial summaries
  cohortSummary <-.summarizeCohorts(0, 
                                    treeTable, shrubTable,
                                    deadTreeTable, deadShrubTable,
                                    cutTreeTable, cutShrubTable)
  speciesSummary<-.summarizeSpecies(0,cohortSummary, xi, SpParams)
  standSummary<-.summarizeStand(0,cohortSummary, xi)

  
  #Simulations
  for(iYear in 1:nYears) {
    year <- yearsUnique[iYear]
    if(verboseDyn) cat(paste0("Simulating year ", year, " (", iYear,"/", nYears,"): "))
    meteoYear <- meteo[years==year,]
    monthsYear <- months[years==year]
    # 1.1 Calls growth model
    if(verboseDyn) cat(paste0(" (a) Growth/mortality"))
    Gi <- growth(xi, meteoYear, latitude = latitude, 
                elevation = elevation, slope = slope, aspect = aspect,
                CO2ByYear = CO2ByYear)

    # 1.2 Store growth results
    growthResults[[iYear]] <- Gi
    
    # 1.3 Retrieve modified growth output
    xo <- Gi$growthOutput

    # 2.2 Update dead tree/shrub tables
    deadTreeTableYear <- .createDeadTreeTable(iYear, year, xo)
    deadShrubTableYear <- .createDeadShrubTable(iYear, year, xo)
    
    # 2.2 Update forest structural variables
    isTree <- is.na(xo$above$Cover)
    forest$treeData$N  <- xo$above$N[isTree]
    forest$treeData$DBH  <- xo$above$DBH[isTree]
    forest$treeData$Height  <- xo$above$H[isTree]
    if(control$shrubDynamics) {
      forest$shrubData$Cover  <- xo$above$Cover[!isTree]
      forest$shrubData$Height  <- xo$above$H[!isTree]
    }
    
    # 2.3 Call management function if required
    cutTreeTableYear <- NULL
    cutShrubTableYear <- NULL
    management_result <- NULL
    planted_forest <- emptyforest()
    if(!is.null(management_function)) {
      management_result <- do.call(management_function, list(x = forest, args= management_args, verbose = FALSE))
      if(verboseDyn) cat(paste0(" & management [", management_result$action,"]"))
      # Update forest and xo objects
      forest$treeData$N <- pmax(0,forest$treeData$N - management_result$N_tree_cut)
      xo$above$N[isTree] <- forest$treeData$N
      forest$shrubData$Cover <- pmax(0,forest$shrubData$Cover - management_result$Cover_shrub_cut)
      xo$above$Cover[!isTree] <- forest$shrubData$Cover
      # Update cut tables
      cutTreeTableYear <- .createCutTreeTable(iYear, year, xo, management_result$N_tree_cut)
      cutShrubTableYear <- .createCutShrubTable(iYear, year, xo, management_result$Cover_shrub_cut)
      # Retrieve plantation information
      planted_forest <- management_result$planted_forest
      if(nrow(planted_forest$treeData)>0) {
        for(i in 1:nrow(planted_forest$treeData)) {
          if(length(grep("[a-z]", planted_forest$treeData$Species[i]))==0) planted_forest$treeData$Species[i] <- SpParams$Name[SpParams$SpIndex == planted_forest$treeData$Species[i]]
          planted_forest$treeData$Z50[i] <- species_parameter(planted_forest$treeData$Species[i], SpParams,"RecrZ50")
          planted_forest$treeData$Z95[i] <- species_parameter(planted_forest$treeData$Species[i], SpParams,"RecrZ95")
          if(is.na(planted_forest$treeData$Z50[i])) planted_forest$treeData$Z50[i] <- 250
          if(is.na(planted_forest$treeData$Z95[i])) planted_forest$treeData$Z95[i] <- 500
        }
      }
      if(nrow(planted_forest$shrubData)>0) {
        for(i in 1:nrow(planted_forest$shrubData)) {
          if(length(grep("[a-z]", planted_forest$shrubData$Species[i]))==0) planted_forest$shrubData$Species[i] <- SpParams$Name[SpParams$SpIndex == planted_forest$shrubData$Species[i]]
          planted_forest$shrubData$Z50[i] <- species_parameter(planted_forest$shrubData$Species[i], SpParams,"RecrZ50")
          planted_forest$shrubData$Z95[i] <- species_parameter(planted_forest$shrubData$Species[i], SpParams,"RecrZ95")
          if(is.na(planted_forest$shrubData$Z50[i])) planted_forest$shrubData$Z50[i] <- 100
          if(is.na(planted_forest$shrubData$Z95[i])) planted_forest$shrubData$Z95[i] <- 300
        }
      }
      
      # Store new management arguments (may have changed)
      management_args <- management_result$management_args
    } 
    
    # 3. Simulate species recruitment and resprouting
    if(verboseDyn && (control$allowRecruitment || control$allowResprouting)) cat(paste0(", (b) Recruitment/resprouting"))
    if(control$allowRecruitment) {
      monthlyMinTemp <- tapply(Gi$weather$MinTemperature, monthsYear, FUN="mean", na.rm=TRUE)
      monthlyMaxTemp <- tapply(Gi$weather$MaxTemperature, monthsYear, FUN="mean", na.rm=TRUE)
      monthlyTemp <- 0.606*monthlyMaxTemp + 0.394*monthlyMinTemp
      minMonthTemp <- min(monthlyTemp, na.rm=TRUE)
      moistureIndex <- sum(Gi$WaterBalance$Precipitation, na.rm=TRUE)/sum(Gi$WaterBalance$PET, na.rm=TRUE)
      recr_forest <- recruitment(forest, SpParams, control, minMonthTemp, moistureIndex, verbose = FALSE)
    } else {
      recr_forest <- emptyforest()
    }
    # 3.2. Simulate species resprouting
    if(control$allowResprouting) {
      resp_forest <- resprouting(forest, xo$internalMortality, SpParams, control,
                                 management_result)
    } else {
      resp_forest <- emptyforest()
    }
    
    # 4. Update inputs for next year 
    nyf <- .nextYearForest(forest, xo, SpParams, control, 
                           planted_forest, recr_forest, resp_forest)

    forest <- nyf$forest
    xi <- nyf$xi

    # 6.1 Store current forest state (after recruitment/resprouting)
    forestStructures[[iYear+1]] <- forest
    
    # 6.2 Process summaries (after recruitment/resprouting)
    treeTableYear <- .createTreeTable(iYear, year, xi)
    shrubTableYear <- .createShrubTable(iYear, year, xi)
    cohSumYear <- .summarizeCohorts(iYear,
                                   treeTableYear, shrubTableYear,
                                   deadTreeTableYear, deadShrubTableYear,
                                   cutTreeTableYear, cutShrubTableYear)
    speciesSummary <- rbind(speciesSummary,  .summarizeSpecies(iYear,cohSumYear, xi, SpParams))
    standSummary <- rbind(standSummary,  .summarizeStand(iYear,cohSumYear, xi))
    cohortSummary <- rbind(cohortSummary, cohSumYear)
    
    # 6.3 Update tree/shrub tables (after recruitment)
    treeTable <- rbind(treeTable, treeTableYear)
    shrubTable <- rbind(shrubTable, shrubTableYear)
    deadTreeTable <- rbind(deadTreeTable, deadTreeTableYear)
    deadShrubTable <- rbind(deadShrubTable, deadShrubTableYear)
    if(!is.null(cutTreeTableYear)) cutTreeTable <- rbind(cutTreeTable, cutTreeTableYear)
    if(!is.null(cutShrubTableYear)) cutShrubTable <- rbind(cutShrubTable, cutShrubTableYear)
    
    if(verboseDyn) cat(paste0("\n"))
  }
  res <- list(
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