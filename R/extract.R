.extractSubdaily<-function(x, output = "E", dates = NULL, long_format = FALSE)  {
  leafTypes= c("LAI", "Abs_PAR", "Abs_SWR","Net_LWR","E","Ag","An","Ci","Gsw","VPD","Temp","Psi","iWUE")  
  sunlitTypes = paste("SunlitLeaves",leafTypes, sep="$")
  shadeTypes = paste("ShadeLeaves",leafTypes, sep="$")
  plantTypes = c("E","Ag","An","dEdP","RootPsi",
                "StemPsi","LeafPsi","StemPLC", "LeafPLC","StemRWC","LeafRWC","StemSympRWC","LeafSympRWC",
                "LeafSympPsi", "StemSympPsi","PWB")
  PWBTYPES = c("PlantLAI","Temperature", "CanopyEnergyBalance", "SoilEnergyBalance",
               "ExtractionInst", plantTypes, sunlitTypes, shadeTypes)
  CBTYPES = c("GrossPhotosynthesis", "MaintenanceRespiration", "GrowthCosts", "RootExudation", "LabileCarbonBalance",
              "SugarLeaf", "SugarSapwood", "StarchLeaf", "StarchSapwood","SugarTransport")
  GROWTHTYPES = c(CBTYPES, PWBTYPES)
  if(is.null(dates)) dates = as.Date(names(x$subdaily))
  
  if(("spwb" %in% class(x)) || ("pwb" %in% class(x))) {
    input = x$spwbInput
    output = match.arg(output, PWBTYPES)
  } else {
    input = x$growthInput
    output = match.arg(output, GROWTHTYPES)
  }
  
  numCohorts = nrow(input$above)
  numDates = length(dates)
  numSteps = input$control$ndailysteps
  h = 0 + (0:(numSteps-1))*(24/numSteps)
  minutes = 60*h%%1
  seconds = round(60*minutes%%1)
  minutes = floor(minutes)
  hours = floor(h)
  times = paste(hours,minutes,seconds, sep=":")
  
  if(output %in% plantTypes) {
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = numCohorts+1))
    for(i in 1:numDates) {
      ori = x$subdaily[[as.character(dates[i])]]$PlantsInst[[output]]
      m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori) 
    }
    colnames(m) = c("datetime", row.names(input$above))
  } else if(output=="PlantLAI") {
    ori1 = x$subdaily[[as.character(dates[1])]]$Plants$LAI
    nc = length(ori1)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = nc+1))
    for(i in 1:numDates) {
      ori = x$subdaily[[as.character(dates[i])]]$Plants$LAI
      m[((i-1)*numSteps+1):(i*numSteps), 2:(nc+1)] = matrix(ori, nrow = numSteps, ncol = nc, byrow = TRUE) 
    }
    colnames(m) = c("datetime", row.names(input$above))
  } else if(output %in% c("Temperature", "CanopyEnergyBalance", "SoilEnergyBalance")) {
    ori1 = x$subdaily[[as.character(dates[1])]]$EnergyBalance[[output]]
    ncols = ncol(ori1)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = ncols+1))
    for(i in 1:numDates) {
      ori = x$subdaily[[as.character(dates[i])]]$EnergyBalance[[output]]
      m[((i-1)*numSteps+1):(i*numSteps), 2:(ncols+1)] = ori 
    }
    colnames(m) = c("datetime", colnames(ori1))
  } else if(output=="ExtractionInst") {
    ori1 = x$subdaily[[as.character(dates[1])]]$ExtractionInst
    ncols = nrow(ori1)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = ncols+1))
    for(i in 1:numDates) {
      ori = x$subdaily[[as.character(dates[i])]]$ExtractionInst
      m[((i-1)*numSteps+1):(i*numSteps), 2:(ncols+1)] = t(ori) 
    }
    colnames(m) = c("datetime", rownames(ori1))
  } else if(output %in% sunlitTypes) {
    leafType = strsplit(output,"[$]")[[1]][2]
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = numCohorts+1))
    for(i in 1:numDates) {
      if(leafType=="iWUE") {
          ori1 = x$subdaily[[as.character(dates[i])]]$SunlitLeavesInst$An
          ori2 = x$subdaily[[as.character(dates[i])]]$SunlitLeavesInst$Gsw
          m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori1/ori2) 
      } else{
        ori = x$subdaily[[as.character(dates[i])]]$SunlitLeavesInst[[leafType]]
        m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori) 
      }
    }
    colnames(m) = c("datetime", row.names(input$above))
  } else if(output %in% shadeTypes) {
    leafType = strsplit(output,"[$]")[[1]][2]
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = numCohorts+1))
    for(i in 1:numDates) {
      if(leafType=="iWUE") {
        ori1 = x$subdaily[[as.character(dates[i])]]$ShadeLeavesInst$An
        ori2 = x$subdaily[[as.character(dates[i])]]$ShadeLeavesInst$Gsw
        m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori1/ori2) 
      } else{
        ori = x$subdaily[[as.character(dates[i])]]$ShadeLeavesInst[[leafType]]
        m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori) 
      }
    }
    colnames(m) = c("datetime", row.names(input$above))
  } else if(output %in% CBTYPES) {
    ori1 = x$subdaily[[as.character(dates[1])]]$LabileCarbonBalanceInst[[output]]
    ncols = nrow(ori1)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = ncols+1))
    for(i in 1:numDates) {
      ori = x$subdaily[[as.character(dates[i])]]$LabileCarbonBalanceInst[[output]]
      m[((i-1)*numSteps+1):(i*numSteps), 2:(ncols+1)] = t(ori) 
    }
    colnames(m) = c("datetime", row.names(ori1))
  }
  m$datetime = as.POSIXct(paste(dates[gl(n=numDates, k=numSteps)], times))
  ## To long format
  if(long_format) {
    cohorts <- names(m)[-1]
    df = data.frame("datetime" = rep(as.POSIXct(m$datetime),length(cohorts)),
                    "cohort" = gl(n = length(cohorts), k=nrow(m), labels = cohorts),
                    "Y" = as.vector(as.matrix(m[,-1])))
    names(df)[3] <- output
    m <- df
  }
  return(m)
}

#' Extracts model outputs
#' 
#' Function \code{extract()} extracts daily or subdaily output and returns it as a tidy data frame.
#' 
#' @param x An object returned by simulation functions \code{\link{spwb}}, \code{\link{pwb}} or \code{\link{growth}}.
#' @param level Level of simulation output, either "forest" (stand-level results), "soillayer" (soil layer-level results), "cohort" (cohort-level results), 
#' "sunlitleaf" or "shadeleaf" (leaf-level results)
#' @param output Section of the model output to be explored. See details.
#' @param vars Variables to be extracted (by default, all of them).
#' @param dates A date vector indicating the subset of simulated days for which output is desired.
#' @param subdaily A flag to indicate that subdaily values are desired (see details).
#' 
#' @details 
#' When \code{subdaily = FALSE}, parameter \code{output} is used to restrict the section in \code{x} where variables are located. For
#' example \code{output = "Plants"} will correspond to variables "LAI", "LAIlive", "Transpiration", "StemPLC",... as returned by a call
#' \code{names(x$Plants)}. 
#' 
#' Option \code{subdaily = TRUE} only works when simulations have been carried using control option 'subdailyResults = TRUE' (see \code{\link{defaultControl}}). 
#' When using \code{subdaily = TRUE}, parameter \code{output} is not taken into account, and options for parameter \code{vars} are the following:
#' \itemize{
#'   \item{Variables for \code{level = "forest"} or \code{level = "soillayer"}: Not allowed. An error is raised.}
#'   \item{Variables for \code{level = "cohort"}: "E","Ag","An","dEdP","RootPsi","StemPsi","LeafPsi","StemPLC","StemRWC","LeafRWC","StemSympRWC","LeafSympRWC","PWB".}
#'   \item{Variables for \code{level = "shadeleaf"} and \code{level="sunlitleaf"}: "Abs_SWR","Abs_PAR","Net_LWR","E","Ag","An","Ci","Gsw","VPD","Temp","Psi","iWUE".}
#' }
#' 
#' @return 
#' Function \code{extract()} returns a data frame:
#' \itemize{
#'   \item{If \code{level = "forest"}, columns are "date" and variable names.}
#'   \item{If \code{level = "soillayer"}, columns are "date", "soillayer" and variable names.} 
#'   \item{If \code{level = "cohort"}, \code{level = "sunlitleaf"} or \code{level = "shadeleaf"}, columns are "date", "cohorts", "species" and variable names.} 
#'   \item{If \code{subdaily = TRUE}, columns are "datetime", "cohorts", "species" and variable names.} 
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @name extract
#' @examples
#' #Load example daily meteorological data
#' data(examplemeteo)
#' 
#' #Load example plot plant data
#' data(exampleforest)
#' 
#' #Default species parameterization
#' data(SpParamsMED)
#' 
#' #Define soil with default soil params (4 layers)
#' examplesoil <- defaultSoilParams(4)
#' 
#' #Initialize control parameters
#' control <- defaultControl("Granier")
#' 
#' #Initialize input
#' x <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
#' 
#' #Call simulation function (ten days)
#' S1<-spwb(x, examplemeteo[1:10, ], latitude = 41.82592, elevation = 100)
#' 
#' #Extracts daily forest-level output as a data frame
#' extract(S1, level = "forest")
#'
#' #Extracts daily soil layer-level output as a data frame
#' extract(S1, level = "soillayer")
#' 
#' #Extracts daily cohort-level output as a data frame
#' extract(S1, level = "cohort")
#' 
#' #Select the output tables/variables to be extracted
#' extract(S1, level ="cohort", output="Plants", vars = c("PlantStress", "StemPLC"))
#' 
#' @seealso \code{\link{summary.spwb}}
#' @export
extract<-function(x, level = "forest", output = NULL, vars = NULL, dates = NULL, subdaily = FALSE)  {
  level <- match.arg(level, c("forest", "soillayer","cohort", "sunlitleaf", "shadeleaf"))
  
  if(inherits(x, "spwb") || inherits(x, "pwb")) {
    cohorts <- x$spwbInput$cohorts
  } else if(inherits(x, "growth")) {
    cohorts <- x$growthInput$cohorts
  }
  cohnames <- row.names(cohorts)
  spnames <- cohorts$Name
  
  if(subdaily) {
    if(level %in% c("forest", "soillayer")) stop("Subdaily results are for levels 'cohort', 'sunlitleaf' and 'shadeleaf' only.")
    if(level=="cohort") {
      cohort_types = c("E","Ag","An","dEdP","RootPsi",
                   "StemPsi","LeafPsi","StemPLC", "LeafPLC","StemRWC","LeafRWC","StemSympRWC","LeafSympRWC",
                   "LeafSympPsi", "StemSympPsi","PWB")
      if(is.null(vars)) {
        vars <- cohort_types
      } else {
        for(i in 1:length(vars)) vars[i] <-match.arg(vars[i], cohort_types)
      }
    }
    if(level %in% c("sunlitleaf", "shadeleaf")) {
      leaf_types= c("LAI", "Abs_PAR", "Abs_SWR","Net_LWR","E","Ag","An","Ci","Gsw","VPD","Temp","Psi","iWUE")  
      if(is.null(vars)) {
        vars <- leaf_types
      } else {
        for(i in 1:length(vars)) vars[i] <-match.arg(vars[i], leaf_types)
      }
    }
    df <- NULL
    for(ivar in 1:length(vars)) {
      var_name <- vars[ivar]
      outputSubdaily <- var_name
      if(level=="sunlitleaf") outputSubdaily <- paste0("SunlitLeaves$", outputSubdaily)
      if(level=="shadeleaf") outputSubdaily <- paste0("ShadeLeaves$", outputSubdaily)
      if(is.null(df)) {
        df <- .extractSubdaily(x, output = outputSubdaily, dates = dates, long_format = TRUE)
        df[["species"]] <-as.character(gl(length(spnames), nrow(df)/length(spnames), labels = spnames))
        df <- df[,c(1,2,4,3)]
        names(df)[4] <- var_name
      } else {
        df_var <- .extractSubdaily(x, output = outputSubdaily, dates = dates, long_format = TRUE)
        df[[var_name]] <- df_var[,3]
      }
    }
    return(df)
  }
  if(is.null(dates)) dates <- row.names(x$WaterBalance)
  else {
    if(!all(dates %in% row.names(x$WaterBalance))) stop("Some dates are outside the range in 'x'")
  }

  
  if(level=="forest") {
    stand_level_names <-c("WaterBalance", "Stand", "Snow",
                          "EnergyBalance", "Temperature","CarbonBalance", "BiomassBalance", "FireHazard")
    if(is.null(output)) {
      output <- stand_level_names
    } else {
      for(i in 1:length(output)) output[i] <-match.arg(output[i], stand_level_names)
    }
    out <- data.frame(date = dates)
    for(n in output) {
      if(n %in% names(x)) {
        M <- x[[n]]
        varnames <- names(M)
        M <- M[rownames(M) %in% dates, varnames, drop = FALSE]
        if(!is.null(vars)) M <- M[,colnames(M) %in% vars, drop = FALSE]
        if(ncol(M)>0) {
          row.names(M) <- NULL
          out <- cbind(out, M)
        }
      }
    }
  } else if (level =="soillayer") {
    layers <- c(1:length(x$spwbInput$soil$widths), "Overall")
    output <- "Soil"
    out <- data.frame(date = rep(dates, length(layers)),
                      soillayer = as.character(gl(length(layers), length(dates), labels = layers)))
    for(n in output) {
      if(n %in% names(x)) {
        P = x[[n]]
        if(is.null(vars)) vars <- names(P)
        for(v in vars) {
          M <- P[[v]]
          M <- M[rownames(M) %in% dates, , drop = FALSE]
          out[[v]] <- as.vector(M)
        }
      }
    }
  } else if (level =="cohort") {
    plant_level_names <-c("Plants", "LabileCarbonBalance","PlantBiomassBalance", 
                          "PlantStructure", "GrowthMortality")
    if(is.null(output)) {
      output <- plant_level_names
    } else {
      for(i in 1:length(output)) output[i] <-match.arg(output[i], plant_level_names)
    }
    out <- data.frame(date = rep(dates, length(cohnames)),
                      cohort = as.character(gl(length(cohnames), length(dates), labels = cohnames)),
                      species = as.character(gl(length(cohnames), length(dates), labels = spnames)))
    for(n in output) {
      if(n %in% names(x)) {
        P = x[[n]]
        if(is.null(vars)) vars <- names(P)
        vars <- vars[vars!="RhizoPsi"]
        for(v in vars) {
          M <- P[[v]]
          M <- M[rownames(M) %in% dates, , drop = FALSE]
          out[[v]] <- as.vector(M)
        }
      }
    }
  } else if (level =="sunlitleaf") {
    leaf_level_names <-c("SunlitLeaves")
    if(is.null(output)) {
      output <- leaf_level_names
    } else {
      for(i in 1:length(output)) output[i] <-match.arg(output[i], leaf_level_names)
    }
    out <- data.frame(date = rep(dates, length(cohnames)),
                      cohort = as.character(gl(length(cohnames), length(dates), labels = cohnames)),
                      species = as.character(gl(length(cohnames), length(dates), labels = spnames)))
    for(n in output) {
      if(n %in% names(x)) {
        P = x[[n]]
        if(is.null(vars)) vars <- names(P)
        for(v in vars) {
          M <- P[[v]]
          M <- M[rownames(M) %in% dates, , drop = FALSE]
          out[[v]] <- as.vector(M)
        }
      }
    }
  } else if (level =="shadeleaf") {
    leaf_level_names <-c("ShadeLeaves")
    if(is.null(output)) {
      output <- leaf_level_names
    } else {
      for(i in 1:length(output)) output[i] <-match.arg(output[i], leaf_level_names)
    }
    out <- data.frame(date = rep(dates, length(cohnames)),
                      cohort = as.character(gl(length(cohnames), length(dates), labels = cohnames)),
                      species = as.character(gl(length(cohnames), length(dates), labels = spnames)))
    for(n in output) {
      if(n %in% names(x)) {
        P = x[[n]]
        if(is.null(vars)) vars <- names(P)
        for(v in vars) {
          M <- P[[v]]
          M <- M[rownames(M) %in% dates, , drop = FALSE]
          out[[v]] <- as.vector(M)
        }
      }
    }
  }
  return(out)
}
