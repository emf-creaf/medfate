.extractSubdaily<-function(x, output = "E", dates = NULL, long_format = FALSE)  {
  leafTypes= c("LAI", "Abs_PAR", "Abs_SWR","Net_LWR","E","Ag","An","Ci","Gsw","VPD","Temp","Psi","iWUE")  
  sunlitTypes = paste("SunlitLeaves",leafTypes, sep="$")
  shadeTypes = paste("ShadeLeaves",leafTypes, sep="$")
  plantTypes = c("E","Ag","An","dEdP","RootPsi",
                "StemPsi","LeafPsi","StemPLC", "LeafPLC","StemRWC","LeafRWC","StemSympRWC","LeafSympRWC",
                "LeafSympPsi", "StemSympPsi","PWB")
  PWBTYPES = c("PlantLAI","Temperature", "SoilTemperature","CanopyEnergyBalance", "SoilEnergyBalance",
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
  } else if(output %in% c("Temperature", "SoilTemperature","CanopyEnergyBalance", "SoilEnergyBalance")) {
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

.getSubdailyOutputUnit<-function(level, variable) {
  u <- NULL
  if(level=="cohort") {
    if(variable %in% c("E", "PWB")) u <- "l/m2"
    else if(variable %in% c("Ag", "An")) u <- "g/m2"
    else if(variable %in% c("dEdP")) u <- "mmol/m2/s/MPa"
    else if(variable %in% c("RootPsi", "StemPsi", "LeafPsi", "LeafSympPsi", "StemSympPsi")) u <- "MPa"
    
  } else if(level %in% c("sunlitleaf", "shadeleaf")){
    if(variable %in% c("LAI")) u <- "m2/m2"
    else if(variable %in% c("Abs_PAR", "Abs_SWR","Net_LWR")) u <- "W/m2"
    else if(variable %in% c("Ag", "An")) u <- "micromol/m2/s"
    else if(variable %in% c("E")) u <- "mmol/m2/s"
    else if(variable %in% c("Ci")) u <- "ppm"
    else if(variable %in% c("Gsw")) u <- "mol/m2/s"
    else if(variable %in% c("VPD")) u <- "kPa"
    else if(variable %in% c("Temp")) u <- "celsius"
    else if(variable %in% c("Psi")) u <- "MPa"
    else if(variable %in% c("iWUE")) u <- "micromol/mol"
  }
  return(u)
}
.addSubdailyOutputUnits<-function(x, level) {
  cnames <- colnames(x)
  cnames <- cnames[!(cnames %in% c("datetime", "cohort", "species", "sunlitleaf", "shadeleaf"))]
  for(i in 1:ncol(x)) {
    u <- .getSubdailyOutputUnit(level, cnames[i])
    if(!is.null(u)) {
      x[[cnames[i]]] <- units::set_units(x[[cnames[i]]], u, mode = "standard")
    }
  }
  return(x)
}

.getDailyOutputUnit<-function(level, output, variable) {
  u <- NULL 
  if(level=="forest") {
    if(output=="WaterBalance") {
      u <- "l/m2"
    } else if(output=="Snow") {
      u <- "l/m2"
    } else if(output=="Stand") {
      if(variable %in% c("LAI", "LAIherb", "LAIlive", "LAIexpanded", "LAIdead")) u <- "m2/m2"
      else if(variable %in% c("Cm")) u <- "l/m2"
      else if(variable %in% c("LgroundPAR", "LgroundSWR")) u <- "%"
    } else if(output=="EnergyBalance") {
      if(variable %in% c("SWRcan", "LWRcan", "LEVcan", "LEFsnow", "Hcan", "Ebalcan", "Hcansoil", "SWRsoil", "LWRsoil", "LEVsoil", "Ebalsoil")) u <- "W/m2"
    } else if(output=="Temperature") {
      if(variable %in% c("Tatm_mean", "Tatm_min", "Tatm_max", "Tcan_mean", "Tcan_min", "Tcan_max", "Tsoil_mean", "Tsoil_min", "Tsoil_max")) u <- "celsius"
    } else if(output=="CarbonBalance") {
      u <- "g/m2"
    } else if(output=="BiomassBalance") {
      u <- "g/m2"
    } else if(output=="FireHazard") {
      if(variable %in% c("Loading_overstory", "Loading_understory")) u <- "kg/m2"
      else if(variable %in% c("CFMC_understory", "CFMC_overstory", "DFMC")) u <- "%"
      else if(variable %in% c("ROS_surface", "ROS_crown")) u <- "m/min"
      else if(variable %in% c("I_b_surface", "I_b_crown")) u <- "kW/m"
      else if(variable %in% c("FL_surface", "FL_crown")) u <- "m"
      else if(variable %in% c("t_r_surface","t_r_crown")) u <- "s"
    }
  } else if(level=="soillayer") {
    if(variable %in% c("ML")) u <- "l/m2"
    else if(variable %in% c("PlantExt", "HydraulicInput")) u <- "l/m2"
    else if(variable %in% c("Psi")) u <- "MPa"
  } else if(level=="cohort") {
    if(variable %in% c("LAI", "LAIlive")) u <- "m2/m2"
    else if(variable %in% c("FPAR", "AbsorbedSWR")) u <- "%"
    else if(variable %in% c("NetLWR")) u <- "W/m2"
    else if(variable %in% c("Transpiration", "PlantWaterBalance")) u <- "l/m2"
    else if(variable %in% c("GrossPhotosynthesis", "NetPhotosynthesis")) u <- "l/m2"
    else if(variable %in% c("PlantPsi", "LeafPsiMin", "LeafPsiMax", "StemPsi", "RootPsi")) u <- "MPa"
    else if(variable %in% c("LFMC", "LeafRWC", "StemRWC")) u <- "%"
    else if(variable %in% c("dEdP")) u <- "mmol / s /m2 /MPa"
    else if(variable %in% c("GrowthCosts", "RootExudation", "LabileCarbonBalance")) u <- "g/g"
    else if(variable %in% c("SugarLeaf", "StarchLeaf", "SugarSapwood", "StarchSapwood")) u <- "mol/l"
    else if(variable %in% c("SugarTransport")) u <- "mol/s"
    else if(variable %in% c("StructuralBiomassBalance", "LabileBiomassBalance", "PlantBiomassBalance", "MortalityBiomassLoss", "CohortBiomassBalance")) u <- "g/m2"
    else if(variable %in% c("LeafBiomass", "SapwoodBiomass", "FineRootBiomass")) u <- "g"
    else if(variable %in% c("LeafArea", "FineRootArea")) u <- "m2"
    else if(variable %in% c("SapwoodArea")) u <- "cm2"
    else if(variable %in% c("HuberValue")) u <- "cm2/m2"
    else if(variable %in% c("RootAreaLeafArea")) u <- "m2/m2"
    else if(variable %in% c("DBH", "Height")) u <- "cm"
    else if(variable %in% c("LAgrowth", "FRAgrowth")) u <- "m2/day"
    else if(variable %in% c("SAgrowth")) u <- "cm2/day"
    else if(variable %in% c("StarvationRate", "DessicationRate", "MortalityRate")) u <- "day-1"
  } else if(level=="sunlitleaf" || level=="shadeleaf") {
    if(variable %in% c("LeafPsiMin", "LeafPsiMax")) u <- "MPa" 
    else if(variable %in% c("TempMin", "TempMax")) u <- "celsius" 
    else if(variable %in% c("GSWMin", "GSWMax")) u <- "mol m-2 s-1" 
  }
  return(u)
}
.addDailyOutputUnits<-function(x, level, output) {
  cnames <- colnames(x)
  cnames <- cnames[!(cnames %in% c("date", "soillayer", "cohort", "species", "sunlitleaf", "shadeleaf"))]
  for(i in 1:ncol(x)) {
    u <- .getDailyOutputUnit(level, output, cnames[i])
    if(!is.null(u)) {
      x[[cnames[i]]] <- units::set_units(x[[cnames[i]]], u, mode = "standard")
    }
  }
  return(x)
}
.extractInner<-function(x, level = "forest", output = NULL, vars = NULL, dates = NULL, subdaily = FALSE, addunits = FALSE)  {
  if(inherits(x, "aspwb")) {
    level <- match.arg(level, "soillayer")
    subdaily <- FALSE
  } else {
    level <- match.arg(level, c("forest", "soillayer","cohort", "sunlitleaf", "shadeleaf"))
  }
  if(inherits(x, "aspwb")) {
    input <- x$aspwbInput
    cohorts <- character(0)
    cohnames <- character(0)
    spnames <- character(0)
  } else {
    if(inherits(x, "spwb") || inherits(x, "pwb")) {
      input <- x$spwbInput
      cohorts <- x$spwbInput$cohorts
    } else if(inherits(x, "growth")) {
      input <- x$growthInput
      cohorts <- x$growthInput$cohorts
    }
    cohnames <- row.names(cohorts)
    spnames <- cohorts$Name
  }  
  
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
    if(addunits) {
      df <- .addSubdailyOutputUnits(df, level)
    }
    return(df)
  }
  if(is.null(dates)) dates <- row.names(x$WaterBalance)
  else {
    if(!all(dates %in% row.names(x$WaterBalance))) stop("Some dates are outside the range in 'x'")
  }
  
  
  if(level=="forest") {
    stand_level_names <-c("WaterBalance", "Stand", "Snow",
                          "EnergyBalance", "Temperature", "CarbonBalance", "BiomassBalance", "FireHazard")
    if(is.null(output)) {
      output <- stand_level_names
    } else {
      for(i in 1:length(output)) output[i] <-match.arg(output[i], stand_level_names)
    }
    out <- data.frame(date = dates)
    if(addunits) out$date <- as.Date(out$date)
    for(n in output) {
      if(n %in% names(x)) {
        M <- x[[n]]
        varnames <- colnames(M)
        M <- M[rownames(M) %in% dates, varnames, drop = FALSE]
        if(!is.null(vars)) M <- M[,colnames(M) %in% vars, drop = FALSE]
        if(ncol(M)>0) {
          row.names(M) <- NULL
          M <- as.data.frame(M)
          if(addunits) M <- .addDailyOutputUnits(M, level, n)
          out <- cbind(out, M)
        }
      }
    }
  } else if (level =="soillayer") {
    layers <- c(1:length(input$soil$widths), "Overall")
    output <- "Soil"
    out <- data.frame(date = rep(dates, length(layers)),
                      soillayer = layers[gl(length(layers), length(dates))])
    if(addunits) out$date <- as.Date(out$date)
    for(n in output) {
      if(n %in% names(x)) {
        P = x[[n]]
        vars_n <- vars
        if(is.null(vars_n)) vars_n <- names(P)
        for(v in vars_n) {
          M <- P[[v]]
          M <- M[rownames(M) %in% dates, , drop = FALSE]
          out[[v]] <- as.vector(M)
        }
      }
    }
    if(addunits) out <- .addDailyOutputUnits(out, level, "")
  } else if (level =="cohort") {
    plant_level_names <-c("Plants", "LabileCarbonBalance","PlantBiomassBalance", 
                          "PlantStructure", "GrowthMortality")
    if(is.null(output)) {
      output <- plant_level_names
    } else {
      for(i in 1:length(output)) output[i] <-match.arg(output[i], plant_level_names)
    }
    out <- data.frame(date = rep(dates, length(cohnames)),
                      cohort = cohnames[gl(length(cohnames), length(dates))],
                      species = spnames[gl(length(cohnames), length(dates))])
    if(addunits) out$date <- as.Date(out$date)
    for(n in output) {
      if(n %in% names(x)) {
        P <- x[[n]]
        vars_n <- vars
        if(is.null(vars_n)) vars_n <- names(P)
        vars_n <- vars_n[vars_n!="RhizoPsi"]
        for(v in vars_n) {
          M <- P[[v]]
          M <- M[rownames(M) %in% dates, , drop = FALSE]
          out[[v]] <- as.vector(M)
        }
      }
    }
    if(addunits) out <- .addDailyOutputUnits(out, level, "")
  } else if (level =="sunlitleaf") {
    leaf_level_names <-c("SunlitLeaves")
    if(is.null(output)) {
      output <- leaf_level_names
    } else {
      for(i in 1:length(output)) output[i] <-match.arg(output[i], leaf_level_names)
    }
    out <- data.frame(date = rep(dates, length(cohnames)),
                      cohort = cohnames[gl(length(cohnames), length(dates))],
                      species = spnames[gl(length(cohnames), length(dates))],
                      leaftype = "sunlit")
    if(addunits) out$date <- as.Date(out$date)
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
    if(addunits) out <- .addDailyOutputUnits(out, level, "")
  } else if (level =="shadeleaf") {
    leaf_level_names <-c("ShadeLeaves")
    if(is.null(output)) {
      output <- leaf_level_names
    } else {
      for(i in 1:length(output)) output[i] <-match.arg(output[i], leaf_level_names)
    }
    out <- data.frame(date = rep(dates, length(cohnames)),
                      cohort = cohnames[gl(length(cohnames), length(dates))],
                      species = spnames[gl(length(cohnames), length(dates))],
                      leaftype = "shade")
    if(addunits) out$date <- as.Date(out$date)
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
    if(addunits) out <- .addDailyOutputUnits(out, level, "")
  }
  return(out)
}

#' Extracts model outputs
#' 
#' Function \code{extract()} extracts daily or subdaily output and returns it as a tidy data frame.
#' 
#' @param x An object returned by simulation functions \code{\link{spwb}}, \code{\link{aspwb}}, \code{\link{pwb}}, \code{\link{growth}} or \code{\link{fordyn}}.
#' @param level Level of simulation output, either "forest" (stand-level results), "soillayer" (soil layer-level results), "cohort" (cohort-level results), 
#' "sunlitleaf" or "shadeleaf" (leaf-level results)
#' @param output Section of the model output to be explored. See details.
#' @param vars Variables to be extracted (by default, all of them).
#' @param dates A date vector indicating the subset of simulated days for which output is desired.
#' @param subdaily A flag to indicate that subdaily values are desired (see details).
#' @param addunits A flag to indicate that variable units should be added whenever possible.
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
#' extract(S1, level = "forest", addunits = TRUE)
#'
#' #Extracts daily soil layer-level output as a data frame
#' extract(S1, level = "soillayer", addunits = TRUE)
#' 
#' #Extracts daily cohort-level output as a data frame
#' extract(S1, level = "cohort", addunits = TRUE)
#' 
#' #Select the output tables/variables to be extracted
#' extract(S1, level ="cohort", output="Plants", vars = c("PlantStress", "StemPLC"))
#' 
#' @seealso \code{\link{summary.spwb}}
#' @export
extract<-function(x, level = "forest", output = NULL, vars = NULL, dates = NULL, subdaily = FALSE, addunits = FALSE)  {
  if(inherits(x, "aspwb") || inherits(x, "pwb") || inherits(x, "spwb")|| inherits(x, "growth"))  {
    return(.extractInner(x, level, output, vars, dates, subdaily, addunits))
  } else if (inherits(x, "fordyn")) {
    n <- length(x$GrowthResults)
    out <- data.frame()
    for(i in 1:n) {
      out <- rbind(out, .extractInner(x$GrowthResults[[i]], level, output, vars, dates, subdaily, addunits))
    }
    return(out)
  }
  stop("Wrong class for 'x'. Should be one of 'aswpb', 'pwb', 'spwb', 'growth' or 'fordyn'")
}
