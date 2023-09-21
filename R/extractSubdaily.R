#' Extracts subdaily output
#' 
#' Given the result of simulations, this function extracts subdaily output corresponding to each simulated day and returns it as a data frame.
#' 
#' @param  x An object returned by simulation functions \code{\link{spwb}}, \code{\link{pwb}} or \code{\link{growth}}.
#' @param output See options in section details.
#' @param dates A date vector indicating the subset of simulated days for which subdaily output is desired.
#' 
#' @details This function only works when simulations have been carried using control option 'subdailyResults = TRUE' (see \code{\link{defaultControl}}). Subdaily simulation results will then be stored as elements of the a list called 'subdaily' in the simulation output. Function \code{extractSubdaily} will assemble subdaily results from this list and return them as a data frame. Options for parameter 'output' are the following:
#' \itemize{
#'   \item{Functions pwb() and spwb(): "E","Ag","An","dEdP","RootPsi","StemPsi","LeafPsi","StemPLC","StemRWC","LeafRWC","StemSympRWC","LeafSympRWC","PWB", "Temperature", "ExtractionInst".}
#'   \item{Additional options for shade and sunlit leaves in pwb() and spbw(): Either "SunlitLeaves$x" or "ShadeLeaves$x" where 'x' is one of the following: "Abs_SWR","Abs_PAR","Net_LWR","E","Ag","An","Ci","Gsw","VPD","Temp","Psi","iWUE".}
#'   \item{Additional options for function growth(): "GrossPhotosynthesis", "MaintenanceRespiration", "GrowthCosts", "LabileCarbonBalance","SugarLeaf", "SugarSapwood", "StarchLeaf", "StarchSapwood","SugarTransport".}
#' }
#' 
#' @return A data frame with a column 'datetime' and as many columns as plant cohorts.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{summary.spwb}}
extractSubdaily<-function(x, output = "E", dates = NULL)  {
  leafTypes= c("Abs_PAR", "Abs_SWR","Net_LWR","E","Ag","An","Ci","Gsw","VPD","Temp","Psi","iWUE")  
  sunlitTypes = paste("SunlitLeaves",leafTypes, sep="$")
  shadeTypes = paste("ShadeLeaves",leafTypes, sep="$")
  plantTypes = c("E","Ag","An","dEdP","RootPsi",
                "StemPsi","LeafPsi","StemPLC","StemRWC","LeafRWC","StemSympRWC","LeafSympRWC",
                "LeafSympPsi", "StemSympPsi","PWB",
                "StemPLC")
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
  return(m)
}
