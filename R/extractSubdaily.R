extractSubdaily<-function(x, output = "E", dates = NULL)  {
  leafTypes= c("Abs_SWR","Net_LWR","E","Ag","An","Ci","Gsw","VPD","Temp","Psi","iWUE")  
  sunlitTypes = paste("SunlitLeaves",leafTypes, sep="$")
  shadeTypes = paste("ShadeLeaves",leafTypes, sep="$")
  plantTypes = c("E","Ag","An","dEdP","RootPsi",
                "StemPsi","LeafPsi","StemPLC","StemRWC","LeafRWC","StemSympRWC","LeafSympRWC","PWB")
  PWBTYPES = c("Temperature", "ExtractionInst", plantTypes, sunlitTypes, shadeTypes)
  CBTYPES = c("GrossPhotosynthesis", "MaintenanceRespiration", "GrowthCosts", "CarbonBalance",
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
  
  numDates = length(dates)
  numSteps = input$control$ndailysteps
  h = 0 + (0:(numSteps-1))*(24/numSteps)
  minutes = 60*h%%1
  seconds = round(60*minutes%%1)
  minutes = floor(minutes)
  hours = floor(h)
  times = paste(hours,minutes,seconds, sep=":")
  
  if(output %in% plantTypes) {
    numCohorts = nrow(x$spwbInput$above)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = numCohorts+1))
    for(i in 1:numDates) {
      ori = x$subdaily[[as.character(dates[i])]]$PlantsInst[[output]]
      m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori) 
    }
    colnames(m) = c("datetime", row.names(input$above))
  } else if(output=="Temperature") {
    ori1 = x$subdaily[[as.character(dates[1])]]$EnergyBalance$Temperature
    ncols = ncol(ori1)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = ncols+1))
    for(i in 1:numDates) {
      ori = x$subdaily[[as.character(dates[i])]]$EnergyBalance$Temperature
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
    numCohorts = nrow(input$above)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = numCohorts+1))
    for(i in 1:numDates) {
      if(leafType=="E") {
        ori1 = x$subdaily[[as.character(dates[i])]]$SunlitLeavesInst$Gsw
        ori2 = x$subdaily[[as.character(dates[i])]]$SunlitLeavesInst$VPD
        m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori1*ori2) 
      } else if(leafType=="iWUE") {
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
    numCohorts = nrow(input$above)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = numCohorts+1))
    for(i in 1:numDates) {
      if(leafType=="E") {
        ori1 = x$subdaily[[as.character(dates[i])]]$ShadeLeavesInst$Gsw
        ori2 = x$subdaily[[as.character(dates[i])]]$ShadeLeavesInst$VPD
        m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori1*ori2) 
      } else if(leafType=="iWUE") {
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
    ori1 = x$subdaily[[as.character(dates[1])]]$PlantCBInst[[output]]
    ncols = nrow(ori1)
    m<-data.frame(matrix(nrow = numDates*numSteps, ncol = ncols+1))
    for(i in 1:numDates) {
      ori = x$subdaily[[as.character(dates[i])]]$PlantCBInst[[output]]
      m[((i-1)*numSteps+1):(i*numSteps), 2:(ncols+1)] = t(ori) 
    }
    colnames(m) = c("datetime", row.names(ori1))
  }
  m$datetime = as.character(as.POSIXct(paste(dates[gl(n=numDates, k=numSteps)], times)))
  return(m)
}