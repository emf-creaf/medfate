spwbgrid<-function(y, SpParams, meteo, dates = NULL,
                  summaryFreq = "years", trackSpecies = numeric(), 
                  control = defaultControl()) {

  #check input
  if(!inherits(y, "SpatialGridLandscape"))
    stop("'y' has to be of class 'SpatialGridLandscape'.")
  if(!inherits(meteo,"SpatialGridMeteorology") && 
     !inherits(meteo,"data.frame")) 
    stop("'meteo' has to be of class 'SpatialGridMeteorology' or 'data.frame'.")
  if(!is.null(dates)) if(!inherits(dates, "Date")) stop("'dates' has to be of class 'Date'.")
  
  elevation = y@data$elevation
  
  if(inherits(meteo,"data.frame")) {
    oneMeteoCell = ("MeanTemperature" %in% names(meteo))
    datesMeteo = as.Date(row.names(meteo))
  } else {
    datesMeteo = meteo@dates
  }
  if(is.null(dates)) {
    dates = datesMeteo
  } else {
    if(sum(dates %in% datesMeteo)<length(dates)) 
      stop("Dates in 'dates' is nnot a subset of dates in 'meteo'.")
  }
  date.factor = cut(dates, breaks=summaryFreq)
  df.int = as.numeric(date.factor)
  nDays = length(dates)
  nCells = length(y@forestlist) 
  nSummary = sum(table(date.factor)>0)
  
  #Print information area
  cat("\n------------  spwbgrid ------------\n")
  cat(paste("Grid cells: ", nCells,", area: ", areaSpatialGrid(y)/10000," ha\n", sep=""))
  cat(paste("Number of days to simulate: ",nDays,"\n", sep=""))
  cat(paste("Number of summaries: ", nSummary,"\n\n", sep=""))
  
    
  if(control$verbose) cat(paste("Preparing spwb input"))
  spwbInputList = vector("list", nCells)
  for(i in 1:nCells) {
    yid = y@forestlist[[i]]
    soil = y@soillist[[i]]
    if((!is.na(yid)) && (!is.na(soil))) {             
      xi = forest2spwbInput(yid, soil, SpParams, control)
      spwbInputList[[i]] = xi
    }  else {
      spwbInputList[[i]] = NA
    }
  }
  cat(paste(" - number of cells with spwbInput == NA: ", sum(is.na(spwbInputList)),"\n\n", sep=""))
  

  #Output matrices
  Runon = matrix(0,nrow=nCells, ncol=nSummary)
  colnames(Runon) = levels(date.factor)[1:nSummary]
  Runoff = Runon
  Infiltration = Runon
  NetRain =Runon
  DeepDrainage = Runon
  Esoil = Runon
  Eplant = Runon
  W1mean = Runon
  W2mean = Runon
  W3mean = Runon
  LandscapeBalance = data.frame(Rainfall = rep(0, nSummary),
                                Interception = rep(0, nSummary),
                                Esoil = rep(0,nSummary),
                                Eplant = rep(0, nSummary),
                                Runoff = rep(0, nSummary),
                                DeepDrainage= rep(0, nSummary))
                                
  
  nTrackSpecies = length(trackSpecies)
  if(nTrackSpecies>0) {
    DI = array(0.0,dim=c(nCells, nSummary, nTrackSpecies))
    Transpiration = array(0.0,dim=c(nCells, nSummary, nTrackSpecies))
  } else {
    DI = NULL
    Transpiration = NULL
  }
  
  ifactor = 1
  gridGDD = rep(0,nCells)
  for(day in 1:nDays) {
    cat(paste("Day #", day))
    i = which(datesMeteo == dates[day]) #date index in meteo data
    doy = as.numeric(format(dates[day],"%j"))
    if(inherits(meteo,"SpatialGridMeteorology")) {
      ml = meteo@data[[i]]
      gridMeanTemperature = as.numeric(ml["MeanTemperature"])
      gridPET = as.numeric(ml["PET"])
      gridPrecipitation = as.numeric(ml["Precipitation"])
      gridRadiation = as.numeric(ml["Radiation"])
    } else {
      if(!oneMeteoCell) { # Read meteo grid from file
        f = paste(meteo$dir[i], meteo$filename[i],sep="/")
        if(!file.exists(f)) stop(paste("Meteorology file '", f,"' does not exist!", sep=""))
        ml = readmeteorologygrid(f)
        gridMeanTemperature = as.numeric(ml["MeanTemperature"])
        gridPET = as.numeric(ml["PET"])
        gridPrecipitation = as.numeric(ml["Precipitation"])
        gridRadiation = as.numeric(ml["Radiation"])
      } else { # repeat values for all cells
        gridMeanTemperature = rep(meteo[i,"MeanTemperature"], nCells)
        gridPET = rep(meteo[i, "PET"], nCells)
        gridPrecipitation = rep(meteo[i,"Precipitation"], nCells)
        gridRadiation = rep(meteo[i, "Radiation"], nCells)
      }
    }
    gridER = rep(.er(doy),nCells) #ER
    gridGDD = gridGDD + pmax(gridMeanTemperature - 5.0, 0.0) #Increase GDD
    if(doy>=365) gridGDD = rep(0, nCells) #Reset GDD 
    df = .spwbgridDay(y@lct, spwbInputList, y@soillist, 
                     y@waterOrder, y@queenNeigh, y@waterQ,
                     gridGDD, gridPET, gridPrecipitation, gridER, gridRadiation, elevation,
                     trackSpecies)      
    if(df.int[day]==ifactor) {
      Runon[,ifactor] = Runon[,ifactor] + df$WaterBalance$Runon
      Runoff[,ifactor] = Runoff[,ifactor] + df$WaterBalance$Runoff
      NetRain[,ifactor] = NetRain[,ifactor] + df$WaterBalance$NetRain
      Infiltration[,ifactor] = Infiltration[,ifactor] + df$WaterBalance$Infiltration
      DeepDrainage[,ifactor] = DeepDrainage[,ifactor] + df$WaterBalance$DeepDrainage
      Esoil[,ifactor] = Esoil[,ifactor] + df$WaterBalance$Esoil
      Eplant[,ifactor] = Eplant[,ifactor] + df$WaterBalance$Eplant
      W1mean[,ifactor] = W1mean[,ifactor] + df$WaterBalance$W1
      W2mean[,ifactor] = W2mean[,ifactor] + df$WaterBalance$W2
      W3mean[,ifactor] = W3mean[,ifactor] + df$WaterBalance$W3
      Transpiration[,ifactor,] = Transpiration[,ifactor,]+df$Transpiration
      DI[,ifactor,] = DI[,ifactor,]+pmax((0.5-df$DDS)/0.5,0.0)
      #Landscape balance
      LandscapeBalance$Rainfall[ifactor]= LandscapeBalance$Rainfall[ifactor] + sum(gridPrecipitation)
      LandscapeBalance$Interception[ifactor]= LandscapeBalance$Interception[ifactor] + (sum(gridPrecipitation)-sum(df$WaterBalance$NetRain, na.rm=T))
      LandscapeBalance$Runoff[ifactor] = LandscapeBalance$Runoff[ifactor] + df$RunoffExport
      LandscapeBalance$DeepDrainage[ifactor] = LandscapeBalance$DeepDrainage[ifactor] + sum(df$WaterBalance$DeepDrainage, na.rm=T)
      LandscapeBalance$Esoil[ifactor] = LandscapeBalance$Esoil[ifactor] + sum(df$WaterBalance$Esoil, na.rm=T)
      LandscapeBalance$Eplant[ifactor] = LandscapeBalance$Eplant[ifactor] + sum(df$WaterBalance$Eplant, na.rm=T)
    } else {
      W1mean[,ifactor] = W1mean[,ifactor]/nDays
      W2mean[,ifactor] = W2mean[,ifactor]/nDays
      W3mean[,ifactor] = W3mean[,ifactor]/nDays
      DI[,ifactor,] = DI[,ifactor,]/nDays
      ifactor = ifactor+1
    }
    if(control$verbose) cat("\n")
  }
  #Average last summary
  W1mean[,ifactor] = W1mean[,ifactor]/nDays
  W2mean[,ifactor] = W2mean[,ifactor]/nDays
  W3mean[,ifactor] = W3mean[,ifactor]/nDays
  DI[,ifactor,] = DI[,ifactor,]/nDays
  cat("\n------------  spwbgrid ------------\n")
    
  l <- list(grid = y@grid, LandscapeBalance = LandscapeBalance, NetRain = NetRain, Runon = Runon, Runoff=Runoff, Infiltration=Infiltration, 
                    DeepDrainage = DeepDrainage, Esoil = Esoil, Eplant = Eplant, W1mean = W1mean,
                    W2mean = W2mean, W3mean = W3mean,  DI = DI, Transpiration = Transpiration)
  class(l)<-c("spwbgrid","list")
  return(l)
}