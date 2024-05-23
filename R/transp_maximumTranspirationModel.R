#' Maximum transpiration vs. LAI
#' 
#' Builds a model of maximum transpiration (Tmax) over potential evapotranspiration (PET) for increasing leaf area index (LAI) values for each plant cohort.
#' 
#' @param x An object of class \code{\link{spwbInput}}, built using the 'Sperry' transpiration mode.
#' @param meteo A data frame with daily meteorological data series.
#' @param latitude Latitude (in degrees).
#' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
#' @param LAI_seq Sequence of stand LAI values to be tested.
#' @param draw Logical flag to indicate plotting of results.
#' 
#' @details This function performs a meta-modelling exercise using the Sperry transpiration model, with the aim to estimate coefficients 
#' for the equation used in the Granier transpiration model (Granier et al. 1999). The model to be fitted is: \code{y ~ a*LAI + b*LAI^2}, 
#' where \code{y} is the ratio between maximum transpiration (Tmax) and Penman's potential evapotranspiration (PET) and \code{LAI} is the stand LAI. 
#' Unlike the original equation of Granier et al. (1999), we fit a zero intercept model so that LAI = 0 translates into zero plant transpiration. 
#' 
#' The function fits the model for each cohort separately, assuming it represents the whole stand. 
#' For each stand LAI value in the input sequence, the function uses simulations with Sperry transpiration and the input weather to estimate \code{y = Tmax/PET} 
#' as a function of stand's LAI (deciduous stands include leaf phenology). 
#' Once simulations have been conducted for each stand LAI value, the function fits a Generalized Linear Model with the above equation, 
#' assuming a Gamma distribution of residuals and an identity link.
#' 
#' The coefficients of the model can be used to parametrize Granier's transpiration, 
#' since coefficients \code{a} and \code{b} in the equation above correspond to parameters \code{Tmax_LAI} and \code{Tmax_LAIsq}, 
#' respectively (see \code{\link{SpParamsMED}}).
#' 
#' @return Returns a list with as many elements as plant cohorts, each element being a \code{\link{glm}} model.
#'
#' @references 
#' Granier A, \enc{Bréda}{Breda} N, Biron P, Villette S (1999) A lumped water balance model to evaluate duration and intensity of drought constraints in forest stands. Ecol Modell 116:269–283. https://doi.org/10.1016/S0304-3800(98)00205-1. 
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{spwb}}, \code{\link{transp_transpirationGranier}}, \code{\link{transp_transpirationSperry}}, \code{\link{SpParamsMED}}
#' 
#' @examples 
#' \donttest{
#' #Load example daily meteorological data
#' data(examplemeteo)
#' 
#' # Load example plot plant data
#' data(exampleforest)
#' 
#' # Load default species parameters
#' data(SpParamsMED)
#'
#' # Define soil with default soil params
#' examplesoil <- defaultSoilParams(4)
#' 
#' # Initialize control parameters for 'Sperry' transpiration mode
#' control <- defaultControl(transpirationMode="Sperry")
#' 
#' # Initialize input
#' x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
#' 
#' # Estimate maximum transpiration ratio models for each cohort
#' # Weather is subset to speed-up results
#' m <- transp_maximumTranspirationModel(x2, examplemeteo[1:10,], 
#'                                      41.82592, elevation = 100, 
#'                                      slope = 0, aspect = 0)
#' 
#' # Inspect the model for first cohort
#' m[[1]]
#' }
#' 
transp_maximumTranspirationModel<-function(x, meteo, latitude, elevation, slope, aspect,
                                           LAI_seq = c(0.1,0.25, seq(0.5, 10, by=0.5)),
                                           draw = TRUE) {
  
  
  if("dates" %in% names(meteo)) {
    dates <- meteo$dates
  } else {
    dates <- row.names(meteo)
  }
  #Calculate PET using penman
  if("PET" %in% names(meteo)) meteo$PET = NULL
  PET <- numeric(nrow(meteo))
  
  for (i in 1:length(meteo[['MinTemperature']])) {
    PET[i] <- meteoland::penman(
      latrad = latitude*pi/180,
      elevation = elevation,
      slorad = slope*pi/180,
      asprad = aspect*pi/180,
      J = meteoland::radiation_dateStringToJulianDays(dates[i]),
      Tmin = meteo[['MinTemperature']][i],
      Tmax = meteo[['MaxTemperature']][i],
      RHmin = meteo[['MinRelativeHumidity']][i],
      RHmax = meteo[['MaxRelativeHumidity']][i],
      R_s = meteo[['Radiation']][i],
      u = meteo[['WindSpeed']][i],
      z = 2
    )
  }
  
  #Exclude days without PET
  meteo = meteo[!is.na(PET), ]
  PET = PET[!is.na(PET)]

  #Subsample days from PET
  # PET_cut = cut(PET,breaks = seq(0, max(PET), length.out = 20))
  # nPET = table(PET_cut)
  # PETw = as.numeric(1/nPET[PET_cut])
  # PETw = PETw/sum(PETw)
  # s = sample(1:nrow(meteo),ndays,replace = TRUE, prob = PETw)
  # s = s[order(PET[s])]
  # meteo = meteo[s,]
  # PET =PET[s]
  
  ncoh = nrow(x$above)
  ndays = nrow(meteo)
  nlai = length(LAI_seq)
  
  cohnames <- row.names(x$cohorts)
  LAItotal <- sum(x$above$LAI_live)
  
  xIni = x
  xIni$control$modifyInput = FALSE
  xIni$control$unlimitedSoilWater = TRUE
  xIni$control$stemCavitationRecovery = "total"
  xIni$control$leafCavitationRecovery = "total"
  xIni$control$verbose = FALSE
  Tmax = matrix(NA, nrow=ndays, ncol = nlai)  
  colnames(Tmax) = LAI_seq
  rownames(Tmax) = row.names(meteo)
  LAI = matrix(NA, nrow=ndays, ncol = nlai)  
  colnames(LAI) = LAI_seq
  rownames(LAI) = row.names(meteo)
  
  s_res = vector("list", nlai)
  pb = txtProgressBar(0, nlai, style=3)
  for(j in 1:nlai) {
    customParams = LAI_seq[j]*(x$above$LAI_live/LAItotal)
    names(customParams) = paste0(cohnames,"/LAI_live")
    xlai = modifyInputParams(xIni, customParams, FALSE)
    s_res[[j]] = spwb(xlai, meteo,
                      latitude = latitude, 
                      elevation = elevation, slope = slope, aspect = aspect)
    Tmax[,j] = s_res[[j]]$WaterBalance$Transpiration
    LAI[,j] = s_res[[j]]$Stand$LAI
    setTxtProgressBar(pb, j)
  }
  TmaxRatio = sweep(Tmax,1,PET,"/")
  Tmaxratiovec = as.vector(TmaxRatio)
  laivec = as.vector(LAI)
  df = data.frame(y=Tmaxratiovec, LAI = laivec, Prec = meteo$Precipitation)
  df = df[df$Prec==0,] #Exclude precipitation days
  df = df[!is.na(df$y),, drop=FALSE] # Exclude missing ratio
  df = df[(df$y > 0.0) & (df$y < 1.0),, drop=FALSE] # Exclude extreme ratio
  mod <- glm(y ~ -1 + LAI + I(LAI^2), 
                   start = c(0.134,-0.006),
                   data =df, family=Gamma(link="identity"))
  if(draw==TRUE) {
    TmaxPETGranier = -0.006*(LAI_seq^2)+0.134*LAI_seq
    plot(LAI_seq, TmaxPETGranier, type="l", col="gray", lwd=2, 
         xlab = "Stand's Leaf Area Index", ylab = "Tmax/PET", ylim=c(0,1))
    df2<-data.frame(LAI = LAI_seq)
    lines(df2$LAI, predict(mod, newdata = df2), col="black", lwd=2)
    legend("topright", legend=c("Granier's equation", "Forest stand"),
           lwd=2, lty=1, col=c("gray", "black"),bty="n", cex =0.8)
  }
  return(mod)
}
