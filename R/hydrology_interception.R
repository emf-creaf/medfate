#' Rainfall interception
#' 
#' Function \code{hydrology_rainInterception} calculates the amount of rainfall intercepted daily by the canopy, given a rainfall and canopy characteristics. 
#' Two canopy interception models are currently available: the sparse Gash (1995) model and the Liu (2001) model. 
#' In both cases the current implementation assumes no trunk interception.
#' 
#' @param Rainfall A numeric vector of (daily) rainfall.
#' @param Cm Canopy water storage capacity.
#' @param p Proportion of throughfall (normally 1 - c, where c is the canopy cover).
#' @param ER The ratio of evaporation rate to rainfall rate.
#' @param model Rainfall interception model (either \code{"Gash1995"} or \code{"Liu2001"}).
#' 
#' @details 
#' Function \code{hydrology_rainInterception} can accept either vectors or scalars as parameters \code{Cm}, \code{p} and \code{ER}. If they are supplied as vectors they should be of the same length as \code{Rainfall}.
#' 
#' Function \code{hydrology_rainfallIntensity} estimates the rainfall intensity (mm/h) for input values of rainfall and seasonal variation in rainfall intensity (mm/h).
#' 
#' @return 
#' Function \code{hydrology_rainInterception} returns a vector of the same length as \code{Rainfall} containing intercepted rain values. 
#' 
#' Function \code{hydrology_rainfallIntensity} returns a scalar with the rainfall intensity.
#' 
#' @references 
#' Liu (2001). Evaluation of the Liu model for predicting rainfall interception in forests world-wide. - Hydrol. Process. 15: 2341-2360.
#' 
#' Gash (1979). An analytical model of rainfall interception by forests. - Quarterly Journal of the Royal Meteorological Society.                                       
#' 
#' Gash et al. (1995). Estimating sparse forest rainfall interception with an analytical model. - Journal of Hydrology.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{spwb}}
#' 
#' @examples 
#' #Load example plot plant data
#' data(exampleforest)
#' 
#' #Default species parameterization
#' data(SpParamsMED)
#' 
#' #Draw rainfall interception for two values of the E/R ratio
#' hydrology_interceptionPlot(exampleforest, SpParamsMED, ER = c(0.05, 0.2))
#' 
#' @name hydrology_interception
hydrology_rainInterception<-function(Rainfall, Cm, p, ER=0.05, model="Gash1995"){
  model <- match.arg(model, c("Liu2001","Gash1995"))
  if(length(ER)==1) ER =rep(ER, length(Rainfall))
  if(length(Cm)==1) Cm =rep(Cm, length(Rainfall))
  if(length(p)==1) p =rep(p, length(Rainfall))  
  
  if(model=="Gash1995") {
    PG = (-Cm/(ER*(1-p)))*log(1-ER) #Rainfall need to saturate the canopy
    PG[Cm==0 | p==1]=0 #Avoid NAs
    sel = Rainfall > PG #Days where the canopy becomes saturated
    I = rep(NA,length(Rainfall))
    I[sel] = (1-p[sel])*PG[sel] + (1-p[sel])*ER[sel]*(Rainfall[sel]-PG[sel]) 
    I[!sel] = (1-p[!sel])*Rainfall[!sel]  
  } else if(model=="Liu2001") {
    I = Cm*(1-exp(-1*(Rainfall)*((1-p)/Cm)))*(1-(ER/(1-p)))+(ER*Rainfall)
  }
  return(I)
}

#' @param x An object of class \code{\link{spwbInput}}.
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsMED}}).
#' @param gdd Growth degree days (in Celsius).
#' @param throughfall Boolean flag to plot relative throughfall instead of percentage of intercepted rainfall.
#' 
#' @rdname hydrology_interception
hydrology_interceptionPlot<-function(x, SpParams, ER = 0.05, gdd = NA, throughfall = FALSE, model = "Gash1995"){
  
  LAI_coh = plant_LAI(x, SpParams, gdd)
  g_coh = plant_parameter(x, SpParams, "g", TRUE)
  
  Cm = sum(LAI_coh*g_coh)
  
  p = .parExtinctionProfile(0,x, SpParams, gdd)/100

  precipitation = seq(0.5,50, by=0.5)
  
  m2<-precipitation-hydrology_rainInterception(precipitation, Cm,p,ER=ER[1], model = model)
  rt = 100*m2/precipitation
  er = rep(ER[1], length(rt))
  if(length(ER)>1) {
    for(i in 2:length(ER)) {
      m2<-precipitation-hydrology_rainInterception(precipitation, Cm,p,ER=ER[i], model = model)
      rt2 = 100*m2/precipitation
      rt = c(rt, rt2)
      er = c(er, rep(ER[i], length(rt2)))
    }
  }
  if(!throughfall) {
    rt = 100 - rt
    ylab="Percentage of intercepted rainfall (%)"
  } else {
    ylab="Relative throughfall (%)"
  }
  xlab="Gross rainfall (mm)"
  ylim=c(0,100)
  df = data.frame(P = precipitation, RT = rt, ER = paste0("ER = ",er))
  g<-ggplot(df, aes(x=.data$P, y=.data$RT))+
    xlab(xlab)+ylab(ylab)+ylim(ylim)+
    theme_bw()
  if(length(ER)==1) {
    g<-g + geom_path()
  } else {
    g<-g + geom_path(aes(col=.data$ER, linetype=.data$ER))+
      scale_color_discrete(name="")+
      scale_linetype_discrete(name="")
  }
  return(g)
}
