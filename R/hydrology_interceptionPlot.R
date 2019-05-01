hydrology_interceptionPlot<-function(x, SpParams, ER = 0.05, gdd = NA, throughfall = FALSE){
  
  LAI_coh = plant_LAI(x, SpParams, gdd)
  g_coh = plant_parameter(x, SpParams, "g")
  
  Cm = sum(LAI_coh*g_coh)
  
  p = .parExtinctionProfile(0,x, SpParams, gdd)/100

  precipitation = seq(0.5,50, by=0.5)
  
  m2<-precipitation-hydrology_rainInterception(precipitation, Cm,p,ER=ER[1])
  rt = 100*m2/precipitation
  if(throughfall) {
    plot(precipitation,rt, type="l", axes=TRUE, ylab="Relative throughfall (%)", 
          xlab="Gross rainfall (mm)",  
          lty=1:length(Cm), col="black", ylim=c(0,100))
  } else {
    plot(precipitation,100-rt, type="l", axes=TRUE, ylab="Percentage of intercepted rainfall (%)", 
         xlab="Gross rainfall (mm)",  
         lty=1:length(Cm), col="black", ylim=c(0,100))
  }
  
  if(length(ER)>1) {
    for(i in 2:length(ER)) {
      m2<-precipitation-hydrology_rainInterception(precipitation, Cm,p,ER=ER[i])
      rt = 100*m2/precipitation
      if(throughfall) {
        lines(precipitation, rt, lty=i)
      } else {
        lines(precipitation, 100-rt, lty=i)
      }
    }
    if(throughfall) legend("bottomright",lty=1:length(ER), legend=paste("ER =",ER), bty="n")
    else legend("topright",lty=1:length(ER), legend=paste("ER =",ER), bty="n")
  }
}
