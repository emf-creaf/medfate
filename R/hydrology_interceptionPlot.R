hydrology_interceptionPlot<-function(x, SpParams, ER = 0.05, gdd = NA, throughfall = FALSE){
  
  LAI_coh = plant_LAI(x, SpParams, gdd)
  g_coh = plant_parameter(x, SpParams, "g")
  
  Cm = sum(LAI_coh*g_coh)
  
  p = .parExtinctionProfile(0,x, SpParams, gdd)/100

  precipitation = seq(0.5,50, by=0.5)
  
  m2<-precipitation-hydrology_rainInterception(precipitation, Cm,p,ER=ER[1])
  rt = 100*m2/precipitation
  er = rep(ER[1], length(rt))
  if(length(ER)>1) {
    for(i in 2:length(ER)) {
      m2<-precipitation-hydrology_rainInterception(precipitation, Cm,p,ER=ER[i])
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
  g<-ggplot(df, aes_string(x="P", y="RT"))+
    xlab(xlab)+ylab(ylab)+ylim(ylim)+
    theme_bw()
  if(length(ER)==1) {
    g<-g + geom_path()
  } else {
    g<-g + geom_path(aes_string(col="ER", linetype="ER"))+
      scale_color_discrete(name="")+
      scale_linetype_discrete(name="")
  }
  return(g)
}
