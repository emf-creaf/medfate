plot.fordyn<-function(x, type="StandBasalArea", 
                      ylim=NULL, xlab = NULL, ylab=NULL, 
                      ...) {
  TYPES =   c("StandBasalArea", "StandLAI", "StandDensity",
              "SpeciesBasalArea", "SpeciesLAI", "SpeciesDensity",
              "CohortBasalArea", "CohortLAI", "CohortDensity")
  type = match.arg(type,TYPES)  
  i_type = which(TYPES %in% type)

  vars = rep(c("TreeBasalAreaLive", "LeafAreaIndex", "TreeDensityLive"),3)
  tables = c(rep("StandSummary",3),rep("SpeciesSummary",3),rep("CohortSummary",3))
  
  if(is.null(ylab)) ylab = .getYLab(type)
  if(is.null(xlab)) xlab = "Step"
  
  out = x[[tables[i_type]]]
  df = data.frame(Step = out[["Step"]], 
                  y = out[[vars[i_type]]])
  if(tables[i_type]=="SpeciesSummary") df$group = as.character(out[["Species"]])
  else if(tables[i_type]=="CohortSummary") df$group = as.character(out[["Cohort"]])
  df = df[!is.na(df$y),]
  g<-ggplot(df, aes_string(x="Step", y="y"))
  if("group" %in% names(df)) {
    g <- g+ geom_line(aes_string(col="group"))+
      scale_color_discrete(name="")
  } else {
    g <- g+ geom_line()
  }
  if(!is.null(ylim)) g <- g+ylim(ylim)
  g<-g+theme_bw()+ylab(ylab)+xlab(xlab)
  return(g)
}