.multiple_x<-function(x, y, xlab = "", ylab=NULL, xlim = NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x)
  df = data.frame("X" = as.vector(x), 
                  "Y" = y,
                  "Cohort" = gl(length(colnames(x)), nrow(x), labels=labels))
  g<-ggplot(df, aes_string(x="X", y="Y"))+
    geom_path(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(xlim)) g <- g+xlim(xlim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}
.multiple_y<-function(x, y, xlab = "", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(y)
  df = data.frame("Y" = as.vector(y), "X" = x,
                  "Cohort" = gl(length(colnames(y)), nrow(y), labels=labels))
  g<-ggplot(df, aes_string(x="X", y="Y"))+
    geom_path(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.multiple_dynamics<-function(x, xlab = "", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x)
  df = data.frame("Y" = as.vector(x), 
                  "Date" = as.Date(rownames(x)),
                  "Cohort" = gl(length(colnames(x)), nrow(x), labels=labels))
  g<-ggplot(df, aes_string(x="Date", y="Y"))+
    geom_line(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.single_dynamics<-function(x, xlab="", ylab=NULL, ylim = NULL) {
  df = data.frame("Y" = x, 
                  "Date" = as.Date(names(x)))
  g<-ggplot(df, aes_string(x = "Date", y= "Y"))+
    geom_line()+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.single_subday_dynamics <-function(x, xlab="Time step", ylab=NULL, ylim = NULL) {
  df = data.frame("Y" = x, 
                  "TimeStep" = 1:length(x))
  g<-ggplot(df, aes_string(x = "TimeStep", y= "Y"))+
    geom_line()+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.multiple_subday_dynamics<-function(x, xlab = "Time step", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x)
  df = data.frame("Y" = as.vector(x), 
                  "TimeStep" = as.numeric(rownames(x)),
                  "Cohort" = gl(length(colnames(x)), nrow(x), labels=labels))
  g<-ggplot(df, aes_string(x="TimeStep", 
                           y="Y"))+
    geom_line(aes_string(col="Cohort", linetype = "Cohort"))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}

.multiple_subday_dynamics_sunlit_shade<-function(x_sl, x_sh, xlab = "Time step", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x_sl)
  df_sl = data.frame("Y" = as.vector(x_sl), 
                     "TimeStep" = as.numeric(rownames(x_sl)),
                     "Cohort" = gl(length(colnames(x_sl)), nrow(x_sl), labels=labels),
                     "LeafType" = "Sunlit", stringsAsFactors = F)
  df_sh = data.frame("Y" = as.vector(x_sh), 
                     "TimeStep" = as.numeric(rownames(x_sh)),
                     "Cohort" = gl(length(colnames(x_sh)), nrow(x_sh), labels=labels),
                     "LeafType" = "Shade", stringsAsFactors = F)
  df = as.data.frame(rbind(df_sl, df_sh), stringsAsFactors = F)
  df$LeafType = factor(df$LeafType, levels =c("Sunlit","Shade"))
  g<-ggplot(df, aes_string(x="TimeStep", y="Y"))+
    geom_line(aes_string(col="Cohort", linetype = "Cohort"))+
    facet_wrap(~LeafType, ncol=2)+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    theme_bw()+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}