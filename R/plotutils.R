.multiple_dynamics<-function(x, xlab = "", ylab=NULL, ylim = NULL, labels = NULL) {
  if(is.null(labels)) labels = colnames(x)
  df = data.frame(Y = as.vector(x), Date = as.Date(rownames(x)),
                  Cohort = gl(length(colnames(x)), nrow(x), labels=labels))
  g<-ggplot(df, aes(x=Date, y=Y))+
    geom_line(aes(col=Cohort, linetype = Cohort))+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}
.single_dynamics<-function(x, xlab="", ylab=NULL, ylim = NULL) {
  df = data.frame(Y = x, Date = as.Date(names(x)))
  g<-ggplot(df, aes(x = Date, y= Y))+
    geom_line()+
    scale_color_discrete(name="")+
    scale_linetype_discrete(name="")+
    xlab(xlab)
  if(!is.null(ylim)) g <- g+ylim(ylim)
  if(!is.null(ylab)) g <- g+ylab(ylab)
  return(g)
}