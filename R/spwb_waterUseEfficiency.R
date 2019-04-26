spwb_waterUseEfficiency<-function(x, type = "iWUE", leaves = "average", freq="days") {
  if(!("spwb" %in% class(x)) && !("pwb" %in% class(x))) {
    stop("'x' should be of class 'spwb' or 'pwb'")
  }
}