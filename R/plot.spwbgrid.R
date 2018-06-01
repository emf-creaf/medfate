plot.spwbgrid<-function(x, type = "Runon", summaryIndex = 1, spIndex = NULL, ...) {
  
  grid = x$grid
  y = x[[type]]
  if(type %in% c("DI","Transpiration")) {
    y = y[, summaryIndex, spIndex]
  } else {
    y = y[,summaryIndex]
  }
  spplot(SpatialGridDataFrame(grid, data.frame(y)),...)
}