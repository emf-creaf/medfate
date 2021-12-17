stand_dominantTreeHeight<-function(x) {
  h = x$treeData$Height
  n = x$treeData$N
  o <-order(h, decreasing=TRUE)
  h = h[o]
  n = n[o]
  ncum = 0
  for(i in 1:length(h)) {
    ncum = ncum + n[i]
    if(ncum>100) return(sum(h[1:i]*n[1:i])/sum(n[1:i]))
  }
  return(sum(h*n)/sum(n))
}
stand_dominantTreeDiameter<-function(x) {
  dbh = x$treeData$DBH
  n = x$treeData$N
  o <-order(dbh, decreasing=TRUE)
  dbh = dbh[o]
  n = n[o]
  ncum = 0
  for(i in 1:length(dbh)) {
    ncum = ncum + n[i]
    if(ncum>100) return(sum(dbh[1:i]*n[1:i])/sum(n[1:i]))
  }
  return(sum(dbh*n)/sum(n))
}
stand_hartBeckingIndex<-function(x) {
  n = x$treeData$N
  return((10000/stand_dominantTreeHeight(x))*sqrt(10000/sum(n)))
}