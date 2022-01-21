.quadraticMeanTreeDiameter<-function(n, dbh){
  if(length(n)<1) return(NA)
  ba = sum(.treeBasalArea(n, dbh))
  return(2.0*sqrt(10000*ba/(pi*sum(n))))
}

.dominantTreeHeight<-function(n, h) {
  if(length(n)<1) return(NA)
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

.dominantTreeDiameter<-function(n, dbh) {
  if(length(n)<1) return(NA)
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
.hartBeckingIndex<-function(n,h) {
  return((10000/.dominantTreeHeight(n,h))*sqrt(10000/sum(n)))
}

stand_dominantTreeDiameter<-function(x) {
  return(.dominantTreeDiameter(n = x$treeData$N, dbh = x$treeData$DBH))
}
stand_dominantTreeHeight<-function(x) {
  return(.dominantTreeHeight(n = x$treeData$N, h = x$treeData$Height))
}
stand_hartBeckingIndex<-function(x) {
  return(.hartBeckingIndex(n = x$treeData$N, h = x$treeData$Height))
}
stand_quadraticMeanTreeDiameter<-function(x) {
  return(.quadraticMeanTreeDiameter(n = x$treeData$N, dbh = x$treeData$DBH))
}