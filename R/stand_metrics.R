.quadraticMeanTreeDiameter<-function(n, dbh, minDBH = 7.5){
  if(length(n)<1) return(NA)
  tba = .treeBasalArea(n, dbh)
  ba = sum(tba[dbh>=minDBH])
  return(2.0*sqrt(10000*ba/(pi*sum(n[dbh>=minDBH]))))
}

.dominantTreeHeight<-function(n, h, dbh, minDBH = 7.5) {
  if(length(n)<1) return(NA)
  o <-order(h, decreasing=TRUE)
  dbh = dbh[o]
  h = h[o]
  n = n[o]
  n = n[dbh>=minDBH]
  h = h[dbh>=minDBH]
  if(length(n)>0) {
    ncum = 0
    for(i in 1:length(h)) {
      ncum = ncum + n[i]
      if(ncum>100) return(sum(h[1:i]*n[1:i])/sum(n[1:i]))
    }
    return(sum(h*n)/sum(n))
  }
  return(NA)
}

.dominantTreeDiameter<-function(n, dbh, minDBH = 7.5) {
  if(length(n)<1) return(NA)
  o <-order(dbh, decreasing=TRUE)
  dbh = dbh[o]
  n = n[o]
  n = n[dbh>=minDBH]
  dbh = dbh[dbh>=minDBH]
  dtd = NA
  if(length(dbh)>0) {
    ncum = 0
    for(i in 1:length(dbh)) {
      ncum = ncum + n[i]
      if(ncum>100) return(sum(dbh[1:i]*n[1:i])/sum(n[1:i]))
    }
    dtd  = sum(dbh*n)/sum(n)
  }
  return(dtd)
}
.hartBeckingIndex<-function(n,h, dbh, minDBH = 7.5) {
  return((10000/.dominantTreeHeight(n,h, dbh, minDBH))*sqrt(10000/sum(n[dbh>=minDBH])))
}

stand_dominantTreeDiameter<-function(x, minDBH = 7.5) {
  return(.dominantTreeDiameter(n = x$treeData$N, dbh = x$treeData$DBH, minDBH = minDBH))
}
stand_dominantTreeHeight<-function(x, minDBH = 7.5) {
  return(.dominantTreeHeight(n = x$treeData$N, h = x$treeData$Height, dbh = x$treeData$DBH, minDBH = minDBH))
}
stand_hartBeckingIndex<-function(x, minDBH = 7.5) {
  return(.hartBeckingIndex(n = x$treeData$N, h = x$treeData$Height, dbh = x$treeData$DBH, minDBH = minDBH))
}
stand_quadraticMeanTreeDiameter<-function(x, minDBH = 7.5) {
  return(.quadraticMeanTreeDiameter(n = x$treeData$N, dbh = x$treeData$DBH, minDBH = 7.5))
}