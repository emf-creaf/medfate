woodformation_relativeGrowthRate<-function(dbh1, dbh2, yeardiff, lower = -2, upper = 8){
  grow<-function(dbh1, x, ny = 10) {
    dbh = dbh1
    for(i in 1:ny) {
      b <- (dbh/2)^2+dbh*x
      if(b<0) return(0)
      dbh = 2*sqrt(b)
    }
    return(dbh)
  }
  r<-uniroot(function(x) dbh2 - grow(dbh1,x,yeardiff), lower=lower, upper = upper)
  return(r$root)
}