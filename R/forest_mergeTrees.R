forest_mergeTrees<-function(x) {
  x2 = x
  td = x2$treeData
  ntree = nrow(td)
  if(ntree>0) {
    BA = plant_basalArea(x)[1:ntree]
    BAsp = tapply(BA, td$Species, FUN = sum)
    Nsp = as.numeric(tapply(td$N, td$Species, FUN = sum))
    Hsp = as.numeric(tapply(td$Height*BA, td$Species, FUN = sum)/BAsp)
    DBHsp =  2*sqrt(10000*as.numeric(BAsp)/(pi*Nsp))
    Z50sp = as.numeric(tapply(td$Z50*BA, td$Species, FUN = sum)/BAsp)
    Z95sp = as.numeric(tapply(td$Z95*BA, td$Species, FUN = sum)/BAsp)
    td2 = data.frame(Species = as.numeric(names(BAsp)), N = Nsp, DBH = DBHsp,
                     Height = Hsp, Z50 = Z50sp, Z95 = Z95sp, row.names = 1:length(BAsp), stringsAsFactors = FALSE)
    x2$treeData = td2
  }
  return(x2)
}