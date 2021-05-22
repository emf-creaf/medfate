forest_mergeTrees<-function(x, byDBHclass = TRUE) {
  mergeTreesSize<-function(x) {
    ntree = nrow(x)
    if(ntree>0) {
      BA = x$N*pi*(x$DBH/200)^2
      BAsp = tapply(BA, x$Species, FUN = sum)
      Nsp = as.numeric(tapply(x$N, x$Species, FUN = sum))
      DBHsp =  2*sqrt(10000*as.numeric(BAsp)/(pi*Nsp))
      y = data.frame(Species = as.numeric(names(BAsp)),
                     N = Nsp, DBH = DBHsp, 
                     row.names = 1:length(BAsp),
                     stringsAsFactors = FALSE)
      y$Height = as.numeric(tapply(x$Height*BA, x$Species, FUN = sum)/BAsp)
      y$Z50 = as.numeric(tapply(x$Z50*BA, x$Species, FUN = sum)/BAsp)
      y$Z95 = as.numeric(tapply(x$Z95*BA, x$Species, FUN = sum)/BAsp)
      return(y)
    }
    return(x)
  }
  mergeTreesBySizeClass<-function(x) {
    sel0a = (x$DBH >= 47.5) & (x$DBH < 52.5)
    sel0b = (x$DBH >= 42.5) & (x$DBH < 47.5)
    sel1a = (x$DBH >= 37.5) & (x$DBH < 42.5)
    sel1b = (x$DBH >= 32.5) & (x$DBH < 37.5)
    sel1c = (x$DBH >= 27.5) & (x$DBH < 32.5)
    sel1d = (x$DBH >= 22.5) & (x$DBH < 27.5)
    sel2a = (x$DBH >= 17.5) & (x$DBH < 22.5)
    sel2b = (x$DBH >= 12.5) & (x$DBH < 17.5)
    sel3a = (x$DBH >= 7.5) & (x$DBH < 12.5)
    sel3b = (x$DBH >= 2.5) & (x$DBH < 7.5)
    sel4 = (x$DBH < 2.5)
    nosel = !(sel0a | sel0b | sel1a | sel1b| sel1c| sel1d| sel2a| sel2b | sel3a | sel3b | sel4)
    y = x[nosel,, drop = FALSE]
    if(sum(sel0a)>0) y <- rbind(y, mergeTreesSize(x[sel0a,, drop = FALSE]))
    if(sum(sel0b)>0) y <- rbind(y, mergeTreesSize(x[sel0b,, drop = FALSE]))
    if(sum(sel1a)>0) y <- rbind(y, mergeTreesSize(x[sel1a,, drop = FALSE]))
    if(sum(sel1b)>0) y <- rbind(y, mergeTreesSize(x[sel1b,, drop = FALSE]))
    if(sum(sel1c)>0) y <- rbind(y, mergeTreesSize(x[sel1c,, drop = FALSE]))
    if(sum(sel1d)>0) y <- rbind(y, mergeTreesSize(x[sel1d,, drop = FALSE]))
    if(sum(sel2a)>0) y <- rbind(y, mergeTreesSize(x[sel2a,, drop = FALSE]))
    if(sum(sel2b)>0) y <- rbind(y, mergeTreesSize(x[sel2b,, drop = FALSE]))
    if(sum(sel3a)>0) y <- rbind(y, mergeTreesSize(x[sel3a,, drop = FALSE]))
    if(sum(sel3b)>0) y <- rbind(y, mergeTreesSize(x[sel3b,, drop = FALSE]))
    if(sum(sel4)>0) y <- rbind(y, mergeTreesSize(x[sel4,, drop = FALSE]))
    return(y)
  }
  td = x$treeData
  ntree = nrow(td)
  x2 = x
  if(ntree>0) {
    if(byDBHclass) {
      td2 = mergeTreesBySizeClass(x$treeData)
    } else {
      BA = plant_basalArea(x)[1:ntree]
      BAsp = tapply(BA, td$Species, FUN = sum)
      Nsp = as.numeric(tapply(td$N, td$Species, FUN = sum))
      Hsp = as.numeric(tapply(td$Height*BA, td$Species, FUN = sum)/BAsp)
      DBHsp =  2*sqrt(10000*as.numeric(BAsp)/(pi*Nsp))
      Z50sp = as.numeric(tapply(td$Z50*BA, td$Species, FUN = sum)/BAsp)
      Z95sp = as.numeric(tapply(td$Z95*BA, td$Species, FUN = sum)/BAsp)
      td2 = data.frame(Species = as.numeric(names(BAsp)), N = Nsp, DBH = DBHsp,
                       Height = Hsp, Z50 = Z50sp, Z95 = Z95sp, row.names = 1:length(BAsp), stringsAsFactors = FALSE)
    }
    x2$treeData = td2
  }
  return(x2)
}
forest_mergeShrubs<-function(x, byHeightclass = TRUE) {
  mergeShrubsSize<-function(x) {
    nshrub = nrow(x)
    if(nshrub>0) {
      Coversp = tapply(x$Cover, x$Species, FUN = sum)
      Heightsp = tapply(x$Height*x$Cover, x$Species, FUN = sum)/Coversp
      Z50sp = tapply(x$Z50*x$Cover, x$Species, FUN = sum)/Coversp
      Z95sp = tapply(x$Z95*x$Cover, x$Species, FUN = sum)/Coversp
      y = data.frame(Species = as.numeric(names(Coversp)),
                     Cover = as.numeric(Coversp),
                     Height = as.numeric(Heightsp),
                     Z50 = as.numeric(Z50sp),
                     Z95 = as.numeric(Z95sp),
                     row.names = 1:length(Coversp),
                     stringsAsFactors = FALSE)
      return(y)
    }
    return(x)
  }
  mergeShrubsBySizeClass<-function(x) {
    sel0a = (x$Height >= 80) & (x$Height < 90)
    sel0b = (x$Height >= 70) & (x$Height < 80)
    sel1a = (x$Height >= 60) & (x$Height < 70)
    sel1b = (x$Height >= 50) & (x$Height < 60)
    sel1c = (x$Height >= 40) & (x$Height < 50)
    sel1d = (x$Height >= 30) & (x$Height < 40)
    sel2a = (x$Height >= 20) & (x$Height < 30)
    sel2b = (x$Height >= 10) & (x$Height < 20)
    sel3 = (x$Height < 10)
    nosel = !(sel0a | sel0b | sel1a | sel1b| sel1c| sel1d| sel2a| sel2b | sel3)
    y = x[nosel,, drop = FALSE]
    if(sum(sel0a)>0) y <- rbind(y, mergeShrubsSize(x[sel0a,, drop = FALSE]))
    if(sum(sel0b)>0) y <- rbind(y, mergeShrubsSize(x[sel0b,, drop = FALSE]))
    if(sum(sel1a)>0) y <- rbind(y, mergeShrubsSize(x[sel1a,, drop = FALSE]))
    if(sum(sel1b)>0) y <- rbind(y, mergeShrubsSize(x[sel1b,, drop = FALSE]))
    if(sum(sel1c)>0) y <- rbind(y, mergeShrubsSize(x[sel1c,, drop = FALSE]))
    if(sum(sel1d)>0) y <- rbind(y, mergeShrubsSize(x[sel1d,, drop = FALSE]))
    if(sum(sel2a)>0) y <- rbind(y, mergeShrubsSize(x[sel2a,, drop = FALSE]))
    if(sum(sel2b)>0) y <- rbind(y, mergeShrubsSize(x[sel2b,, drop = FALSE]))
    if(sum(sel3)>0) y <- rbind(y, mergeShrubsSize(x[sel3,, drop = FALSE]))
    return(y)
  }
  sd = x$shrubData
  nshrub = nrow(sd)
  x2 = x
  if(nshrub>0) {
    if(byHeightclass) {
      sd2 = mergeShrubsBySizeClass(x$shrubData)
    } else {
      sd2 = mergeShrubsSize(x$shrubData)
    }
    x2$shrubData = sd2
  }
  return(x2)
}
