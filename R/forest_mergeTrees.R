#' Forest merge functions
#' 
#' Functions to merge cohorts of a \code{\link{forest}} object.
#' 
#' @param x An object of class \code{\link{forest}}.
#' @param byDBHclass Logical flag to indicate that 5-cm tree DBH classes should be kept separated.
#' 
#' @details Tree DBH classes are defined in 5-cm intervals, whereas shrub height classes are defined in 10-cm intervals.
#' Tree DBH and shrub height classes are defined up to a specific size (i.e. larger plants are not merged) 
#' corresponding to 52.5 cm and 90 cm, respectively.
#' 
#' @return Another \code{\link{forest}} object with merged trees or shrubs, depending on the function.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{spwb}}, \code{\link{forest}},  \code{\link{fordyn}}, \code{\link{summary.forest}}
#' 
#' @name forest_mergeTrees
forest_mergeTrees<-function(x, byDBHclass = TRUE) {
  mergeTreesSize<-function(x) {
    ntree <- nrow(x)
    if(ntree>0) {
      BA <- x$N*pi*(x$DBH/200)^2
      BAsp <- tapply(BA, x$Species, FUN = sum)
      Nsp <- as.numeric(tapply(x$N, x$Species, FUN = sum))
      DBHsp <-  2*sqrt(10000*as.numeric(BAsp)/(pi*Nsp))
      y <- data.frame(Species = names(BAsp),
                     row.names = 1:length(BAsp),
                     stringsAsFactors = FALSE)
      y$DBH <- DBHsp 
      y$Height <- as.numeric(tapply(x$Height*BA, x$Species, FUN = sum)/BAsp)
      y$N <- Nsp 
      y$Z50 <- as.numeric(tapply(x$Z50*BA, x$Species, FUN = sum)/BAsp)
      y$Z95 <- as.numeric(tapply(x$Z95*BA, x$Species, FUN = sum)/BAsp)
      return(y)
    }
    return(x)
  }
  mergeTreesBySizeClass<-function(x) {
    sel0a <- (x$DBH >= 97.5) & (x$DBH < 102.5)
    sel0b <- (x$DBH >= 92.5) & (x$DBH < 97.5)
    sel1a <- (x$DBH >= 87.5) & (x$DBH < 92.5)
    sel1b <- (x$DBH >= 82.5) & (x$DBH < 87.5)
    sel2a <- (x$DBH >= 77.5) & (x$DBH < 82.5)
    sel2b <- (x$DBH >= 72.5) & (x$DBH < 77.5)
    sel3a <- (x$DBH >= 67.5) & (x$DBH < 72.5)
    sel3b <- (x$DBH >= 62.5) & (x$DBH < 67.5)
    sel4a <- (x$DBH >= 57.5) & (x$DBH < 62.5)
    sel4b <- (x$DBH >= 52.5) & (x$DBH < 57.5)
    sel5a <- (x$DBH >= 47.5) & (x$DBH < 52.5)
    sel5b <- (x$DBH >= 42.5) & (x$DBH < 47.5)
    sel6a <- (x$DBH >= 37.5) & (x$DBH < 42.5)
    sel6b <- (x$DBH >= 32.5) & (x$DBH < 37.5)
    sel7a <- (x$DBH >= 27.5) & (x$DBH < 32.5)
    sel7b <- (x$DBH >= 22.5) & (x$DBH < 27.5)
    sel8a <- (x$DBH >= 17.5) & (x$DBH < 22.5)
    sel8b <- (x$DBH >= 12.5) & (x$DBH < 17.5)
    sel9a <- (x$DBH >= 7.5) & (x$DBH < 12.5)
    sel9b <- (x$DBH >= 2.5) & (x$DBH < 7.5)
    sel10 <- (x$DBH < 2.5)
    nosel <- !(sel0a | sel0b | sel1a | sel1b | sel2a| sel2b | sel3a | sel3b | 
               sel4a | sel4b | sel5a | sel5b | sel6a | sel6b | sel7a | sel7b | 
               sel8a | sel8b | sel9a | sel9b | sel10)
    y <- x[nosel, c("Species", "DBH", "Height","N","Z50", "Z95"), drop = FALSE]
    if(sum(sel0a)>0) y <- rbind(y, mergeTreesSize(x[sel0a,, drop = FALSE]))
    if(sum(sel0b)>0) y <- rbind(y, mergeTreesSize(x[sel0b,, drop = FALSE]))
    if(sum(sel1a)>0) y <- rbind(y, mergeTreesSize(x[sel1a,, drop = FALSE]))
    if(sum(sel1b)>0) y <- rbind(y, mergeTreesSize(x[sel1b,, drop = FALSE]))
    if(sum(sel2a)>0) y <- rbind(y, mergeTreesSize(x[sel2a,, drop = FALSE]))
    if(sum(sel2b)>0) y <- rbind(y, mergeTreesSize(x[sel2b,, drop = FALSE]))
    if(sum(sel3a)>0) y <- rbind(y, mergeTreesSize(x[sel3a,, drop = FALSE]))
    if(sum(sel3b)>0) y <- rbind(y, mergeTreesSize(x[sel3b,, drop = FALSE]))
    if(sum(sel4a)>0) y <- rbind(y, mergeTreesSize(x[sel4a,, drop = FALSE]))
    if(sum(sel4b)>0) y <- rbind(y, mergeTreesSize(x[sel4b,, drop = FALSE]))
    if(sum(sel5a)>0) y <- rbind(y, mergeTreesSize(x[sel5a,, drop = FALSE]))
    if(sum(sel5b)>0) y <- rbind(y, mergeTreesSize(x[sel5b,, drop = FALSE]))
    if(sum(sel6a)>0) y <- rbind(y, mergeTreesSize(x[sel6a,, drop = FALSE]))
    if(sum(sel6b)>0) y <- rbind(y, mergeTreesSize(x[sel6b,, drop = FALSE]))
    if(sum(sel7a)>0) y <- rbind(y, mergeTreesSize(x[sel7a,, drop = FALSE]))
    if(sum(sel7b)>0) y <- rbind(y, mergeTreesSize(x[sel7b,, drop = FALSE]))
    if(sum(sel8a)>0) y <- rbind(y, mergeTreesSize(x[sel8a,, drop = FALSE]))
    if(sum(sel8b)>0) y <- rbind(y, mergeTreesSize(x[sel8b,, drop = FALSE]))
    if(sum(sel9a)>0) y <- rbind(y, mergeTreesSize(x[sel9a,, drop = FALSE]))
    if(sum(sel9b)>0) y <- rbind(y, mergeTreesSize(x[sel9b,, drop = FALSE]))
    if(sum(sel10)>0) y <- rbind(y, mergeTreesSize(x[sel10,, drop = FALSE]))
    return(y)
  }
  td <- x$treeData
  ntree <- nrow(td)
  x2 <- x
  if(ntree>0) {
    if(byDBHclass) {
      td2 <- mergeTreesBySizeClass(x$treeData)
    } else {
      BA <- td$N*pi*(td$DBH/200)^2
      BAsp <- tapply(BA, td$Species, FUN = sum)
      Nsp <- as.numeric(tapply(td$N, td$Species, FUN = sum))
      Hsp <- as.numeric(tapply(td$Height*BA, td$Species, FUN = sum)/BAsp)
      DBHsp <-  2*sqrt(10000*as.numeric(BAsp)/(pi*Nsp))
      Z50sp <- as.numeric(tapply(td$Z50*BA, td$Species, FUN = sum)/BAsp)
      Z95sp <- as.numeric(tapply(td$Z95*BA, td$Species, FUN = sum)/BAsp)
      td2 <- data.frame(Species = names(BAsp), 
                        DBH = DBHsp,
                        Height = Hsp, 
                        N = Nsp, 
                        Z50 = Z50sp, 
                        Z95 = Z95sp, 
                        row.names = 1:length(BAsp), stringsAsFactors = FALSE)
    }
    x2$treeData = td2
  }
  return(x2)
}

#' @rdname forest_mergeTrees
#' 
#' @param byHeightclass Boolean flag to indicate that 10-cm shrub height classes should be kept separated.
#' 
forest_mergeShrubs<-function(x, byHeightclass = TRUE) {
  mergeShrubsSize<-function(x) {
    nshrub = nrow(x)
    if(nshrub>0) {
      Coversp = tapply(x$Cover, x$Species, FUN = sum)
      Heightsp = tapply(x$Height*x$Cover, x$Species, FUN = sum)/Coversp
      Z50sp = tapply(x$Z50*x$Cover, x$Species, FUN = sum)/Coversp
      Z95sp = tapply(x$Z95*x$Cover, x$Species, FUN = sum)/Coversp
      y = data.frame(Species = names(Coversp),
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
    sel0a = (x$Height >= 210) & (x$Height < 220)
    sel0b = (x$Height >= 200) & (x$Height < 210)
    sel1a = (x$Height >= 180) & (x$Height < 190)
    sel1b = (x$Height >= 170) & (x$Height < 180)
    sel2a = (x$Height >= 160) & (x$Height < 170)
    sel2b = (x$Height >= 150) & (x$Height < 160)
    sel3a = (x$Height >= 140) & (x$Height < 150)
    sel3b = (x$Height >= 130) & (x$Height < 140)
    sel4b = (x$Height >= 120) & (x$Height < 130)
    sel4a = (x$Height >= 110) & (x$Height < 120)
    sel5a = (x$Height >= 100) & (x$Height < 110)
    sel5b = (x$Height >= 90) & (x$Height < 100)
    sel6a = (x$Height >= 80) & (x$Height < 90)
    sel6b = (x$Height >= 70) & (x$Height < 80)
    sel7a = (x$Height >= 60) & (x$Height < 70)
    sel7b = (x$Height >= 50) & (x$Height < 60)
    sel8a = (x$Height >= 40) & (x$Height < 50)
    sel8b = (x$Height >= 30) & (x$Height < 40)
    sel9a = (x$Height >= 20) & (x$Height < 30)
    sel9b = (x$Height >= 10) & (x$Height < 20)
    sel10 = (x$Height < 10)
    nosel <- !(sel0a | sel0b | sel1a | sel1b | sel2a| sel2b | sel3a | sel3b | 
                 sel4a | sel4b | sel5a | sel5b | sel6a | sel6b | sel7a | sel7b | 
                 sel8a | sel8b | sel9a | sel9b | sel10)
    y = x[nosel,, drop = FALSE]
    if(sum(sel0a)>0) y <- rbind(y, mergeShrubsSize(x[sel0a,, drop = FALSE]))
    if(sum(sel0b)>0) y <- rbind(y, mergeShrubsSize(x[sel0b,, drop = FALSE]))
    if(sum(sel1a)>0) y <- rbind(y, mergeShrubsSize(x[sel1a,, drop = FALSE]))
    if(sum(sel1b)>0) y <- rbind(y, mergeShrubsSize(x[sel1b,, drop = FALSE]))
    if(sum(sel2a)>0) y <- rbind(y, mergeShrubsSize(x[sel2a,, drop = FALSE]))
    if(sum(sel2b)>0) y <- rbind(y, mergeShrubsSize(x[sel2b,, drop = FALSE]))
    if(sum(sel3a)>0) y <- rbind(y, mergeShrubsSize(x[sel3a,, drop = FALSE]))
    if(sum(sel3b)>0) y <- rbind(y, mergeShrubsSize(x[sel3b,, drop = FALSE]))
    if(sum(sel4a)>0) y <- rbind(y, mergeShrubsSize(x[sel4a,, drop = FALSE]))
    if(sum(sel4b)>0) y <- rbind(y, mergeShrubsSize(x[sel4b,, drop = FALSE]))
    if(sum(sel5a)>0) y <- rbind(y, mergeShrubsSize(x[sel5a,, drop = FALSE]))
    if(sum(sel5b)>0) y <- rbind(y, mergeShrubsSize(x[sel5b,, drop = FALSE]))
    if(sum(sel6a)>0) y <- rbind(y, mergeShrubsSize(x[sel6a,, drop = FALSE]))
    if(sum(sel6b)>0) y <- rbind(y, mergeShrubsSize(x[sel6b,, drop = FALSE]))
    if(sum(sel7a)>0) y <- rbind(y, mergeShrubsSize(x[sel7a,, drop = FALSE]))
    if(sum(sel7b)>0) y <- rbind(y, mergeShrubsSize(x[sel7b,, drop = FALSE]))
    if(sum(sel8a)>0) y <- rbind(y, mergeShrubsSize(x[sel8a,, drop = FALSE]))
    if(sum(sel8b)>0) y <- rbind(y, mergeShrubsSize(x[sel8b,, drop = FALSE]))
    if(sum(sel9a)>0) y <- rbind(y, mergeShrubsSize(x[sel9a,, drop = FALSE]))
    if(sum(sel9b)>0) y <- rbind(y, mergeShrubsSize(x[sel9b,, drop = FALSE]))
    if(sum(sel10)>0) y <- rbind(y, mergeShrubsSize(x[sel10,, drop = FALSE]))
    return(y)
  }
  sd <- x$shrubData
  nshrub <- nrow(sd)
  x2 <- x
  if(nshrub>0) {
    if(byHeightclass) {
      sd2 <- mergeShrubsBySizeClass(x$shrubData)
    } else {
      sd2 <- mergeShrubsSize(x$shrubData)
    }
    x2$shrubData <- sd2
  }
  return(x2)
}
