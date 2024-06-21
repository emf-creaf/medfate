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

#' Stand values
#' 
#' Functions to calculate stand attributes of a \code{\link{forest}} object.
#' 
#' @param x An object of class \code{\link{forest}}.
#' @param minDBH Minimum diameter at breast height (in cm) to include in estimation.
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
#' @param gdd Growth degree days (to account for leaf phenology effects).
#' @param includeDead A flag to indicate that standing dead fuels (dead branches) are included.
#' @param bounded A boolean flag to indicate that extreme values should be prevented (maximum tree LAI = 7 and maximum shrub LAI = 3)
#' 
#' @return 
#' \itemize{
#' \item{\code{stand_basalArea}: Stand basal area (m2/ha).}
#' \item{\code{stand_treeDensity}: Stand tree density (in ind/ha).}
#' \item{\code{stand_dominantTreeDiameter}: Dominant tree diameter, i.e the average diameter of the 100 widest trees (in cm).}
#' \item{\code{stand_quadraticMeanTreeDiameter}: Quadratic mean tree diameter, i.e. the diameter value corresponding to the current basal area and density.}
#' \item{\code{stand_meanTreeHeight}: Mean tree height (in cm).}
#' \item{\code{stand_dominantTreeHeight}: Dominant tree height, i.e the average height of the 100 tallest trees (in cm).}
#' \item{\code{stand_dominantTreeSpecies}: Dominant tree species name, determined in terms of basal area (and considering all tree sizes).}
#' \item{\code{stand_hartBeckingIndex}: Hart-Becking index.}
#' \item{\code{stand_foliarBiomass}: Standing biomass of leaves (in kg/m2).}
#' \item{\code{stand_fuel}: Stand fine fuel load (in kg/m2).}
#' \item{\code{stand_LAI}: Stand leaf area index (m2/m2).}
#' \item{\code{stand_shrubVolume}: Stand shrub phytovolume (m3/m2).}
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{forest}}, \code{\link{plant_basalArea}}, \code{\link{summary.forest}}
#' 
#' @examples
#' #Default species parameterization
#' data(SpParamsMED)
#'   
#' #Load example plot
#' data(exampleforest)
#'     
#' #A short way to obtain total basal area
#' stand_basalArea(exampleforest)
#'     
#' @name stand_values
#' @keywords internal
stand_dominantTreeDiameter<-function(x, minDBH = 7.5) {
  if(nrow(x$treeData)>0) {
    return(.dominantTreeDiameter(n = x$treeData$N, dbh = x$treeData$DBH, minDBH = minDBH))
  }
  return(as.numeric(NA))
}

#' @rdname stand_values
#' @keywords internal
stand_treeDensity<-function(x, minDBH = 7.5) {
  N <- as.numeric(NA)
  if(nrow(x$treeData)>0) N <- sum(x$treeData$N[x$treeData$DBH>minDBH])
  return(N)
}

#' @rdname stand_values
#' @keywords internal
stand_meanTreeHeight<-function(x, minDBH = 7.5) {
  MTH <- as.numeric(NA)
  if(nrow(x$treeData)>0) {
    h = x$treeData$Height[x$treeData$DBH>minDBH]
    n = x$treeData$N[x$treeData$DBH>minDBH]
    if(length(h)>0) {
      MTH <- sum(h*n)/sum(n)
    }
  }
  return(MTH)
}

#' @rdname stand_values
#' @keywords internal
stand_dominantTreeHeight<-function(x, minDBH = 7.5) {
  if(nrow(x$treeData)>0) {
    return(.dominantTreeHeight(n = x$treeData$N, h = x$treeData$Height, dbh = x$treeData$DBH, minDBH = minDBH))
  }
  return( as.numeric(NA))
}

#' @rdname stand_values
#' @keywords internal
stand_hartBeckingIndex<-function(x, minDBH = 7.5) {
  if(nrow(x$treeData)>0) {
    return(.hartBeckingIndex(n = x$treeData$N, h = x$treeData$Height, dbh = x$treeData$DBH, minDBH = minDBH))
  }
  return( as.numeric(NA))
}

#' @rdname stand_values
#' @keywords internal
stand_quadraticMeanTreeDiameter<-function(x, minDBH = 7.5) {
  if(nrow(x$treeData)>0) {
    return(.quadraticMeanTreeDiameter(n = x$treeData$N, dbh = x$treeData$DBH, minDBH = 7.5))
  }
  return( as.numeric(NA))
}

#' @rdname stand_values
#' @keywords internal
stand_dominantTreeSpecies <-function(x, SpParams) {
  s<-species_basalArea(x, SpParams)
  s<-s[!is.na(s)]
  s<-s[s>0.0]
  if(length(s)>0) {
    return(names(s)[which.max(s)])
  }
  return(as.character(NA))
}
