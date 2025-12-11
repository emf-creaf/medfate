#' Creation of an empty forest
#' 
#' Creates an empty \code{\link{forest}} object.
#' 
#' @param ntree,nshrub,nherb Number of tree, shrub and herb cohorts, respectively.
#' @param nseedling Number of species in the seedling bank.
#' @param nseed Number of species in the seed bank.
#' @param nlitter Number of items in the litter compartment.
#' @param addcolumns A character vector with additional columns. Currently allowed are (see \code{\link{forest}}):
#' \itemize{
#'  \item{\code{Z100}: A numeric vector with maximum rooting depth in mm.}
#'  \item{\code{LAI}: Leaf area index (m2/m2).}
#'  \item{\code{FoliarBiomass}: Standing dry biomass of leaves (kg/m2).}
#'  \item{\code{FuelLoading}: Fine fuel loading (kg/m2).}
#'  \item{\code{CrownRatio}: The ratio between crown length and total height (between 0 and 1)}
#'  \item{\code{Age}: A numeric vector indicating age of cohorts in years.}
#'  \item{\code{ObsID} A character vector to label specific cohorts in simulations of forest dynamics.}
#' }
#' 
#' @details
#' List elements \code{treeData} and \code{shrubData} are always created, regardless of the number of cohorts. 
#' In contrast, list elements \code{herbData}, \code{seedBank} and \code{seedlingBank} are only created if \code{nherb}, \code{nseed} and \code{nseedling} are non-zero, respectively.
#' 
#' 
#' @return An empty \code{\link{forest}} object.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{forest}}, \code{\link{tree2forest}}, \code{\link{summary.forest}}, \code{\link{forest_mapWoodyTables}},  \code{\link{forest_mergeTrees}},  
#' \code{\link{plot.forest}}
#' 
#' @examples 
#' # Initializes forest with 2 tree cohorts
#' emptyforest(ntree = 2)
#' 
#' # Initializes forest with 2 tree cohorts and 1 shrub cohort
#' emptyforest(ntree = 2, nshrub = 1)
#' 
#' # Initializes with extra column for leaf area index (LAI)
#' emptyforest(ntree = 2, nshrub = 1, addcolumns = "LAI")
#' 
#' # Initializes forest with 2 tree cohorts, 1 shrub cohort and 1 herbaceous cohort
#' emptyforest(ntree = 2, nshrub = 1, nherb = 1)
#' @name emptyforest
emptyforest <- function(ntree = 0, nshrub = 0, nherb = 0, nseedling = 0, nseed = 0, nlitter = 0,
                        addcolumns = NULL) {
  l <- list()
  l$treeData <- data.frame(Species=as.character(rep(NA, ntree)),
                           DBH=as.numeric(rep(NA, ntree)), 
                           Height=as.numeric(rep(NA, ntree)), 
                           N=as.numeric(rep(NA, ntree)),
                           Z50 = as.numeric(rep(NA, ntree)), 
                           Z95=as.numeric(rep(NA, ntree)))
  l$shrubData <- data.frame(Species=as.character(rep(NA, nshrub)), 
                            Height=as.numeric(rep(NA, nshrub)), 
                            Cover = as.numeric(rep(NA, nshrub)), 
                            Z50 = as.numeric(rep(NA, nshrub)), 
                            Z95=as.numeric(rep(NA, nshrub)))

  if(nherb>0) {
    l$herbData <- data.frame(Species=as.character(rep(NA, nshrub)), 
                              Height=as.numeric(rep(NA, nshrub)), 
                              Cover = as.numeric(rep(NA, nshrub)), 
                              Z50 = as.numeric(rep(NA, nshrub)), 
                              Z95=as.numeric(rep(NA, nshrub)))
  }
  if(nseedling > 0) {
    l$seedlingBank <- data.frame(Species = as.character(rep(NA, nseedling)),
                                 Percent = as.numeric(rep(NA, nseedling)),
                                 Age = as.numeric(rep(NA, nseedling)),
                                 Z50 = as.numeric(rep(NA, nseedling)),
                                 Z95 = as.numeric(rep(NA, nseedling)))
  }
  if(nseed > 0 ) {
    l$seedBank <- data.frame(Species = as.character(rep(NA, nseed)),
                             Percent = as.numeric(rep(NA, nseed)))
  }
  if(nlitter > 0) {
    l$litterData <- data.frame(Species = as.character(rep(NA, nlitter)),
                               Type = as.character(rep(NA, nlitter)),
                               Necromass = as.numeric(rep(NA, nlitter)))
  }
  
  if(!is.null(addcolumns)) {
    if(!is.character(addcolumns)) stop("`addcolumns` should be a character vector")
    addcolumns <- match.arg(addcolumns, c("Z100", "LAI", "FoliarBiomass", "FuelLoading", "CrownRatio", "Age", "ObsID") , several.ok = TRUE)
    if("Z100" %in% addcolumns) {
      l$treeData$Z100 <- as.numeric(rep(NA, ntree))
      l$shrubData$Z100 <- as.numeric(rep(NA, nshrub))
      if(nherb>0) l$herbData$Z100 <- as.numeric(rep(NA, nherb))
      if(nseedling>0) l$seedlingBank$Z100 <- as.numeric(rep(NA, nseedling))
    }
    if("Age" %in% addcolumns) {
      l$treeData$Age <- as.numeric(rep(NA, ntree))
      l$shrubData$Age <- as.numeric(rep(NA, nshrub))
      if(nherb>0) l$herbData$Age <- as.numeric(rep(NA, nherb))
    }
    if("LAI" %in% addcolumns) {
      l$treeData$LAI <- as.numeric(rep(NA, ntree))
      l$shrubData$LAI <- as.numeric(rep(NA, nshrub))
      if(nherb>0) l$herbData$LAI <- as.numeric(rep(NA, nherb))
    }
    if("FoliarBiomass" %in% addcolumns) {
      l$treeData$FoliarBiomass <- as.numeric(rep(NA, ntree))
      l$shrubData$FoliarBiomass <- as.numeric(rep(NA, nshrub))
      if(nherb>0) l$herbData$FoliarBiomass <- as.numeric(rep(NA, nherb))
    }
    if("FuelLoading" %in% addcolumns) {
      l$treeData$FuelLoading <- as.numeric(rep(NA, ntree))
      l$shrubData$FuelLoading <- as.numeric(rep(NA, nshrub))
      if(nherb>0) l$herbData$FuelLoading <- as.numeric(rep(NA, nherb))
    }
    if("CrownRatio" %in% addcolumns) {
      l$treeData$CrownRatio <- as.numeric(rep(NA, ntree))
      l$shrubData$CrownRatio <- as.numeric(rep(NA, nshrub))
      if(nherb>0) l$herbData$CrownRatio <- as.numeric(rep(NA, nherb))
    }
    if("ObsID" %in% addcolumns) {
      l$treeData$ObsID <- as.character(rep(NA, ntree))
      l$shrubData$ObsID <- as.character(rep(NA, nshrub))
      if(nherb>0) l$herbData$ObsID <- as.numeric(rep(NA, nherb))
    }
  }
  
  class(l)<-c("forest","list")
  return(l)
}

