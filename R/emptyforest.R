#' Creation of an empty forest
#' 
#' Creates an empty \code{\link{forest}} object.
#' 
#' @param ntree,nshrub Number of tree and shrub cohorts, respectively.
#' @param nseedling Number of species in the seedling bank.
#' @param nseed Number of species in the seed bank.
#' @param addcolumns A character vector with additional columns. Currently allowed are:
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
#' @return An empty \code{\link{forest}} object.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{forest}}, \code{\link{tree2forest}}, \code{\link{summary.forest}}, \code{\link{forest_mapWoodyTables}},  \code{\link{forest_mergeTrees}},  
#' \code{\link{plot.forest}}
#' 
#' @examples 
#' # Initializes forest with 2 tree cohorts and 1 shrub cohort
#' emptyforest(ntree = 2, nshrub = 1)
#' 
#' # Initializes with extra column for leaf area index (LAI)
#' emptyforest(ntree = 2, nshrub = 1, addcolumns = "LAI")
#' 
#' @name emptyforest
emptyforest <- function(ntree = 0, nshrub = 0, nseedling = 0, nseed = 0, 
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
  l$herbCover <- NA
  l$herbHeight <- NA
  l$seedlingBank <- data.frame(Species = as.character(rep(NA, nseedling)),
                               Percent = as.numeric(rep(NA, nseedling)),
                               Age = as.numeric(rep(NA, nseedling)),
                               Z50 = as.numeric(rep(NA, nseedling)),
                               Z95 = as.numeric(rep(NA, nseedling)))
  l$seedBank <- data.frame(Species = as.character(rep(NA, nseed)),
                           Percent = as.numeric(rep(NA, nseed)))
  
  if(!is.null(addcolumns)) {
    if(!is.character(addcolumns)) stop("`addcolumns` should be a character vector")
    addcolumns <- match.arg(addcolumns, c("Z100", "LAI", "FoliarBiomass", "FuelLoading", "CrownRatio", "Age", "ObsID") , several.ok = TRUE)
    if("Z100" %in% addcolumns) {
      l$treeData$Z100 <- as.numeric(rep(NA, ntree))
      l$shrubData$Z100 <- as.numeric(rep(NA, nshrub))
      l$seedlingBank$Z100 <- as.numeric(rep(NA, nshrub))
    }
    if("Age" %in% addcolumns) {
      l$treeData$Age <- as.numeric(rep(NA, ntree))
      l$shrubData$Age <- as.numeric(rep(NA, nshrub))
    }
    if("LAI" %in% addcolumns) {
      l$treeData$LAI <- as.numeric(rep(NA, ntree))
      l$shrubData$LAI <- as.numeric(rep(NA, nshrub))
    }
    if("FoliarBiomass" %in% addcolumns) {
      l$treeData$FoliarBiomass <- as.numeric(rep(NA, ntree))
      l$shrubData$FoliarBiomass <- as.numeric(rep(NA, nshrub))
    }
    if("FuelLoading" %in% addcolumns) {
      l$treeData$FuelLoading <- as.numeric(rep(NA, ntree))
      l$shrubData$FuelLoading <- as.numeric(rep(NA, nshrub))
    }
    if("CrownRatio" %in% addcolumns) {
      l$treeData$CrownRatio <- as.numeric(rep(NA, ntree))
      l$shrubData$CrownRatio <- as.numeric(rep(NA, nshrub))
    }
    if("ObsID" %in% addcolumns) {
      l$treeData$ObsID <- as.character(rep(NA, ntree))
      l$shrubData$ObsID <- as.character(rep(NA, nshrub))
    }
  }
  
  class(l)<-c("forest","list")
  return(l)
}

