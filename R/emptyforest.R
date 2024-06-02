#' Creation of an empty forest
#' 
#' Creates an empty \code{\link{forest}} object.
#' 
#' @param ntree,nshrub Number of tree and shrub cohorts, respectively.
#' @param nseed Number of species in the seed bank.
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
#' @name emptyforest
emptyforest <- function(ntree = 0, nshrub = 0, nseed = 0) {
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
  l$herbCover <- NA;
  l$herbHeight <- NA;
  l$seedBank <- data.frame(Species = as.character(rep(NA, nseed)),
                           Percent = as.numeric(rep(NA, nseed)))
  class(l)<-c("forest","list")
  return(l)
}

