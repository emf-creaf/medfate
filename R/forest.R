#' Forest description
#' 
#' Description of a forest stand
#' 
#' @param object An object of class \code{forest} has the following structure:
#' \itemize{
#'   \item{\code{treeData}: A data frame of tree cohorts (in rows) and the following columns:
#'       \itemize{
#'         \item{\code{Species}: Species name or a non-negative integer for tree species identity (i.e., 0,1,2,...) matching SpParams.}
#'         \item{\code{Height}: Total height (in cm).}
#'         \item{\code{DBH}: Diameter at breast height (in cm).}
#'         \item{\code{N}: Density (number of individuals/hectare).}
#'         \item{\code{Z50}: Depth (in mm) corresponding to 50\% of fine roots.}
#'         \item{\code{Z95}: Depth (in mm) corresponding to 95\% of fine roots.}
#'      }
#'   }
#'   \item{\code{shrubData}: A data frame of shrub cohorts (in rows) and the following columns:
#'       \itemize{
#'         \item{\code{Species}: Species name or a non-negative integer for shrub species identity (i.e., 0,1,2,...) matching SpParams.}
#'         \item{\code{Height}: Total height (in cm).}
#'         \item{\code{Cover}: Percent cover.}
#'         \item{\code{Z50}: Depth (in mm) corresponding to 50\% of fine roots.}
#'         \item{\code{Z95}: Depth (in mm) corresponding to 95\% of fine roots.}
#'       }
#'   }
#'   \item{\code{herbCover}: Percent cover of the herb layer.}
#'   \item{\code{herbHeight}: Mean height (in cm) of the herb layer.}
#' }
#' 
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
#' @param mode Calculation mode, either "MED" or "US".
#' @param x The object returned by \code{summary.forest}.
#' @param digits Minimal number of significant digits.
#' @param ... Additional parameters for functions \code{\link{summary}} and \code{\link{print}}.
#' @param ntree,nshrub Number of tree and shrub cohorts, respectively.
#' 
#' @details Function \code{summary.forest} can be used to summarize a \code{forest} object in the console. 
#' Function \code{emptyforest} creates an empty \code{forest} object.
#' 
#' @return Function \code{summary.forest} returns a list with several structural attributes, such as the basal area and LAI of the forest. 
#' Function \code{emptyforest} returns an empty \code{forest} object.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{exampleforestMED}}, \code{\link{forest_mergeTrees}},  \code{\link{plot.forest}}
#' 
#' @examples 
#' data(exampleforestMED)
#' data(SpParamsMED)
#' 
#' summary(exampleforestMED, SpParamsMED)
#' 
#' @name forest

#' @rdname forest
emptyforest <- function(ntree = 0, nshrub = 0) {
  l <- list()
  l$treeData <- data.frame(Species=character(ntree),DBH=numeric(ntree), 
                          Height=numeric(ntree), N=numeric(ntree),
                          Z50 = numeric(ntree), Z95=numeric(ntree))
  l$shrubData <- data.frame(Species=character(nshrub), Height=numeric(nshrub), 
                           Cover = numeric(nshrub), 
                           Z50 = numeric(nshrub), Z95=numeric(nshrub))
  l$herbCover <- 0;
  l$herbHeight <- 0;
  class(l)<-c("forest","list")
  return(l)
}

#' @rdname forest
summary.forest<-function(object, SpParams, mode = "MED", ...) {

  ntree <- nrow(object$treeData)
  nshrub <- nrow(object$shrubData)
  
  s <- list()
  s["BA"] <- stand_basalArea(object)
  s["N"] <- sum(plant_density(object, SpParams, mode = mode), na.rm=TRUE)
  coh_cov <- plant_cover(object, SpParams, mode = mode)
  covt <- 0
  covsh <- 0
  if(ntree>0) covt <- min(100,sum(coh_cov[1:ntree], na.rm=T))
  if(nshrub>0) covsh <- min(100,sum(coh_cov[(1+ntree):(1+ntree+nshrub)], na.rm=T))
  s["Tree_cover"] <- covt
  s["Shrub_cover"] <- covsh
  LAIc <- plant_LAI(object, SpParams, mode= mode)
  LAIt <- 0
  LAIsh <- 0
  if(ntree>0) LAIt <- sum(LAIc[1:ntree], na.rm=T)
  if(nshrub>0) LAIsh <- sum(LAIc[(1+ntree):(1+ntree+nshrub)], na.rm=T)
  s["LAI_trees"] <- LAIt
  s["LAI_shrubs"] <- LAIsh
  s["LAI"] <- LAIsh + LAIt
  Fuelc <- plant_fuel(object, SpParams, mode= mode)
  Fuelt <- 0
  Fuelsh <- 0
  if(ntree>0) Fuelt <- sum(Fuelc[1:ntree], na.rm=T)
  if(nshrub>0) Fuelsh <- sum(Fuelc[(1+ntree):(1+ntree+nshrub)], na.rm=T)
  s["Fuel_trees"] <- Fuelt
  s["Fuel_shrubs"] <- Fuelsh
  s["Fuel"] <- Fuelt + Fuelsh
  s["Phytovolume"] <- sum(plant_phytovolume(object, SpParams),na.rm=TRUE)
  s["PARground"] <- light_PARground(object, SpParams, mode = mode)
  s["SWRground"] <- light_SWRground(object, SpParams, mode = mode)
  class(s)<-c("summary.forest","list")
  return(s)
}

#' @rdname forest
print.summary.forest<-function(x, digits=getOption("digits"),...) {
  cat(paste("Tree density (ind/ha):", x[["N"]],"\n"))
  cat(paste("Tree BA (m2/ha):", round(x[["BA"]],digits),"\n"))
  cat(paste("Cover (%) trees (open ground):", round(x[["Tree_cover"]],digits), " shrubs:", round(x[["Shrub_cover"]],digits),"\n"))
  cat(paste("Shrub crown phytovolume (m3/m2):", round(x[["Phytovolume"]],digits),"\n"))
  cat(paste("LAI (m2/m2) total:", round(x[["LAI"]], digits)," trees:", round(x[["LAI_trees"]], digits),
            " shrubs:", round(x[["LAI_shrubs"]], digits),"\n"))
  cat(paste("Live fine fuel (kg/m2) total:", round(x[["Fuel"]], digits)," trees:", round(x[["Fuel_trees"]], digits),
            " shrubs:", round(x[["Fuel_shrubs"]], digits),"\n"))
  cat(paste("PAR ground (%):", round(x[["PARground"]],digits)," SWR ground (%):", round(x[["SWRground"]],digits),"\n"))
}