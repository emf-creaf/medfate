#' Summary of forest structure
#' 
#' Displays a summary of forest structure
#' 
#' @param object An object of class \code{\link{forest}}
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
#' @param x The object returned by \code{summary.forest}.
#' @param digits Minimal number of significant digits.
#' @param ... Additional parameters for functions \code{\link{summary}} and \code{\link{print}}.
#' 
#' @details Function \code{summary.forest} can be used to summarize a \code{forest} object in the console. 
#' 
#' @return Function \code{summary.forest} returns a list with several structural attributes, such as the basal area and LAI of the forest. 
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{forest}}, \code{\link{forest_mapWoodyTables}},  \code{\link{forest_mergeTrees}},  
#' \code{\link{plot.forest}}, \code{\link{tree2forest}}
#' 
#' @examples 
#' # Summary of example forest
#' summary(exampleforest, SpParamsMED)
#' 
#' @name summary.forest
summary.forest<-function(object, SpParams, ...) {

  # Checks
  if(length(object$herbCover)!=1) object$herbCover = NA
  if(length(object$herbHeight)!=1) object$herbHeight = NA
  
  ntree <- nrow(object$treeData)
  nshrub <- nrow(object$shrubData)
  
  selTree <- c(rep(TRUE, ntree), rep(FALSE, nshrub))
  selAdult <- c(object$treeData$DBH > 5.0, rep(FALSE, nshrub))
  selSapling <- c(object$treeData$DBH <= 5.0, rep(FALSE, nshrub))
  selShrub <- c(rep(FALSE, ntree), rep(TRUE, nshrub))

  coh_N <- plant_density(object, SpParams)
  coh_cov <- plant_cover(object, SpParams)
  coh_lai <- plant_LAI(object, SpParams)
  coh_fuel <- plant_fuelLoading(object, SpParams)
  
  coh_ba <- plant_basalArea(object, SpParams)
  
  s <- list()
  s["Tree_BA"] <- sum(coh_ba[selTree], na.rm=TRUE)
  s["Adult_BA"] <- sum(coh_ba[selAdult], na.rm=TRUE)
  s["Sapling_BA"] <- sum(coh_ba[selSapling], na.rm=TRUE)
  
  s["Tree_density"] <- sum(coh_N[selTree], na.rm=TRUE)
  s["Adult_density"] <- sum(coh_N[selAdult], na.rm=TRUE)
  s["Sapling_density"] <- sum(coh_N[selSapling], na.rm=TRUE)
  s["Shrub_density"] <- sum(coh_N[selShrub], na.rm=TRUE)


  s["Tree_cover"] <- pmin(100,sum(coh_cov[selTree], na.rm=TRUE))
  s["Adult_cover"] <- pmin(100,sum(coh_cov[selAdult], na.rm=TRUE))
  s["Sapling_cover"] <- pmin(100,sum(coh_cov[selSapling], na.rm=TRUE))
  s["Shrub_cover"] <- pmin(100,sum(coh_cov[selShrub], na.rm=TRUE))
  hc <- object$herbCover
  if(is.na(hc)) hc <- 0
  s["Herb_cover"] <- hc

  s["Tree_lai"] <- sum(coh_lai[selTree], na.rm=TRUE)
  s["Adult_lai"] <- sum(coh_lai[selAdult], na.rm=TRUE)
  s["Sapling_lai"] <- sum(coh_lai[selSapling], na.rm=TRUE)
  s["Shrub_lai"] <- sum(coh_lai[selShrub], na.rm=TRUE)
  s["Herb_lai"] <- herb_LAI(object, SpParams)
  s["Total_lai"] <- s[["Tree_lai"]] + s[["Shrub_lai"]] + s[["Herb_lai"]]
  
  s["Tree_fuel"] <- sum(coh_fuel[selTree], na.rm=TRUE)
  s["Adult_fuel"] <- sum(coh_fuel[selAdult], na.rm=TRUE)
  s["Sapling_fuel"] <- sum(coh_fuel[selSapling], na.rm=TRUE)
  s["Shrub_fuel"] <- sum(coh_fuel[selShrub], na.rm=TRUE)
  s["Herb_fuel"] <- herb_fuelLoading(object, SpParams)
  s["Total_fuel"] <- s[["Tree_fuel"]] + s[["Shrub_fuel"]] + s[["Herb_fuel"]]
  
  s["PARground"] <- NA
  s["SWRground"] <- NA
  
  if(all(!is.na(object$treeData$Height)) && all(!is.na(object$shrubData$Height))) {
    s["PARground"] <- light_PARground(object, SpParams)
    s["SWRground"] <- light_SWRground(object, SpParams)
  }
  class(s)<-c("summary.forest","list")
  return(s)
}

#' @rdname summary.forest
print.summary.forest<-function(x, digits=getOption("digits"),...) {
  cat(paste("Tree BA (m2/ha):", round(x[["Tree_BA"]],digits)," adult trees:", round(x[["Adult_BA"]], digits), 
            " saplings:", round(x[["Sapling_BA"]], digits),"\n"))
  cat(paste("Density (ind/ha) adult trees:", round(x[["Adult_density"]], digits), 
            " saplings:", round(x[["Sapling_density"]], digits), 
            " shrubs (estimated):", round(x[["Shrub_density"]],digits),"\n"))
  cat(paste("Cover (%) adult trees:", round(x[["Adult_cover"]],digits), 
            " saplings:", round(x[["Sapling_cover"]], digits), 
            " shrubs:", round(x[["Shrub_cover"]],digits),
            " herbs:", round(x[["Herb_cover"]],digits),
            "\n"))
  cat(paste("LAI (m2/m2) total:", round(x[["Total_lai"]], digits),
            " adult trees:", round(x[["Adult_lai"]], digits),
            " saplings:", round(x[["Sapling_lai"]], digits), 
            " shrubs:", round(x[["Shrub_lai"]], digits),
            " herbs:", round(x[["Herb_lai"]], digits),
            "\n"))
  cat(paste("Fuel loading (kg/m2) total:", round(x[["Total_fuel"]], digits),
            " adult trees:", round(x[["Adult_fuel"]], digits),
            " saplings:", round(x[["Sapling_fuel"]], digits), 
            " shrubs:", round(x[["Shrub_fuel"]], digits),
            " herbs:", round(x[["Herb_fuel"]],digits),
            "\n"))
  cat(paste("PAR ground (%):", round(x[["PARground"]],digits)," SWR ground (%):", round(x[["SWRground"]],digits),"\n"))
}
