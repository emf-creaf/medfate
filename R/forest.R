#' Forest description
#' 
#' Description of a forest stand
#' 
#' @param object An object of class \code{forest} has the following structure (see details):
#' \itemize{
#'   \item{\code{treeData}: A data frame of tree cohorts (in rows) and the following columns:
#'       \itemize{
#'         \item{\code{Species}: String with species (taxon) name or a non-negative integer for tree species identity (i.e., 0,1,2,...) matching SpParams.}
#'         \item{\code{Height}: Total tree height (in cm).}
#'         \item{\code{DBH}: Tree diameter at breast height (in cm).}
#'         \item{\code{N}: Density (number of individuals/hectare) that the measured tree represents.}
#'         \item{\code{Z50}: Depth (in mm) corresponding to 50\% of fine roots.}
#'         \item{\code{Z95}: Depth (in mm) corresponding to 95\% of fine roots.}
#'      }
#'   }
#'   \item{\code{shrubData}: A data frame of shrub cohorts (in rows) and the following columns:
#'       \itemize{
#'         \item{\code{Species}: String with species (taxon) name or a non-negative integer for shrub species identity (i.e., 0,1,2,...) matching SpParams.}
#'         \item{\code{Height}: Average total height of plants (in cm).}
#'         \item{\code{Cover}: Percent cover.}
#'         \item{\code{Z50}: Depth (in mm) corresponding to 50\% of fine roots.}
#'         \item{\code{Z95}: Depth (in mm) corresponding to 95\% of fine roots.}
#'       }
#'   }
#'   \item{\code{herbCover}: Percent cover of the herb layer (optional).}
#'   \item{\code{herbHeight}: Mean height (in cm) of the herb layer (optional).}
#'   \item{\code{seedBank}: A data frame containing seed bank information with the following columns:
#'       \itemize{
#'         \item{\code{Species}: String with species (taxon) name or a non-negative integer for tree species identity (i.e., 0,1,2,...) matching SpParams.}
#'         \item{\code{Percent}: Amount of seeds in relation to full seed bank (in \%).}
#'      }
#'   }
#' }
#' 
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
#' @param x The object returned by \code{summary.forest}.
#' @param digits Minimal number of significant digits.
#' @param ... Additional parameters for functions \code{\link{summary}} and \code{\link{print}}.
#' @param ntree,nshrub Number of tree and shrub cohorts, respectively.
#' @param nseed Number of species in the seed bank.
#' 
#' @details Function \code{summary.forest} can be used to summarize a \code{forest} object in the console. 
#' Function \code{emptyforest} creates an empty \code{forest} object.
#' 
#' The structure presented above for \code{forest} objects corresponds to the required data elements. 
#' A \code{forest} object can contain additional information when this is available. Data frames \code{treeData} 
#' and \code{shrubData} can contain additional columns:
#' \itemize{
#'   \item{\code{LAI}: Leaf area index (m2/m2)}
#'   \item{\code{FoliarBiomass}: Standing dry biomass of leaves (kg/m2)}
#'   \item{\code{FuelLoading}: Fine fuel loading (kg/m2)}
#'   \item{\code{CrownRatio}: The ratio between crown length and total height (between 0 and 1)}
#' }
#' Similarly, one can define \code{forest} list elements \code{herbLAI}, \code{herbFoliarBiomass} or \code{herbFuelLoading}.
#' All these values are used to override allometry-based estimates of those variables when initializing
#' inputs for functions \code{\link{spwb}} or \code{\link{spwb_day}}. Note that leaf area index, foliar biomass and
#' fuel loading are related entities, and they are treated as such in medfate. Therefore, users are expected to supply 
#' one or the other, and not all of them at the same time.
#' 
#' 
#' @return Function \code{summary.forest} returns a list with several structural attributes, such as the basal area and LAI of the forest. 
#' Function \code{emptyforest} returns an empty \code{forest} object.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{exampleforest}}, \code{\link{forest_mapWoodyTables}},  \code{\link{forest_mergeTrees}},  
#' \code{\link{plot.forest}}, \code{\link{tree2forest}}
#' 
#' @examples 
#' data(exampleforest)
#' data(SpParamsMED)
#' 
#' # Prints forest as a list of data items
#' exampleforest
#' 
#' # Summary of example forest
#' summary(exampleforest, SpParamsMED)
#' 
#' @name forest

#' @rdname forest
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

#' @rdname forest
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
    light_PARground(object, SpParams)
    light_SWRground(object, SpParams)
  }
  class(s)<-c("summary.forest","list")
  return(s)
}

#' @rdname forest
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
