#' Single-cohort forests
#'
#' Creates a \code{\link{forest}} object with a single plant cohort
#'
#'
#' @param Species String with species (taxon) name or a non-negative integer for species identity (i.e., 0,1,2,...) matching SpParams.
#' @param Height Plant height (cm).
#' @param LAI Leaf area index (m2/m2)
#' @param N Tree density (ind/ha)
#' @param DBH Tree DBH (cm).
#' @param CrownRatio Crown ratio (fraction of total height)
#' @param FoliarBiomass Standing dry biomass of leaves (kg/m2)
#' @param FuelLoading Fine fuel loading (kg/m2)
#' @param Z50 Depth (in mm) corresponding to 50\% of fine roots.
#' @param Z95 Depth (in mm) corresponding to 95\% of fine roots.
#'
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{forest}}, \code{\link{forest_mergeTrees}},  \code{\link{plot.forest}}
#' @return An object of class \code{\link{forest}}
#' @export
#'
#' @name tree2forest
#' @examples
#' 
#' oak_forest <-tree2forest("Quercus ilex", Height= 200, LAI = 2)
#' 
tree2forest<-function(Species, Height, LAI = NA, N = NA, DBH = NA, Z50 = NA, Z95 = NA,
                      CrownRatio = NA, FoliarBiomass = NA, FuelLoading = NA) {
  f<- emptyforest(ntree = 1) 
  f$treeData[["Species"]] <- Species
  f$treeData[["Height"]] <- Height
  f$treeData[["N"]] <- N
  f$treeData[["DBH"]] <- DBH
  f$treeData[["Z50"]] <- Z50
  f$treeData[["Z95"]] <- Z95
  if(!is.na(LAI)) {
    f$treeData[["LAI"]] <- LAI
  }
  if(!is.na(CrownRatio)) {
    f$treeData[["CrownRatio"]] <- CrownRatio
  }
  if(!is.na(FoliarBiomass)) {
    f$treeData[["FoliarBiomass"]] <- FoliarBiomass
  }
  if(!is.na(FuelLoading)) {
    f$treeData[["FuelLoading"]] <- FuelLoading
  }
  return(f)
}

#' @rdname tree2forest
#' @param Cover Percent cover
#' @export
shrub2forest<-function(Species, Height, LAI = NA, Cover = NA, Z50 = NA, Z95 = NA,
                       CrownRatio = NA, FoliarBiomass = NA, FuelLoading = NA) {
  f <- emptyforest(nshrub = 1) 
  f$shrubData[["Species"]] <- Species
  f$shrubData[["Height"]] <- Height
  f$shrubData[["Cover"]] <- Cover
  f$shrubData[["Z50"]] <- Z50
  f$shrubData[["Z95"]] <- Z95
  if(!is.na(LAI)) {
    f$shrubData[["LAI"]] <- LAI
  }
  if(!is.na(CrownRatio)) {
    f$shrubData[["CrownRatio"]] <- CrownRatio
  }
  if(!is.na(FoliarBiomass)) {
    f$treeData[["FoliarBiomass"]] <- FoliarBiomass
  }
  if(!is.na(FuelLoading)) {
    f$treeData[["FuelLoading"]] <- FuelLoading
  }
  return(f)
}