#' Map forest plot data
#'
#' Mapping functions to facilitate building forest objects from forest plot data
#'
#' @param x A data frame with tree records in rows and attributes in columns. Tree records can correspond to individual trees or groups of trees with an associated density.
#' @param y A data frame with shrub records in rows and attributes in columns. Records can correspond to individual shrubs (with crown dimensions and height) or groups of shrubs with an associated cover estimate.
#' @param mapping_x A named character vector to specify mappings of columns in \code{x} into attributes of \code{treeData} data frames. Accepted names (and the corresponding specifications for the columns in \code{x} are:
#' @param mapping_y A named character vector to specify mappings of columns in \code{y} into attributes of \code{shrubData} data frames. Accepted names (and the corresponding specifications for the columns in \code{y}) are:
#' \itemize{
#' \item{"Species": Species code (should follow codes in \code{SpParams}).}
#' \item{"Species.name": Species name. In this case, the species code will be drawn by matching names with species names in \code{SpParams}.}
#' \item{"N": Tree density (in ind./ha).}
#' \item{"Cover": Shrub cover (in \%).}
#' \item{"D1": Shrub largest crown diameter (in cm).}
#' \item{"D2": Shrub crown diameter orthogonal to the largest one (in cm).}
#' \item{"plot.size": Plot size (in m2) to which each record refers to. This is used to calculate tree density (stems per hectare) when not supplied or shrub cover when shrub data is given at the individual level.}
#' \item{"DBH": Diameter at breast height (in cm).}
#' \item{"Height": Tree or shrub height (in cm).}
#' \item{"Z50": Depth (in mm) corresponding to 50\% of fine roots.}
#' \item{"Z95": Depth (in mm) corresponding to 95\% of fine roots.}
#' }
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}) from which valid species names are drawn.
#' @param plot_size_x The size of tree plot sampled area (in m2). Alternatively, 'plot_size_x'
#' can be a column in \code{x} and specified in \code{mapping_x} to indicate that trees
#' have been measured in different subplots and, therefore, they represent different
#' densities per hectare.
#' @param plot_size_y The size of shrub plot sampled area (in m2). Alternatively, 'plot_size_y'
#' can be a column in \code{y} and specified in \code{mapping_y} to indicate that shrubs
#' have been measured in different subplots and, therefore, they represent different
#' cover values.
#'
#' @return Functions \code{forest_mapTreeTable} and \code{forest_mapShrubTable} return a data frame with the structure of \code{treeData} and \code{shrubData} from \code{\link{forest}} objects. Function \code{forest_mapWoodyTable} returns directly a \code{\link{forest}} object.
#'
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#'
#' @seealso \code{\link{forest}}, \code{\link{poblet_trees}}, \code{\link{forest_mergeTrees}}, \code{\link{tree2forest}}
#'
#' @name forest_mapWoodyTables
#'
#' @examples
#'
#' # Load species parameters
#' data(SpParamsMED)
#'
#' # Create an empty forest object
#' f <- emptyforest()
#'
#' # (1) Mapping tree data
#' # Load Poblet tree data
#' data(poblet_trees)
#'
#' # Subset control plot
#' x <- subset(poblet_trees, Plot.Code=="POBL_CTL")
#'
#' # Estimate sampled area (15-m radius plot)
#' sampled_area <- pi*15^2
#'
#' # Define mapping
#' mapping_x <- c("Species.name" = "Species", "DBH" = "Diameter.cm")
#'
#' # Map tree data for plot 'POBL_CTL'
#' f$treeData <- forest_mapTreeTable(x,
#'                     mapping_x = mapping_x, SpParams = SpParamsMED,
#'                     plot_size_x = sampled_area)
#'
#' # (2) Mapping shrub individual data
#' #
#' # Create the individual shrub data frame
#' species <- c("Erica arborea","Cistus albidus", "Erica arborea", "Cistus albidus", "Cistus albidus")
#' H <- c(200,50,100,40,30)
#' D1 <- c(140,40,100, 35,30)
#' D2 <- D1
#' y <- data.frame(species, H, D1, D2)
#'
#' # Define mapping (D1 and D2 map to variables with the same name)
#' mapping_y <- c("Species.name"= "species", "Height" ="H", "D1", "D2")
#'
#' # Map individual shrub data to cover data (here each individual becomes a cohort)
#' # assuming that the sampled area was 4 m2
#' f$shrubData <- forest_mapShrubTable(y,
#'                      mapping_y = mapping_y, SpParams = SpParamsMED,
#'                      plot_size_y = 4)
#'
#' # (3) Print forest attributes
#' summary(f, SpParamsMED)
#'
#' # (4) Forest initialization in a single step
#' f <- forest_mapWoodyTables(x, y,
#'                            mapping_x = mapping_x, mapping_y = mapping_y,
#'                            SpParams = SpParamsMED,
#'                            plot_size_x = sampled_area, plot_size_y = 4)
#' summary(f, SpParamsMED)
forest_mapTreeTable<-function(x, mapping_x, SpParams, plot_size_x = NULL) {
  n = nrow(x)
  treeData = data.frame(
    Species = rep(NA, n),
    N = rep(NA, n),
    Height = rep(NA, n),
    DBH = rep(NA, n),
    Z50 = rep(NA, n),
    Z95 = rep(NA, n))
  
  # Fill any empty mapping names with correspoinding mapping values
  names(mapping_x)[names(mapping_x)==""] = mapping_x[names(mapping_x)==""]
  
  if("Height" %in% names(mapping_x)) {
    treeData$Height = x[[mapping_x[["Height"]]]]
  }
  if("DBH" %in% names(mapping_x)) {
    treeData$DBH = x[[mapping_x[["DBH"]]]]
  }
  if("N" %in% names(mapping_x)) {
    treeData$N = x[[mapping_x[["N"]]]]
  } else {
    treeData$N = 1
  }
  if("Z50" %in% names(mapping_x)) {
    treeData$Z50 = x[[mapping_x[["Z50"]]]]
  }
  if("Z95" %in% names(mapping_x)) {
    treeData$Z95 = x[[mapping_x[["Z95"]]]]
  }
  if("Species" %in% names(mapping_x)) {
    treeData$Species = x[[mapping_x[["Species"]]]]
  }
  if("plot_size_x" %in% names(mapping_x)) {
    plot_size_x = x[[mapping_x[["plot_size_x"]]]]
  }
  if(!is.null(plot_size_x)) {
    treeData$N = treeData$N*(10000/plot_size_x)
  }
  if("Species.name" %in% names(mapping_x)) {
    non_recognized = character(0)
    Species.name = x[[mapping_x[["Species.name"]]]]
    for(i in 1:n) {
      indices = which(SpParams$Name==Species.name[i])
      if(length(indices)>0) {
        treeData$Species[i] = SpParams$Name[indices]
      } else {
        non_recognized = unique(c(non_recognized, Species.name[i]))
      }
    }
    if(length(non_recognized)>0) {
      warning(paste0("Taxon names that were not matched: ", paste0(non_recognized, collapse=","),"."))
    }
  }
  return(treeData)
}

#' @rdname forest_mapWoodyTables
#' @export
forest_mapShrubTable<-function(y, mapping_y, SpParams, plot_size_y = NULL) {
  n = nrow(y)
  shrubData = data.frame(
    Species = rep(NA, n),
    Height = rep(NA, n),
    Cover = rep(NA, n),
    Z50 = rep(NA, n),
    Z95 = rep(NA, n))
  
  # Fill any empty mapping names with correspoinding mapping values
  names(mapping_y)[names(mapping_y)==""] = mapping_y[names(mapping_y)==""]
  
  if("Height" %in% names(mapping_y)) {
    shrubData$Height = y[[mapping_y[["Height"]]]]
  }
  if("plot_size_y" %in% names(mapping_y)) {
    plot_size_y = y[[mapping_y[["plot_size_y"]]]]
  }
  if("Cover" %in% names(mapping_y)) {
    shrubData$Cover = y[[mapping_y[["Cover"]]]]
  } else if("D1" %in% names(mapping_y)) {
    if(is.null(plot_size_y)) stop("You must supply plot size when mapping crown diameter data.")
    D1 = y[[mapping_y[["D1"]]]] # D1 in cm
    D2 = D1
    if("D2" %in% names(mapping_y)) {
      D2 = y[[mapping_y[["D2"]]]] # D2 in cm
    }
    area = pi*((D1/200)*(D2/200)) #area in m2
    shrubData$Cover = pmin(100,100*area/plot_size_y)
  }
  if("Z50" %in% names(mapping_y)) {
    shrubData$Z50 = y[[mapping_y[["Z50"]]]]
  }
  if("Z95" %in% names(mapping_y)) {
    shrubData$Z95 = y[[mapping_y[["Z95"]]]]
  }
  if("Species" %in% names(mapping_y)) {
    shrubData$Species = y[[mapping_y[["Species"]]]]
  }
  if("Species.name" %in% names(mapping_y)) {
    non_recognized = character(0)
    Species.name = y[[mapping_y[["Species.name"]]]]
    for(i in 1:n) {
      indices = which(SpParams$Name==Species.name[i])
      if(length(indices)>0) {
        shrubData$Species[i] = SpParams$Name[indices]
      } else {
        non_recognized = unique(c(non_recognized, Species.name[i]))
      }
    }
    if(length(non_recognized)>0) {
      warning(paste0("Taxon names that were not matched: ", paste0(non_recognized, collapse=",")))
    }
  }
  return(shrubData)
}

#' @rdname forest_mapWoodyTables
#' @export
forest_mapWoodyTables<-function(x = NULL, y = NULL, mapping_x = NULL, mapping_y = NULL, SpParams, plot_size_x=NULL, plot_size_y = NULL) {
  f = emptyforest()
  if(!is.null(x)) {
    if(is.null(mapping_x)) stop("You need to specify a mapping for 'x'")
    f$treeData = forest_mapTreeTable(x, mapping_x = mapping_x,  SpParams=SpParams, plot_size_x = plot_size_x)
  }
  if(!is.null(y)) {
    if(is.null(mapping_y)) stop("You need to specify a mapping for 'y'")
    f$shrubData = forest_mapShrubTable(y, mapping_y = mapping_y,  SpParams=SpParams, plot_size_y = plot_size_y)
  }
  return(f)
}