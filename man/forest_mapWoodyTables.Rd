\name{forest_mapWoodyTables}
\alias{forest_mapWoodyTables}
\alias{forest_mapTreeTable}
\alias{forest_mapShrubTable}
\title{
Map forest plot data
}
\description{
Mapping functions to facilitate building forest objects from forest plot data 
}
\usage{
forest_mapWoodyTables(x, y, mapping, SpParams, plot.size = NULL)
forest_mapTreeTable(x, mapping, SpParams, plot.size = NULL)
forest_mapShrubTable(y, mapping, SpParams)
}
\arguments{
  \item{x}{A data frame with tree records in rows and attributes in columns.}
  \item{y}{A data frame with shrub records in rows and attributes in columns.}
  \item{mapping}{A named character vector to specify mappings of columns in \code{x} or \code{y} into attributes of \code{treeData} and \code{shrubData} data frames. Accepted names (and the corresponding specifications for the columns in \code{x} and \code{y}) are:
    \itemize{
      \item{"Species": Species code (should follow codes in \code{SpParams}).}
      \item{"Species.name": Species name. In this case, the species code will be drawn by matching names with species names in \code{SpParams}.}
      \item{"N": Tree density (in ind./ha).}
      \item{"plot.size": Plot size (in m2) to which each plot record refers to.}
      \item{"DBH": Diameter at breast height (in cm).}
      \item{"Height": Tree or shrub height (in cm).}
      \item{"Z50": Depth (in mm) corresponding to 50\% of fine roots.}
      \item{"Z95": Depth (in mm) corresponding to 95\% of fine roots.}
    }
  }
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{plot.size}{The size of plot sampled area (in m2). Alternatively, 'plot.size' can a column in \code{x} and specified in \code{mapping} to indicate that trees have been measured in different subplots and, therefore, they represent different densities per hectare.}
}
\value{
Functions \code{forest_mapTreeTable} and \code{forest_mapShrubTable} return a data frame with the structure of \code{treeData} and \code{shrubData} from \code{\link{forest}} objects. Function \code{forest_mapWoodyTable} returns directly a \code{\link{forest}} object.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
\code{\link{forest}}
}
