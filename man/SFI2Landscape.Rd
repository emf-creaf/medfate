\encoding{UTF-8}
\name{SFI2SPL}
\alias{SFI2SPL}
\alias{SFI2SGL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Landscape Classes from SFI
}
\description{
Functions to create spatial landscape classes from Spanish Forest Inventory (SFI) data (see details). Internally, these functions make calls to \code{\link{extractSFIforest}}.
}
\usage{
SFI2SPL(SFItreeData, SFIshrubData, SFIherbData = NULL,
                    SpatialPointsIDs, elevation, slope, aspect, SpParams, 
                    SoilParamData = NULL, SFIcodes=NULL, control = defaultControl())
SFI2SGL(landTopo, SFItreeData, SFIshrubData, SpatialPointsIDs, 
                SpParams, SoilParamData = NULL, lctInput, 
                forestLCTs, shrublandLCTs, grasslandLCTs, agricultureLCTs= numeric(0), 
                rockLCTs= numeric(0), staticLCTs= numeric(0),
                SFIcodes=NULL, FMcodes = NULL, control = defaultControl())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{landTopo}{An object of class \code{\link{SpatialGridTopography}}.}
  \item{SFItreeData}{A data frame with measured tree data.}
  \item{SFIshrubData}{A data frame with measured shrub data.}
  \item{SFIherbData}{A data frame with percent cover and mean height of the herb layer.}
  \item{SpatialPointsIDs}{An object of class \code{\link{SpatialPoints-class}} containing the coordinates of the forest plots. Coordinates must include row names corresponding to plot IDs.}
  \item{elevation}{A numeric vector with elevation values (in m).}
  \item{slope}{A numeric vector with slope values (in degrees).}
  \item{aspect}{A numeric vector with aspect values (in degrees from North).}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{SoilParamData}{A data frame with soil parameters for each forest stand. The rows must match \code{SpatialPointsIDs} or \code{SpatialGridIDs}, depending on the function called. If \code{NULL}, then the parameters given by \code{\link{defaultSoilParams}} are used.}
  \item{SFIcodes}{A string vector (of length equal to the number of rows in \code{SpParams} of the SFI species codes that correspond to the model species codification. Each string may contain different coma-separated codes in order to merge SFI species into a single model species.}
  \item{FMcodes}{A string vector with the fuel model code corresponding to each forest plot. Names of this string vector should be plot IDs.}
  \item{control}{A list of control parameters (see \code{\link{defaultControl}}).}
  \item{lctInput}{A character vector with the land cover type input for initialization of unsurveyed cells.}
  \item{forestLCTs, shrublandLCTs, grasslandLCTs}{Integer vectors with the land cover types that are to be considered 'wildland' (i.e. forests, shrublands or grasslands). Initialization of forest and shrublands cells will make use of SFI plot data, but in shrublands only shrubs are used. Initialization of grasslands specifies an \code{\link{emptyforest}}.}
  \item{agricultureLCTs, rockLCTs, staticLCTs}{Integer vectors with the land cover types that are to be considered 'agriculture', 'rock' (i.e. outcrops) or 'static' (i.e. urban, water masses). Agricultural areas are like grasslands, but they can never be vegetated with wildland. Rock cells allow surface water flow over them, but do not contain vegetation nor soil and cannot be burned. Static cells do not burn and water entering them does not flow to any other cell.}
}
\value{
An object of class \code{\link{SpatialPointsLandscape-class}} or \code{\link{SpatialGridLandscape-class}} depending on the function.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}

\seealso{\code{\link{forest}}, \code{\link{extractSFIforest}}}
