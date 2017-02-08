\encoding{UTF-8}
\name{swbInput}
\alias{swbInput}
\alias{forest2aboveground}
\alias{forest2swbInput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Input variables for soil water balance
}
\description{
Function \code{forest2swbInput} takes an object of class \code{\link{forest}} to calculate input variables for  \code{\link{swb}}. Function \code{swbInput} does the same from input vectors. Function \code{forest2abovegroundInput} does the same as \code{forest2swbInput} but only calculates aboveground variables.
}
\usage{
forest2swbInput(x, soil, SpParams, control, gdd = NA)
forest2aboveground(x, SpParams, gdd = NA)
swbInput(SP, LAI, H, CR,  V,  soil, SpParams, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{soil}{An object of class \code{\link{soil}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsMED}}).}
  \item{control}{A list with default control parameters (see \code{\link{defaultControl}}).}
  \item{gdd}{Growth degree days to account for leaf phenology effects (in Celsius).}
  \item{SP}{An integer vector of plant cohort species identity (first species is 0).}
  \item{LAI}{A numeric vector (same length as \code{SP}) with plant cohort leaf area index (one-side leaf area relative to plot area).}
  \item{H}{A numeric vector (same length as \code{SP}) with plant cohort total heights (in cm).}
  \item{CR}{A numeric vector (same length as \code{SP}) with plant cohort crown ratio values (between 0 and 1).}
  \item{V}{A numeric matrix (with as many columns as soil layers and as many rows as the length as \code{SP}) containing the proportion of roots of each plant in each soil layer.}
  
}
\details{
Functions \code{forest2swbInput} and \code{forest2abovegroundInput} extracts height and species identity from plant cohorts of \code{x}, and calculate leaf area index and crown ratio.\code{forest2swbInput} also calculates the distribution of fine roots across soil. Both \code{forest2swbInput} and \code{swbInput} find parameter values for each plant cohort according to the parameters of its species as specified in \code{SpParams}. If \code{control$transpirationMode = "Sperry"} the functions also estimate the maximum conductance of rhizosphere, root xylem and stem xylem elements.
}
\value{
Function \code{forest2aboveground} returns a data frame with the following columns:
\itemize{
  \item{\code{SP}: Species identity (an integer).}
  \item{\code{N}: Cohort density (ind/ha) (see function \code{\link{plant.Density}}).}
  \item{\code{DBH}: Tree diameter at breast height (cm).}
  \item{\code{H}: Plant height (cm).}
  \item{\code{CR}: Crown ratio (crown length to total height).}
  \item{\code{LAI}: Leaf area index (m2/m2).}
  \item{\code{LA_live}: Live leaf area per individual (m2/ind).}
  \item{\code{LA_dead}: Dead leaf area per individual (m2/ind).}
  \item{\code{SA}: Live sapwood area per individual (cm2/ind).}
}

Functions \code{forest2swbInput} and \code{swbInput} return list of class \code{swbInput} with the following elements:
\itemize{
  \item{\code{verbose}: Boolean flag to indicate extra console output.}
  \item{\code{TranspirationMode}: Transpiration mode (either 'Simple' or 'Sperry').}
  \item{\code{above}: A data frame with columns \code{SP}, \code{H}, \code{CR} and \code{LAI} (see function \code{forest2aboveground}).}
  \item{\code{V}: A matrix with the proportion of fine roots of each cohort (in rows) in each soil layer (in columns).}
  \item{\code{k}: PAR extinction coefficient.}
  \item{\code{g}: Canopy water retention capacity per LAI unit (mm/LAI).}
  \item{\code{CR}: Crown ratio (crown length to total height).}
  \item{\code{Sgdd}: Growth degree days needed for leaf budburst (for winter decideous species).}
  \item{\code{Transpiration}: Plant cohort transpiration of the current day (mm of water = L/m2; filled with zeroes before simulations).}
  \item{\code{Photosynthesis}: Plant cohort photosynthesis of the current day (g C/m2; filled with zeroes before simulations).}
}
If \code{control$transpirationMode = "Simple"} there is an additional element:
\itemize{
  \item{\code{Psi_Extract}: Water potential corresponding to 50\% relative conductance (in MPa).}
  \item{\code{WUE}: Water use efficiency for carbon assimilation (g C/mm water).}
}
If \code{control$transpirationMode = "Sperry"} there are the following additional elements:
\itemize{
  \item{\code{Gwmax}: Maximum stomatal conductance to water vapor (in mol H2O·m-2·s-1).}
  \item{\code{VGrhizo_kmax}: A matrix with maximum rhizosphere conductance values of each cohort (in rows) in each soil layer (in columns).}
  \item{\code{VGroot_kmax}: A matrix with maximum root xylem conductance values of each cohort (in rows) in each soil layer (in columns).}
  \item{\code{VGstem_kmax}: Maximum stem xylem conductance values of each cohort.}
  \item{\code{VCstem_c}, \code{VCstem_d}: Parameters of the root xylem vulnerability curve.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{\code{\link{swb}}, \code{\link{soil}},  \code{\link{forest}}, \code{\link{SpParamsMED}}, \code{\link{defaultSoilParams}}}

\examples{
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Aboveground parameters
forest2aboveground(exampleforest, SpParamsMED)

#Initialize soil with default soil params
examplesoil = soil(defaultSoilParams())

#Initialize control parameters
control = defaultControl()

#Prepare input for simple hydraulic mode
forest2swbInput(exampleforest,examplesoil,SpParamsMED, control)

#Prepare input for Sperry hydraulic mode
control$transpirationMode = "Sperry"

forest2swbInput(exampleforest,examplesoil,SpParamsMED, control)
}
