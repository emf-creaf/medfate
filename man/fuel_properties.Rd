\encoding{UTF-8}
\name{fuel_properties}
\alias{fuel_FCCS}
\alias{fuel_stratification}
\alias{fuel_windAdjustmentFactor}
\title{
Fuel stratification and fuel characteristics
}
\description{
Function \code{fuel_stratification} provides a stratification of the stand into understory and canopy strata. Function \code{fuel_FCCS} calculates fuel characteristics from a \code{forest} object following an adaptation of the protocols described for the Fuel Characteristics Classification System (Prichard et al. 2013). Function \code{fuel_windAdjustmentFactor} determines the adjustment factor of wind for surface fires, according to Andrews (2012). 
}
\usage{
fuel_stratification(object, SpParams, gdd = NA, mode = "MED",
                      heightProfileStep = 10.0, maxHeightProfile = 5000.0,
                      bulkDensityThreshold = 0.05)
fuel_FCCS(object, ShrubCover, CanopyCover, SpParams, cohortFMC = as.numeric(c()),
          gdd = NA, mode = "MED", heightProfileStep = 10, maxHeightProfile = 5000, 
          bulkDensityThreshold = 0.05, depthMode = "crownaverage")
fuel_windAdjustmentFactor(topShrubHeight, bottomCanopyHeight, topCanopyHeight,
                          canopyCover)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{\link{forest}}}
  \item{ShrubCover}{Total shrub cover (in percent) of the stand.}
  \item{CanopyCover}{Total canopy cover (in percent) of the stand.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}
  \item{cohortFMC}{A numeric vector of (actual) fuel moisture content by cohort (e.g. can be taken as a row of the matrix returned by \code{\link{moisture_cohortFMC}}).}
  \item{gdd}{Growth degree-days.}
  \item{mode}{Calculation mode, either "MED" or "US".}
  \item{heightProfileStep}{Precision for the fuel bulk density profile.}
  \item{maxHeightProfile}{Maximum height for the fuel bulk density profile.}
  \item{bulkDensityThreshold}{Minimum fuel bulk density to delimit fuel strata.}
  \item{depthMode}{Specifies how fuel depth (and therefore canopy and understory bulk density) should be estimated: 
     \itemize{
       \item{\code{"crownaverage"}: As weighed average of crown lengths using loadings as weights.}
       \item{\code{"profile"}: As the difference of base and top heights in bulk density profiles.}  
       \item{\code{"absoluteprofile"}: As the difference of absolute base and absolute top heights in bulk density profiles.}  
     }
  }
  \item{topShrubHeight}{Shrub stratum top height (in m).}
  \item{bottomCanopyHeight}{Canopy base height (in m).}
  \item{topCanopyHeight}{Canopy top height (in m).}
  \item{canopyCover}{Canopy percent cover.}
}
\details{
Details are described in a vignette.
}
\value{
Function \code{fuel_FCCS} returns a data frame with five rows corresponding to  fuel layers: \code{canopy}, \code{shrub}, \code{herb}, \code{woody} and \code{litter}. Columns correspond fuel properties:
\itemize{
\item{\code{w}: Fine fuel loading (in kg/m2).}
\item{\code{cover}: Percent cover.}
\item{\code{hbc}: Height to base of crowns (in m).}
\item{\code{htc}: Height to top of crowns (in m).}
\item{\code{delta}: Fuel depth (in m).}
\item{\code{rhob}: Fuel bulk density (in kg/m3).}
\item{\code{rhop}: Fuel particle density (in kg/m3).}
\item{\code{PV}: Particle volume (in m3/m2).}
\item{\code{beta}: Packing ratio (unitless).}
\item{\code{betarel}: Relative packing ratio (unitless).}
\item{\code{etabetarel}: Reaction efficiency (unitless).}
\item{\code{sigma}: Surface area-to-volume ratio (m2/m3).}
\item{\code{pDead}: Proportion of dead fuels.}
\item{\code{FAI}: Fuel area index (unitless).}
\item{\code{h}: High heat content (in kJ/kg).}
\item{\code{etaF}: Flammability modifier (between 1 and 2).}
\item{\code{RV}: Reactive volume (in m3/m2).}
\item{\code{MinFMC}: Minimum fuel moisture content (as percent over dry weight).}
\item{\code{MaxFMC}: Maximum fuel moisture content (as percent over dry weight).}
}

Function \code{fuel_stratification} returns a list with the following items:
\itemize{
\item{\code{surfaceLayerBaseHeight}: Base height of crowns of shrubs in the surface layer (in cm).}
\item{\code{surfaceLayerTopHeight}: Top height of crowns of shrubs in the surface layer (in cm).}
\item{\code{understoryLAI}: Cumulated LAI of the understory layer (i.e. leaf area comprised between surface layer base and top heights).}
\item{\code{canopyBaseHeight}: Base height of tree crowns in the canopy (in cm).}
\item{\code{canopyTopHeight}: Top height of tree crowns in the canopy (in cm).}
\item{\code{canopyLAI}: Cumulated LAI of the canopy (i.e. leaf area comprised between canopy base and top heights).}
}

Function \code{fuel_cohortFineFMC} returns a list with three matrices (for leaves, twigs and fine fuels). Each of them contains live moisture content values for each day (in rows) and plant cohort (in columns).

Function \code{fuel_windAdjustmentFactor} returns a value between 0 and 1.
}
\references{
Andrews, P. L. 2012. Modeling wind adjustment factor and midflame wind speed for Rothermel’s surface fire spread model. USDA Forest Service - General Technical Report RMRS-GTR:1–39.

Prichard, S. J., D. V Sandberg, R. D. Ottmar, E. Eberhardt, A. Andreu, P. Eagle, and K. Swedin. 2013. Classification System Version 3.0: Technical Documentation.

Reinhardt, E., D. Lutes, and J. Scott. 2006. FuelCalc: A method for estimating fuel characteristics. Pages 273–282.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}

\seealso{
\code{\link{fire_FCCS}}, \code{\link{spwb}}
}
\examples{
#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Show stratification of fuels
fuel_stratification(exampleforestMED, SpParamsMED)

#Calculate fuel properties according to FCCS
fccs = fuel_FCCS(exampleforestMED, 50,100, SpParamsMED)
fccs

fuel_windAdjustmentFactor(fccs$htc[2], fccs$hbc[1], fccs$htc[1], fccs$cover[1])
}
