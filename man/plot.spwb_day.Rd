\encoding{UTF-8}
\name{plot.spwb_day}
\alias{plot.growth_day}
\alias{plot.spwb_day}
\alias{plot.pwb_day}
\title{Plots simulation results for one day}
\description{
Functions to plot the subdaily simulation results of \code{\link{spwb_day}}, \code{\link{growth_day}} or the transpiration calculations of \code{\link{transp_transpirationSperry}}.
}
\usage{
\method{plot}{growth_day}(x, type="PlantTranspiration", bySpecies = FALSE, 
                        xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...)
\method{plot}{spwb_day}(x, type="PlantTranspiration", bySpecies = FALSE, 
                        xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...)
\method{plot}{pwb_day}(x, type="PlantTranspiration", bySpecies = FALSE, 
                       xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{spwb_day}, \code{growth_day} or \code{pwb_day}.}
  \item{type}{The information to be plotted (see details).}
  \item{bySpecies}{Allows aggregating output by species, before drawing plots. Aggregation can involve a sum (as for plant lai or transpiration) or a LAI-weighted mean (as for plant stress or plant water potential).}
  \item{xlim}{Range of values for x.}
  \item{ylim}{Range of values for y.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{...}{Additional parameters for function \code{plot}.}

}
\details{
The following plots are currently available for \code{spwb_day} and \code{pwb_day}:
  \itemize{
    \item{\code{"LeafPsi"}:}{Leaf water potential (for shade and sunlit leaves).}
    \item{\code{"LeafPsiAverage"}:}{Average leaf water potential.}
    \item{\code{"RootPsi"}:}{Root crown water potential.}
    \item{\code{"StemPsi"}:}{(Upper) stem water potential.}
    \item{\code{"StemPLC"}:}{(Average) percentage of loss conductance in the stem conduits.}
    \item{\code{"StemRWC"}:}{(Average) relative water content in the stem.}
    \item{\code{"LeafRWC"}:}{Relative water content in the leaf.}
    \item{\code{"StemSympRWC"}:}{(Average) relative water content in the stem symplasm.}
    \item{\code{"LeafSympRWC"}:}{Relative water content in the leaf symplasm.}
    \item{\code{"SoilPlantConductance"}:}{Overall soil plant conductance (calculated as the derivative of the supply function).}
    \item{\code{"PlantExtraction"}:}{ Water extracted from each soil layer.}
    \item{\code{"PlantTranspiration"}:}{ Plant cohort transpiration per ground area.}
    \item{\code{"TranspirationPerLeaf"}:}{ Plant cohort transpiration per leaf area.}
    \item{\code{"PlantGrossPhotosynthesis"}:}{ Plant cohort gross photosynthesis per ground area.}
    \item{\code{"GrossPhotosynthesisPerLeaf"}:}{ Plant cohort gross photosynthesis per leaf area.}
    \item{\code{"PlantNetPhotosynthesis"}:}{ Plant cohort net photosynthesis per ground area.}
    \item{\code{"NetPhotosynthesisPerLeaf"}:}{ Plant cohort net photosynthesis per leaf area.}
    \item{\code{"PlantAbsorbedSWR"}:}{ Absorbed short wave radiation per ground area  (differentiates sunlit and shade leaves).}
    \item{\code{"LeafTranspiration"}:}{ Instantaneous transpiration per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafGrossPhotosynthesis"}:}{ Instantaneous gross photosynthesis per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafNetPhotosynthesis"}:}{ Instantaneous net photosynthesis per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafAbsorbedSWR"}:}{ Absorbed short wave radiation per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafAbsorbedLWR"}:}{ Absorbed long wave radiation per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafCi"}:}{ Leaf intercellular CO2 concentration (differentiates sunlit and shade leaves).}
    \item{\code{"LeafIntrinsicWUE"}:}{ Leaf intrinsic water use efficiency, i.e. the ratio between instantaneous photosynthesis and stomatal conductance (differentiates sunlit and shade leaves).}
    \item{\code{"LeafVPD"}:}{ Leaf vapour pressure deficit (differentiates sunlit and shade leaves).}
    \item{\code{"LeafStomatalConductance"}:}{ Leaf stomatal conductance to water vapour (differentiates sunlit and shade leaves).}
    \item{\code{"LeafTemperature"}:}{ Leaf temperature (differentiates sunlit and shade leaves).}
    \item{\code{"Temperature"}:}{ Above-canopy, inside-canopy and soil temperature.}
    \item{\code{"CanopyEnergyBalance"}:}{ Canopy energy balance components.}
    \item{\code{"SoilEnergyBalance"}:}{ Soil energy balance components.}
    \item{\code{"PlantWaterBalance"}:}{ Difference between water extraction from the soil and transpired water per ground area.}
    \item{\code{"WaterBalancePerLeaf"}:}{ Difference between water extraction from the soil and transpired water per leaf area.}
   }

}
\note{
Only for soil plant water balance simulations using \code{transpirationMode = "Sperry"}. This function can be used to display subdaily dynamics of corresponding to single days on \code{\link{spwb}} runs, if control option \code{subdailyResults} is set to \code{TRUE}. See also option \code{subdaily} in \code{\link{plot.spwb}}.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{spwb_day}}, \code{\link{plot.spwb}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control = defaultControl()
control$ndailysteps = 24  

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2), W=c(0.5,0.5))


#Switch to 'Sperry' transpiration mode
control$transpirationMode="Sperry"

#Simulate one day only
x2 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
d = 100
sd2<-spwb_day(x2, examplesoil, rownames(examplemeteo)[d], 
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, 
             slope= 0, aspect = 0, prec = examplemeteo$Precipitation[d])

#Display transpiration for subdaily steps
plot(sd2, "PlantTranspiration")
}