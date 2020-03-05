\encoding{UTF-8}
\name{plot.spwb_day}
\alias{plot.spwb_day}
\alias{plot.pwb_day}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Displays simulation results for one day}
\description{
Functions to plot the results of the soil water balance model for one day (see \code{\link{spwb_day}}) or the transpiration calculations with Sperry (see \code{\link{transp_transpirationSperry}}).
}
\usage{
\method{plot}{spwb_day}(x, type="PlantTranspiration", bySpecies = FALSE, 
                        xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...)
\method{plot}{pwb_day}(x, type="PlantTranspiration", bySpecies = FALSE, 
                       xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{spwb_day} or \code{pwb_day}.}
  \item{type}{The information to be plotted:
  \itemize{
    \item{\code{"LeafPsi"}:}{Leaf water potential.}
    \item{\code{"StemPsi"}:}{(Upper) stem water potential.}
    \item{\code{"RootPsi"}:}{Root crown water potential.}
    \item{\code{"StemPLC"}:}{(Average) percentage of loss conductance in the stem conduits.}
    \item{\code{"StemRWC"}:}{(Average) relative water content in the stem symplasm.}
    \item{\code{"LeafRWC"}:}{Relative water content in the leaf symplasm.}
    \item{\code{"SoilPlantConductance"}:}{Overall soil plant conductance (calculated as the derivative of the supply function).}
    \item{\code{"PlantExtraction"}:}{ Water extracted from each soil layer.}
    \item{\code{"PlantTranspiration"}:}{ Plant cohort transpiration per ground area.}
    \item{\code{"PlantTranspirationPerLeaf"}:}{ Plant cohort transpiration per leaf area.}
    \item{\code{"PlantPhotosynthesis"}:}{ Plant cohort daily net photosynthesis per ground area.}
    \item{\code{"PlantPhotosynthesisPerLeaf"}:}{ Plant cohort daily net photosynthesis per leaf area.}
    \item{\code{"PlantWaterBalance"}:}{ Difference between water extraction from the soil and transpired water per ground area.}
    \item{\code{"PlantAbsorbedSWR"}:}{ Absorbed short wave radiation per ground area  (differentiates sunlit and shade leaves).}
    \item{\code{"LeafTranspiration"}:}{ Instantaneous transpiration per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafPhotosynthesis"}:}{ Instantaneous net photosynthesis per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafAbsorbedSWR"}:}{ Absorbed short wave radiation per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafVPD"}:}{ Leaf vapour pressure deficit (differentiates sunlit and shade leaves).}
    \item{\code{"LeafStomatalConductance"}:}{ Leaf stomatal conductance (differentiates sunlit and shade leaves).}
    \item{\code{"LeafTemperature"}:}{ Leaf temperature (differentiates sunlit and shade leaves).}
    \item{\code{"LeafCi"}:}{ Leaf intercellular CO2 concentration (differentiates sunlit and shade leaves).}
    \item{\code{"LeafIntrinsicWUE"}:}{ Leaf intrinsic water use efficiency, i.e. the ratio between instantaneous photosynthesis and stomatal conductance (differentiates sunlit and shade leaves).}
    \item{\code{"Temperature"}:}{ Above-canopy, inside-canopy and soil temperature.}
    \item{\code{"CanopyEnergyBalance"}:}{ Canopy energy balance components.}
    \item{\code{"SoilEnergyBalance"}:}{ Soil energy balance components.}
   }
  }
  \item{bySpecies}{Allows aggregating output by species, before drawing plots. Aggregation can involve a sum (as for plant lai or transpiration) or a LAI-weighted mean (as for plant stress or plant water potential).}
  \item{xlim}{Range of values for x.}
  \item{ylim}{Range of values for y.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{...}{Additional parameters for function \code{plot}.}

}
\note{
Only for soil plant water balance simulations using \code{transpirationMode = "Sperry"}. 
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