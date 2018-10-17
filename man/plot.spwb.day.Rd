\encoding{UTF-8}
\name{plot.spwb.day}
\alias{plot.spwb.day}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Displays simulation results for one day}
\description{
Function \code{plot} plots the results of the soil water balance model for one day (see \code{\link{spwb.day}}).
}
\usage{
\method{plot}{spwb.day}(x, type="PlantTranspiration", bySpecies = FALSE, xlab = NULL, ylab = NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{spwb}}
  \item{type}{The information to be plotted:
  \itemize{
    \item{\code{"LeafPsi"}:}{Leaf water potential.}
    \item{\code{"StemPsi"}:}{(Upper) stem water potential.}
    \item{\code{"RootPsi"}:}{Root crown water potential.}
    \item{\code{"StemPLC"}:}{(Average) percentage of loss conductance in the stem conduits.}
    \item{\code{"StemRWC"}:}{(Average) relative water content in the stem symplasm.}
    \item{\code{"LeafRWC"}:}{Relative water content in the leaf symplasm.}
    \item{\code{"PlantTranspiration"}:}{ Plant cohort transpiration per ground area.}
    \item{\code{"LeafTranspiration"}:}{ Plant cohort transpiration per leaf area.}
    \item{\code{"PlantWaterBalance"}:}{ Difference between water extraction from the soil and transpired water per ground area.}
    \item{\code{"PlantPhotosynthesis"}:}{ Plant cohort net photosynthesis per ground area.}
    \item{\code{"LeafPhotosynthesis"}:}{Net photosynthesis per leaf area.}
    \item{\code{"PlantAbsorbedSWR"}:}{Cohort's absorbed short wave radiation per ground area  (differentiates sunlit and shade leaves).}
    \item{\code{"LeafAbsorbedSWR"}:}{Cohort's absorbed short wave radiation per leaf area (differentiates sunlit and shade leaves).}
    \item{\code{"LeafVPD"}:}{Leaf vapour pressure deficit (differentiates sunlit and shade leaves).}
    \item{\code{"LeafStomatalConductance"}:}{Leaf stomatal conductance (differentiates sunlit and shade leaves).}
    \item{\code{"LeafTemperature"}:}{Leaf temperature (differentiates sunlit and shade leaves).}
    \item{\code{"Temperature"}:}{Above-canopy, inside-canopy and soil temperature.}
    \item{\code{"CanopyEnergyBalance"}:}{ Canopy energy balance components.}
    \item{\code{"SoilEnergyBalance"}:}{ Soil energy balance components.}
   }
  }
  \item{bySpecies}{Allows aggregating output by species, before drawing plots. Aggregation can involve a sum (as for plant lai or transpiration) or a LAI-weighted mean (as for plant stress or plant water potential).}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{...}{Additional parameters for function \code{plot}.}

}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{spwb.day}}, \code{\link{plot.spwb}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control = defaultControl()
control$ndailysteps = 24  

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2), W=c(0.5,0.5))


#Switch to 'Complex' transpiration mode
control$transpirationMode="Complex"

#Simulate one day only
x2 = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)
d = 100
sd2<-spwb.day(x2, examplesoil, rownames(examplemeteo)[d], 
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, 
             slope= 0, aspect = 0, prec = examplemeteo$Precipitation[d])

#Display transpiration for subdaily steps
plot(sd2, "PlantTranspiration")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.