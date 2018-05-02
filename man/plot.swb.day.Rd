\encoding{UTF-8}
\name{plot.swb.day}
\alias{plot.swb.day}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Displays simulation results for one day}
\description{
Function \code{plot} plots the results of the soil water balance model for one day (see \code{\link{swb.day}}).
}
\usage{
\method{plot}{swb.day}(x, type="PlantPsi", bySpecies = FALSE, xlab = NULL, ylab = NULL, ...)
}
\arguments{
  \item{x, object}{An object of class \code{swb}}
  \item{type}{The information to be plotted:
  \itemize{
    \item{\code{"PlantPsi"}:}{ Plant cohort water potential.}
    \item{\code{"PlantTranspiration"}:}{ Plant cohort transpiration per ground area.}
    \item{\code{"LeafTranspiration"}:}{ Plant cohort transpiration per leaf area.}
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

}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swb.day}}, \code{\link{plot.swb}}
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
x2 = forest2swbInput(exampleforest,examplesoil, SpParamsMED, control)
d = 100
sd2<-swb.day(x2, examplesoil, rownames(examplemeteo)[d], d , 
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, 
             slope= 0, aspect = 0, rain = examplemeteo$Precipitation[d])

#Display transpiration for subdaily steps
plot(sd2, "PlantTranspiration")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.