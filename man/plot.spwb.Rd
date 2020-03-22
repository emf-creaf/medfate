\encoding{UTF-8}
\name{plot.spwb}
\alias{plot.spwb}
\alias{plot.pwb}
\alias{plot.growth}
\alias{summary.pwb}
\alias{summary.spwb}
\alias{summary.growth}

\title{Displays simulation results}
\description{
Function \code{plot} plots the results of the soil plant water balance model (see \code{\link{spwb}}), plant water balance model (see \code{\link{pwb}}) or the forest growth model (see \code{\link{growth}}), whereas function \code{summary} summarizes the model's output in different temporal steps (i.e. weekly, annual, ...).
}
\usage{
\method{plot}{spwb}(x, type="PET_Precipitation", bySpecies = FALSE,
                    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...)
\method{plot}{pwb}(x, type="PlantTranspiration", bySpecies = FALSE,
                   xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...)                   
\method{plot}{growth}(x, type="PET_Precipitation", bySpecies = FALSE, 
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...)
\method{summary}{pwb}(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...)
\method{summary}{spwb}(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...)
\method{summary}{growth}(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...)
}
\arguments{
  \item{x, object}{An object of class \code{spwb}.}
  \item{type}{The information to be plotted:
  \itemize{
    \item{\code{"PET_Precipitation"}:}{ Potential evapotranspiration and Precipitation.}
    \item{\code{"PET_NetRain"}:}{ Potential evapotranspiration and Net rainfall.}
    \item{\code{"Snow"}:}{ Snow precipitation and snowpack dynamics.}
    \item{\code{"Export"}:}{ Water exported through deep drainage and surface runoff.}
    \item{\code{"Evapotranspiration"}:}{ Plant transpiration and soil evaporation.}
    \item{\code{"SoilPsi"}:}{ Soil water potential.}
    \item{\code{"SoilRWC"}:}{ Soil relative water content (in percent of field capacity).}
    \item{\code{"SoilTheta"}:}{ Soil moisture water content (in percent volume).}
    \item{\code{"SoilVol"}:}{ Soil water volumetric content (in mm).}
    \item{\code{"PlantExtraction"}:}{ Water extracted by plants from each soil layer.}
    \item{\code{"HydraulicRedistribution"}:}{Water added to each soil layer coming from other soil layers, transported through the plant hydraulic network (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"WTD"}:}{ Water table depth.}
    \item{\code{"LAI"}:}{ Expanded and dead leaf area index of the whole stand.}
    \item{\code{"PlantLAI"}:}{ Plant cohort leaf area index (expanded leaves).}
    \item{\code{"SoilPlantConductance"}:}{Average instantaneous overall soil plant conductance (calculated as the derivative of the supply function).}
    \item{\code{"PlantStress"}:}{ Plant cohort average daily drought stress.}
    \item{\code{"PlantPsi"}:}{ Plant cohort water potential (only for \code{transpirationMode = "Granier"}).}
    \item{\code{"LeafPsi"}:}{ Midday leaf water potential (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"StemPsi"}:}{ Midday (upper) stem water potential (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"RootPsi"}:}{ Midday root crown water potential (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"PlantTranspiration"}:}{ Plant cohort transpiration.}
    \item{\code{"PlantTranspirationPerLeaf"}:}{ Plant cohort transpiration per leaf area.}
    \item{\code{"PlantPhotosynthesis"}:}{ Plant cohort photosynthesis.}
    \item{\code{"PlantPhotosynthesisPerLeaf"}:}{ Plant cohort photosynthesis per leaf area.}
    \item{\code{"PlantWUE"}:}{ Plant cohort daily water use efficiency (photosynthesis over transpiration).}
    \item{\code{"PlantAbsorbedSWR"}:}{ Plant cohort absorbed short wave radiation (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"PlantAbsorbedSWRPerLeaf"}:}{ Plant cohort absorbed short wave radiation per leaf area (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"PlantAbsorbedLWR"}:}{ Plant cohort absorbed long wave radiation (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"PlantAbsorbedLWRPerLeaf"}:}{ Plant cohort absorbed long wave radiation per leaf area (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"AirTemperature"}:}{ Minimum/maximum/mean daily temperatures above canopy (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"CanopyTemperature"}:}{ Minimum/maximum/mean daily temperatures inside canopy (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"SoilTemperature"}:}{ Minimum/maximum/mean daily temperatures inside the first soil layer (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"CanopyEnergyBalance"}:}{ Canopy energy balance components (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"SoilEnergyBalance"}:}{ Soil energy balance components (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"PlantRespiration"}:}{ Plant cohort respiration (only for \code{plot.growth}).}
    \item{\code{"PlantRespirationPerLeaf"}:}{ Plant cohort respiration per leaf area (only for \code{plot.growth}).}
    \item{\code{"PlantRespirationPerIndividual"}:}{ Respiration per individual (only for \code{plot.growth}).}
    \item{\code{"PlantCBalance"}:}{ Plant cohort carbon balance (only for \code{plot.growth}).}
    \item{\code{"PlantCBalancePerLeaf"}:}{ Plant cohort carbon balance per leaf area (only for \code{plot.growth}).}
    \item{\code{"PlantCBalancePerIndividual"}:}{ Carbon balance per individual (only for \code{plot.growth}).}
    \item{\code{"PlantCstorageFast"}:}{ Amount of fast-dynamics carbon reserves  (only for \code{plot.growth}).}
    \item{\code{"PlantCstorageSlow"}:}{ Amount of slow-dynamics carbon reserves  (only for \code{plot.growth}).}
    \item{\code{"PlantSA"}:}{ Amount of sapwood area in an individual  (only for \code{plot.growth}).}
    \item{\code{"PlantSAgrowth"}:}{ Amount of newly-created sapwood area  (only for \code{plot.growth}).}
    \item{\code{"PlantRelativeSAgrowth"}:}{ Amount of newly-created sapwood area per sapwood area (only for \code{plot.growth}).}
    \item{\code{"PlantLAIlive"}:}{ Plant cohort leaf area index of live leaves (only for \code{plot.growth}).}
    \item{\code{"PlantLAIdead"}:}{ Plant cohort leaf area index of dead leaves (only for \code{plot.growth}).}
   }
  }
  \item{bySpecies}{Allows aggregating output by species, before calculating summaries or drawing plots (only has an effect with some values of \code{type}). Aggregation can involve a sum (as for plant lai or transpiration) or a LAI-weighted mean (as for plant stress or plant water potential).}
  \item{xlim}{Range of values for x.}
  \item{ylim}{Range of values for y.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
  \item{output}{The data table to be summarized. Accepted values are \code{"DailyBalance"}, \code{"PlantStress"}, \code{"PlantPsi"}, \code{"PlantTranspiration"}, \code{"PlantPhotosynthesis"}  and \code{"SoilWaterBalance"}, \code{"Temperature"} and \code{"EnergyBalance"}.}
  \item{FUN}{The function to summarize results (e.g., \code{sum}, \code{mean}, ...)}
  \item{...}{Additional parameters for function \code{plot} or \code{summary}.}
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{spwb}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl()

#Initialize input
x = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)

#Call simulation function
S1<-spwb(x, examplesoil, examplemeteo, elevation = 100)
        
#Plot results
plot(S1)

#Monthly summary (averages) of soil status
summary(S1, freq="months",FUN=mean, output="Soil")

}