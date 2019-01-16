\encoding{UTF-8}
\name{plot.spwb}
\alias{plot.spwb}
\alias{plot.growth}
\alias{summary.spwb}
\alias{summary.growth}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Displays simulation results}
\description{
Function \code{plot} plots the results of the soil water balance model (see \code{\link{spwb}}) or the forest growth model (see \code{\link{growth}}), whereas function \code{summary} summarizes the model's output in different temporal steps (i.e. weekly, annual, ...).
}
\usage{
\method{plot}{spwb}(x, type="PET_Precipitation", bySpecies = FALSE,
                   yearAxis=FALSE, xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                   add = FALSE, ...)
\method{plot}{growth}(x, type="PET_Precipitation", bySpecies = FALSE, yearAxis=FALSE, 
    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,...)
\method{summary}{spwb}(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...)
\method{summary}{growth}(object, freq="years", output="WaterBalance", FUN=sum, bySpecies = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, object}{An object of class \code{spwb}.}
  \item{type}{The information to be plotted:
  \itemize{
    \item{\code{"PET_Precipitation"}:}{ Potential evapotranspiration and Precipitation.}
    \item{\code{"PET_NetRain"}:}{ Potential evapotranspiration and Net rainfall.}
    \item{\code{"Snow"}:}{ Snow precipitation and snowpack dynamics.}
    \item{\code{"Export"}:}{ Water exported through deep drainage and surface runoff.}
    \item{\code{"ET"}:}{ Plant transpiration and soil evaporation.}
    \item{\code{"SoilPsi"}:}{ Soil water potential.}
    \item{\code{"SoilTheta"}:}{ Soil relative water content.}
    \item{\code{"SoilVol"}:}{ Soil water volumetric content.}
    \item{\code{"WTD"}:}{ Water table depth.}
    \item{\code{"LAI"}:}{ Expanded and dead leaf area index of the whole stand.}
    \item{\code{"PlantLAI"}:}{ Plant cohort leaf area index (expanded leaves).}
    \item{\code{"SoilPlantConductance"}:}{Average instantaneous overall soil plant conductance (calculated as the derivative of the supply function).}
    \item{\code{"PlantStress"}:}{ Plant cohort average daily drought stress.}
    \item{\code{"PlantPsi"}:}{ Plant cohort water potential (only for \code{transpirationMode = "Simple"}).}
    \item{\code{"LeafPsi"}:}{ Midday leaf water potential (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"StemPsi"}:}{ Midday (upper) stem water potential (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"RootPsi"}:}{ Midday root crown water potential (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"PlantTranspiration"}:}{ Plant cohort transpiration.}
    \item{\code{"PlantTranspirationLeaf"}:}{ Plant cohort transpiration per leaf area.}
    \item{\code{"PlantPhotosynthesis"}:}{ Plant cohort photosynthesis.}
    \item{\code{"PlantPhotosynthesisLeaf"}:}{ Plant cohort photosynthesis per leaf area.}
    \item{\code{"PlantAbsorbedSWR"}:}{ Plant cohort absorbed short wave radiation (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"PlantAbsorbedSWRLeaf"}:}{ Plant cohort absorbed short wave radiation per leaf area (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"PlantAbsorbedLWR"}:}{ Plant cohort absorbed long wave radiation (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"PlantAbsorbedLWRLeaf"}:}{ Plant cohort absorbed long wave radiation per leaf area (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"AirTemperature"}:}{ Minimum/maximum/mean daily temperatures above canopy (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"CanopyTemperature"}:}{ Minimum/maximum/mean daily temperatures inside canopy (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"SoilTemperature"}:}{ Minimum/maximum/mean daily temperatures inside the first soil layer (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"CanopyEnergyBalance"}:}{ Canopy energy balance components (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"SoilEnergyBalance"}:}{ Soil energy balance components (only for \code{transpirationMode = "Complex"}).}
    \item{\code{"PlantRespiration"}:}{ Plant cohort respiration (only for \code{plot.growth}).}
    \item{\code{"PlantRespirationLeaf"}:}{ Plant cohort respiration per leaf area (only for \code{plot.growth}).}
    \item{\code{"PlantCBalance"}:}{ Plant cohort carbon balance (only for \code{plot.growth}).}
    \item{\code{"PlantCBalanceLeaf"}:}{ Plant cohort carbon balance per leaf area (only for \code{plot.growth}).}
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
  \item{yearAxis}{A boolean to indicate whether the units of the x-axis are years (by default they are dates).}
  \item{xlim}{Range of values for x.}
  \item{ylim}{Range of values for y.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{add}{Boolean flag to add new elements to the current plot.}
  \item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
  \item{output}{The data table to be summarized. Accepted values are \code{"DailyBalance"}, \code{"PlantStress"}, \code{"PlantPsi"}, \code{"PlantTranspiration"}, \code{"PlantPhotosynthesis"}  and \code{"SoilWaterBalance"}, \code{"Temperature"} and \code{"EnergyBalance"}.}
  \item{FUN}{The function to summarize results (e.g., \code{sum}, \code{mean}, ...)}
  \item{...}{Additional parameters for functions \code{plot} or \code{summary}.}
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{spwb}}, \code{\link{spwbpoints}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params (2 layers)
examplesoil = soil(defaultSoilParams(2))

#Initialize control parameters
control = defaultControl()

#Initialize input
x = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
S1<-spwb(x, examplesoil, examplemeteo, elevation = 100)
        
#Plot results
plot(S1)

#Monthly summary (averages) of soil status
summary(S1, freq="months",FUN=mean, output="Soil")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.