\encoding{UTF-8}
\name{plot.spwb}
\alias{plot.spwb}
\alias{plot.pwb}
\alias{plot.growth}
\alias{plot.fordyn}

\title{Plots simulation results}
\description{
Function \code{plot} produces temporal graphs of the results of the soil plant water balance model (see \code{\link{spwb}}), plant water balance model (see \code{\link{pwb}}), the forest growth model (see \code{\link{growth}}) or the forest dynamics model (see \code{\link{fordyn}}).
}
\usage{
\method{plot}{spwb}(x, type="PET_Precipitation", 
                    cohorts = NULL, bySpecies = FALSE,
                    dates = NULL, subdaily = FALSE,
                    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                    summary.freq = NULL, ...)
\method{plot}{pwb}(x, type="PlantTranspiration", 
                   cohorts = NULL, bySpecies = FALSE,
                   dates = NULL, subdaily = FALSE,
                   xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                   summary.freq = NULL, ...)
\method{plot}{growth}(x, type="PET_Precipitation", 
                      cohorts = NULL, bySpecies = FALSE, 
                      dates = NULL, subdaily = FALSE,
                      xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                      summary.freq = NULL, ...)
\method{plot}{fordyn}(x, type="StandBasalArea",
                      ylim=NULL, xlab=NULL, ylab=NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{spwb}, \code{pwb}, \code{growth} or \code{fordyn}.}
  \item{type}{The information to be plotted (see details)}
  \item{cohorts}{An integer, boolean or character vector to select the plant cohorts to be plotted.}
  \item{bySpecies}{Allows aggregating output by species, before drawing plots (only has an effect with some values of \code{type}). Aggregation can involve a sum (as for plant lai or transpiration) or a LAI-weighted mean (as for plant stress or plant water potential).}
  \item{dates}{A Date vector with a subset of dates to be plotted.}
  \item{subdaily}{Whether subdaily results should be shown, only for simulations using \code{transpirationMode = "Sperry"} and having set \code{subdailyResults = TRUE} in the simulation control object. If \code{subdaily = TRUE}, then the valid strings for \code{type} are listed in \code{\link{plot.spwb_day}}.}
  \item{xlim}{Range of values for x.}
  \item{ylim}{Range of values for y.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{summary.freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
  \item{...}{Additional parameters for function \code{plot} (not used).}
}
\details{
The following plots are currently available for \code{\link{spwb}} (most of them also for \code{\link{pwb}}):
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
    \item{\code{"TranspirationPerLeaf"}:}{ Plant cohort transpiration per leaf area.}
    \item{\code{"PlantGrossPhotosynthesis"}:}{ Plant cohort photosynthesis.}
    \item{\code{"GrossPhotosynthesisPerLeaf"}:}{ Plant cohort photosynthesis per leaf area.}
    \item{\code{"PlantNetPhotosynthesis"}:}{ Plant cohort net photosynthesis (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"NetPhotosynthesisPerLeaf"}:}{ Plant cohort net photosynthesis per leaf area (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"PlantWUE"}:}{ Plant cohort daily water use efficiency (gross photosynthesis over transpiration; only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"PlantAbsorbedSWR"}:}{ Plant cohort absorbed short wave radiation (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"AbsorbedSWRPerLeaf"}:}{ Plant cohort absorbed short wave radiation per leaf area (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"PlantNetLWR"}:}{ Plant cohort net long wave radiation (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"NetLWRPerLeaf"}:}{ Plant cohort net long wave radiation per leaf area (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"AirTemperature"}:}{ Minimum/maximum/mean daily temperatures above canopy (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"CanopyTemperature"}:}{ Minimum/maximum/mean daily temperatures inside canopy (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"SoilTemperature"}:}{ Minimum/maximum/mean daily temperatures inside the first soil layer (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"CanopyEnergyBalance"}:}{ Canopy energy balance components (only for \code{transpirationMode = "Sperry"}).}
    \item{\code{"SoilEnergyBalance"}:}{ Soil energy balance components (only for \code{transpirationMode = "Sperry"}).}
  }
The following are only available for \code{\link{growth}}:
 \itemize{
    \item{\code{"GrossPhotosynthesis"}:}{ Gross photosynthesis rate per dry weight.}
    \item{\code{"MaintenanceRespiration"}:}{ Maintenance respiration cost per dry weight.}
    \item{\code{"CarbonBalance"}:}{ Carbon balance per dry weight.}
    \item{\code{"SugarLeaf"}:}{ Sugar concentration in leaves.}
    \item{\code{"StarchLeaf"}:}{ Starch concentration in leaves.}
    \item{\code{"SugarSapwood"}:}{ Sugar concentration in sapwood.}
    \item{\code{"StarchSapwood"}:}{ Starch concentration in sapwood.}
    \item{\code{"SugarTransport"}:}{ Phloem sugar transport rate.}
    \item{\code{"RootExudation"}:}{ Root exudation rate per dry weight.}
    \item{\code{"LeafPI0"}:}{ Leaf osmotic concentration at full turgor.}
    \item{\code{"StemPI0"}:}{ Sapwood osmotic concentration at full turgor.}
    \item{\code{"SapwoodArea"}:}{ Sapwood area per individual.}
    \item{\code{"LeafArea"}:}{ Leaf area per individual.}
    \item{\code{"FineRootArea"}:}{ Fine root area per individual.}
    \item{\code{"SapwoodBiomass"}:}{ Sapwood dry biomass per individual.}
    \item{\code{"LeafBiomass"}:}{ Leaf dry biomass per individual.}
    \item{\code{"FineRootBiomass"}:}{ Fine root dry biomass per individual.}
    \item{\code{"LabileBiomass"}:}{ Labile C biomass per individual.}
    \item{\code{"TotalLivingBiomass"}:}{ Total dry biomass per individual.}
    \item{\code{"SAgrowth"}:}{ Sapwood area growth rate.}
    \item{\code{"LAgrowth"}:}{ Leaf area growth rate.}
    \item{\code{"FRAgrowth"}:}{ Fine root area growth rate.}
    \item{\code{"HuberValue"}:}{ Ratio of leaf area to sapwood area.}
    \item{\code{"RootAreaLeafArea"}:}{ Ratio of fine root area to leaf area.}
   }
 The following are only available for \code{\link{fordyn}}:
 \itemize{
    \item{\code{"StandBasalArea"}:}{Stand basal area of living trees.}
    \item{\code{"StandLAI"}:}{Stand leaf area index.}
    \item{\code{"StandDensity"}:}{Stand density of living trees.}
    \item{\code{"SpeciesBasalArea"}:}{Basal area of living trees by species.}
    \item{\code{"SpeciesLAI"}:}{Leaf area index by species.}
    \item{\code{"SpeciesDensity"}:}{Density of living trees by species.}
    \item{\code{"CohortBasalArea"}:}{Basal area of living trees by plant cohort.}
    \item{\code{"CohortLAI"}:}{Leaf area index by plant cohort.}
    \item{\code{"CohortDensity"}:}{Density of living trees by plant cohort.}
   }
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
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
control = defaultControl("Granier")

#Initialize input
x = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)

#Call simulation function
S1<-spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
        
#Plot results
plot(S1)

}
