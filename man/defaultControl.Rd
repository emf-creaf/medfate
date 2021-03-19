\encoding{UTF-8}
\name{defaultControl}
\alias{defaultControl}

\title{Default control parameters for models}
\description{
Creates a list with global default parameters for simulation models.
}
\usage{
defaultControl(transpirationMode = "Granier") 
}
\arguments{
  \item{transpirationMode}{Transpiration model (either 'Granier' or 'Sperry'). See \code{\link{spwbInput}}.}
}
\details{
The function returns a list with default parameters. Users can change those defaults that need to be set to other values and use the list as input for model functions. The relevant parameters are different for each model function.
}
\value{
A list, with the following options (default values in parentheses):
\itemize{
\bold{General}:
\itemize{
   \item{\code{verbose (=TRUE)}: Boolean flag to indicate console output during calculations.}
   \item{\code{modifyInput (=TRUE)}: Boolean flag to indicate that simulations will modify input object. If set to FALSE, simulations will not modify the input R object but return the current (modified) state variables within the output.}
   \item{\code{subdailyResults (=FALSE)}: Boolean flag to force subdaily results to be stored (as a list called 'subdaily' of \code{\link{spwb_day}} objects, one by simulated date) in calls to \code{\link{spwb}}.}
   \item{\code{soilFunctions (="SX")}: Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).}
   \item{\code{defaultWindSpeed (= 2.5)}: Default wind speed value (in m/s) to be used when missing from data. }
}
\bold{Water balance}:
\itemize{
  \item{\code{snowpack (=TRUE)}: Boolean flag to indicate the simulation of snow accumulation and melting.}
  \item{\code{rockyLayerDrainage (=TRUE)}: Boolean flag to indicate the simulation of drainage from rocky layers (> 95\% of rocks).}
  \item{\code{unlimitedSoilWater (=FALSE)}: Boolean flag to indicate the simulation of plant transpiration assuming that soil water is always at field capacity.}
  \item{\code{leafPhenology (=TRUE)}: Boolean flag to indicate the simulation of leaf phenology for winter-deciduous species.}
  \item{\code{unfoldingDD (=300)}: Degree-days for complete leaf unfolding after budburst has occurred.}
  \item{\code{plantWaterPools (=FALSE)}: Boolean flag to indicate the simulation of water balance assuming that each cohort has its own water pools.}
  \item{\code{transpirationMode (="Granier")}: Transpiration model (either 'Granier' or 'Sperry'). See \code{\link{spwbInput}}.}
  \item{\code{verticalLayerSize (= 100)}: Size of vertical layers (in cm) for photosynthesis calculation.}
  \item{\code{boundaryLayerSize (= 2000)}: Size of the boundary layer (in cm) over the canopy (relevant for multi-layer canopy energy balance).}
  \item{\code{windMeasurementHeight (=200)}: Height (in cm) over the canopy corresponding to wind measurements.}
  \item{\code{cavitationRefill (= "total")}: A string indicating how refilling of embolized conduits is done ("none" - no refilling, "annual" - every first day of the year, "rate" - following a rate of new sapwood formation, "total" - instantaneous complete refilling).}
  \item{\code{refillMaximumRate (= 0.05)}: Maximum rate of daily refilling of embolized conduits as sapwood area per leaf area (in cm2·m-2·day-1).}
}
\bold{Water balance with 'Sperry'}:
\itemize{
   \item{\code{gainModifier, costModifier (= 1)}: Modifiers (exponents) of the gain and cost functions defined in Sperry et al. (2016).}
   \item{\code{cochard (=FALSE)}: Boolean flag to indicate the stomatal closure when leaf water potential is below turgor loss point so that transpiration depends on cuticular properties (Martin-StPaul et al. 2017).}
   \item{\code{ndailysteps (= 24)}: Number of steps into which each day is divided for determination of stomatal conductance, transpiration and photosynthesis (24 equals 1-hour intervals).}
   \item{\code{nsubsteps (= 3600)}: Number of substeps into which each step is divided for multi-layer canopy energy balance solving.}
   \item{\code{capacitance (=FALSE)}: Whether the effect of plant water compartments is considered in simulations.}
   \item{\code{klatstem (= 0.01)}: Stem symplastic-apoplastic lateral conductance (in mmol·s-1·m-2·MPa-1). Only relevant when \code{capacitance = TRUE}.}
   \item{\code{klatleaf (= 0.01)}: Leaf symplastic-apoplastic lateral conductance (in mmol·s-1·m-2·MPa-1). Only relevant when \code{capacitance = TRUE}.}
   \item{\code{taper (= TRUE)}: Whether taper of xylem conduits is accounted for when calculating aboveground stem conductance from xylem conductivity.}
   \item{\code{numericParams}: A list with the following elements:
      \itemize{
        \item{\code{maxNsteps (= 400)}: Maximum number of steps in supply function.}
        \item{\code{ntrial (= 200)}: Number of iteration trials when finding root of equation system.}
        \item{\code{psiTol (= 0.0001)}: Tolerance value for water potential.}
        \item{\code{ETol (= 0.0001)}: Tolerance value for flow.}
      }
   }
   \item{\code{thermalCapacityLAI (=1000000)}: Thermal canopy capacitance per LAI unit.}
   \item{\code{Catm (=386)}: Atmospheric CO2 concentration (in micromol·mol-1 = ppm).}
   \item{\code{fracLeafResistance (= NA)}: Fraction of plant total resistance (leaf+stem+root) that corresponds to leaves. This fraction is used if \code{VCleaf_kmax = NA}.}
   \item{\code{fracRootResistance (=0.40)}: Fraction of plant total resistance (leaf+stem+root) that corresponds to root system.}
   \item{\code{averageFracRhizosphereResistance (=0.15)}: Fraction to total continuum (leaf+stem+root+rhizosphere) resistance that corresponds to rhizosphere (averaged across soil water potential values).}
   \item{\code{multiLayerBalance (=FALSE)}: Flag to indicate multiple canopy energy balance. If \code{FALSE}, canopy is considered a single layer for energy balance.}
}
\bold{Growth}:
\itemize{
   \item{\code{allowDessication (=TRUE)}: Boolean flag to indicate that mortality by dessication is allowed.}
   \item{\code{allowStarvation (=TRUE)}: Boolean flag to indicate that mortality by starvation is allowed.}
   \item{\code{allocationStrategy (="Plant_kmax")}: Strategy for allocation (either "Plant_kmax", for constant maximum plant conductance, or "Al2As" for constant Huber value).}
   \item{\code{nonStomatalPhotosynthesisLimitation (=TRUE)}: A flag to apply limitations of pho-tosynthesis when leaf sugar levels become too high.}
   \item{\code{phloemConductanceFactor (= 0.2))}: Factor to transform stem xylem conductance to stem phloem conductance (only for transpirationMode = "Sperry").}
   \item{\code{nonSugarConcentration (= 0.3)}: Non-sugar (inorganic) solute concentration  (mol·l-1) in cells.}
   \item{\code{equilibriumOsmoticConcentration [=  c(leaf = 0.8, sapwood = 0.6)]}: Equilibrium osmotic concentrations (mol·l-1) for leaf and sapwood cells. The difference helps phloem transport. }
   \item{\code{minimumSugarForGrowth [=   c(leaf = 0.25, sapwood=0.25, fineroot = 0.25)]}: Minimum concentration of sugars (mol·l-1) for growth to occur in different tissues.}
   \item{\code{respirationRates [=   c(leaf = 0.00260274, sapwood = 6.849315e-05, fineroot = 0.002054795)]}: Maintenance respiration rates for different tissues (g gluc · g dry-1 · day-1). Values for sapwood and fine roots are used for all species, whereas the value for leaves is the default used when leaf respiration is not specified via SpParams (\code{RERleaf}).}
   \item{\code{turnoverRates [=   c(sapwood = 0.0001261398, fineroot = 0.001897231)]}: Turnover rates (day-1) for sapwood and fineroots. Default are equivalent to 4.5\% and 50\% annual turnover for sapwood and fine roots, respectively.}
   \item{\code{constructionCosts [=   c(leaf = 1.5, sapwood = 1.47, fineroot = 1.30)]}: Construction costs, including respiration and structural carbon, per dry weight of new tissue (g gluc · g dry -1).}
   \item{\code{maximumRelativeGrowthRates [=   c(leaf = 0.01, sapwood = 0.002, fineroot = 0.1)]}: Maximum relative growth rates for leaves (m2 leaf ·cm-2 sapwood· day-1), sapwood (cm2 sapwood ·cm-2 sapwood· day-1) and fine roots (g dw · g dw -1 · day -1). Values for leaves and fine roots are used for all species, whereas the value for sapwood is the default used when relative growth rate is not specified via SpParams (\code{RGRsapwoodmax}).}
}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{spwbInput}}, \code{\link{spwb}}
}
