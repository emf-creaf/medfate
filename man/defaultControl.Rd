\encoding{UTF-8}
\name{defaultControl}
\alias{defaultControl}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default control parameters for models}
\description{
Creates a list with global default parameters for simulation models.
}
\usage{
defaultControl() 
}
\details{
The function returns a list with default parameters. Users can change those defaults that need to be set to other values and use the list as input for model functions. The relevant parameters are different for each model function.
}
\value{
A list, with the following options:
\itemize{
\item{\code{verbose (=TRUE)}: Boolean flag to indicate console output during calculations.}
\item{\code{subdailyResults (=FALSE)}: Boolean flag to force subdaily results to be stored (as a list called 'subdaily' of \code{\link{spwb_day}} objects, one by simulated date) in calls to \code{\link{spwb}}.}
\item{\code{soilFunctions (="SX")}: Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).}
\item{\code{snowpack (=TRUE)}: Boolean flag to indicate the simulation of snow accumulation and melting.}
\item{\code{drainage (=TRUE)}: Boolean flag to indicate the simulation of deep drainage.}
\item{\code{leafPhenology (=TRUE)}: Boolean flag to indicate the simulation of leaf phenology for winter-deciduous species.}
\item{\code{unlimitedSoilWater (=FALSE)}: Boolean flag to indicate the simulation of plant transpiration assuming that soil water is always at field capacity.}
\item{\code{plantWaterPools (=FALSE)}: Boolean flag to indicate the simulation of water balance assuming that each cohort has its own water pools (mixing between pools of different cohorts is regulated via \code{poolOverlapFactor}).}
\item{\code{poolOverlapFactor (= 0.2)}: Maximum daily rate of mixing between water pools of different cohorts (actual rate depends on the proportion of fine roots, so the actual mixing rate is normally small for deep layers).}
\item{\code{transpirationMode (="Granier")}: Transpiration model (either 'Granier' or 'Sperry'). See \code{\link{spwbInput}}.}
\item{\code{cuticularTranspiration (=TRUE)}: Boolean flag to indicate the simulation of cuticular transpiration (at night or when stomata are fully closed) according to minimum stomatal conductance.}
\item{\code{gainModifier, costModifier (= 1)}: Modifiers (exponents) of the gain and cost functions defined in Sperry et al. (2016).}
\item{\code{ndailysteps (= 24)}: Number of steps into which each day is divided for determination of stomatal conductance, transpiration and photosynthesis (24 equals 1-hour intervals).}
\item{\code{verticalLayerSize (= 100)}: The size of vertical layers (in cm) for photosynthesis calculation.}
%\item{\code{nStemSegments (= 1)}: Number of segments within the stem.}
\item{\code{capacitance (=FALSE)}: Whether the effect of plant water compartments is considered in simulations.}
\item{\code{cavitationRefill (= "total")}: A string indicating how refilling of embolized conduits is done ("none" - no refilling, "annual" - every first day of the year, "rate" - following a rate of new sapwood formation, "total" - instantaneous complete refilling).}
\item{\code{refillMaximumRate (= 0.05)}: Maximum rate of daily refilling of embolized conduits as sapwood area per leaf area (in cm2·m-2·day-1).}
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
\item{\code{defaultWindSpeed (= 2.5)}: Default wind speed value (in m/s) to be used when missing from data. }
\item{\code{Catm (=386)}: Atmospheric CO2 concentration (in micromol·mol-1 = ppm).}
\item{\code{fracLeafResistance (= NA)}: Fraction of plant total resistance (leaf+stem+root) that corresponds to leaves. This fraction is used if \code{VCleaf_kmax = NA}.}
\item{\code{fracRootResistance (=0.40)}: Fraction of plant total resistance (leaf+stem+root) that corresponds to root system. If \code{fracRootResistance = NA} then root conductance is calculated from root xylem conductivity and root distribution.}
\item{\code{averageFracRhizosphereResistance (=0.15)}: Fraction to total continuum (leaf+stem+root+rhizosphere) resistance that corresponds to rhizosphere (averaged across soil water potential values).}
\item{\code{storagePool (="none")}: Whether carbon storage pools are considered (either "none", "one" or "two").}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{
 \code{\link{spwbInput}}, \code{\link{spwb}}
}