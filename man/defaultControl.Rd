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
   \item{\code{verbose (=TRUE)}: Boolean flag to indicate console output during calculations. In function \code{fordyn} \code{verbose} is always set to FALSE.}
   \item{\code{modifyInput (=TRUE)}: Boolean flag to indicate that simulations will modify input object. If set to FALSE, simulations will not modify the input R object but return the current (modified) state variables within the output. In function \code{fordyn} \code{modifyInput} is always set to FALSE.}
   \item{\code{subdailyResults (=FALSE)}: Boolean flag to force subdaily results to be stored (as a list called 'subdaily' of \code{\link{spwb_day}} objects, one by simulated date) in calls to \code{\link{spwb}}. In function \code{fordyn} \code{subdailyResults} is always set to FALSE.}
}
\bold{Water balance}:
\itemize{
   \item{\code{transpirationMode (="Granier")}: Transpiration model (either 'Granier' or 'Sperry'). See \code{\link{spwbInput}}.}
   \item{\code{soilFunctions (="SX")}: Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten). If \code{transpirationMode = "Sperry"} then soilFunctions is set by default to \code{'VG'}.}
   \item{\code{defaultWindSpeed (= 2.5)}: Default wind speed value (in m/s) to be used when missing from data. }
   \item{\code{snowpack (=TRUE)}: Boolean flag to indicate the simulation of snow accumulation and melting.}
   \item{\code{leafPhenology (=TRUE)}: Boolean flag to indicate the simulation of leaf phenology for winter-deciduous species.}
  \item{\code{rockyLayerDrainage (=TRUE)}: Boolean flag to indicate the simulation of drainage from rocky layers (> 95\% of rocks).}
  \item{\code{unlimitedSoilWater (=FALSE)}: Boolean flag to indicate the simulation of plant transpiration assuming that soil water is always at field capacity.}
  \item{\code{unfoldingDD (=300)}: Degree-days for complete leaf unfolding after budburst has occurred.}
  \item{\code{plantWaterPools (=FALSE)}: Boolean flag to indicate the simulation of water balance assuming that each cohort has its own water pools.}
  \item{\code{verticalLayerSize (= 100)}: Size of vertical layers (in cm) for the calculation of light extinction (and photosynthesis).}
  \item{\code{windMeasurementHeight (=200)}: Height (in cm) over the canopy corresponding to wind measurements.}
  \item{\code{cavitationRefill (= "total")}: A string indicating how refilling of embolized conduits is done:
    \itemize{
      \item{"none" - no refilling.}
      \item{"annual" - every first day of the year.}
      \item{"rate" - following a rate of new sapwood formation (only available if \code{transpirationMode = "Sperry"})}
      \item{"total" - instantaneous complete refilling.}
    }
  }
}
\bold{Water balance with 'Sperry'}:
\itemize{
   \item{\code{ndailysteps (= 24)}: Number of steps into which each day is divided for determination of stomatal conductance, transpiration and photosynthesis (24 equals 1-hour intervals).}
   \item{\code{nsubsteps (= 3600)}: Number of substeps into which each step is divided for multi-layer canopy energy balance solving.}
   \item{\code{capacitance (=FALSE)}: Whether the effect of plant water compartments is considered in simulations.}
   \item{\code{multiLayerBalance (=FALSE)}: Flag to indicate multiple canopy energy balance. If \code{FALSE}, canopy is considered a single layer for energy balance.}
   \item{\code{cochard (=FALSE)}: Boolean flag to indicate the stomatal closure when leaf water potential is below turgor loss point so that transpiration depends on cuticular properties (Martin-StPaul et al. 2017).}
   \item{\code{taper (= TRUE)}: Whether taper of xylem conduits is accounted for when calculating aboveground stem conductance from xylem conductivity.}
   \item{\code{gainModifier, costModifier (= 1)}: Modifiers (exponents) of the gain and cost functions defined in Sperry et al. (2016).}
   \item{\code{klatstem (= 0.01)}: Stem symplastic-apoplastic lateral conductance (in mmol·s-1·m-2·MPa-1). Only relevant when \code{capacitance = TRUE}.}
   \item{\code{klatleaf (= 0.01)}: Leaf symplastic-apoplastic lateral conductance (in mmol·s-1·m-2·MPa-1). Only relevant when \code{capacitance = TRUE}.}
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
   \item{\code{boundaryLayerSize (= 2000)}: Size of the boundary layer (in cm) over the canopy (relevant for multi-layer canopy energy balance).}
   \item{\code{refillMaximumRate (= 0.05)}: Maximum rate of daily refilling of embolized conduits as sapwood area per leaf area (in cm2·m-2·day-1).}
}
\bold{Growth/mortality}:
\itemize{
   \item{\code{allowDessication (=TRUE)}: Boolean flag to indicate that mortality by dessication is allowed.}
   \item{\code{allowStarvation (=TRUE)}: Boolean flag to indicate that mortality by starvation is allowed.}
   \item{\code{allowDefoliation (=TRUE)}: Boolean flag to indicate that complete drought-driven defoliation is allowed.}
   \item{\code{sinkLimitation (=TRUE): Boolean flag to indicate that temperature and turgor limitations to growth are applied.}}
   \item{\code{ shrubDynamics [= FALSE]}: Boolean flag to allow the application of demographic processes to shrubs.}
   \item{\code{allocationStrategy (="Plant_kmax")}: Strategy for allocation (either "Plant_kmax", for constant maximum plant conductance, or "Al2As" for constant Huber value).}
   \item{\code{nonStomatalPhotosynthesisLimitation (=TRUE)}: A flag to apply limitations of pho-tosynthesis when leaf sugar levels become too high.}
   \item{\code{phloemConductanceFactor (= 0.2))}: Factor to transform stem xylem conductance to stem phloem conductance (only for transpirationMode = "Sperry").}
   \item{\code{nonSugarConcentration (= 0.25)}: Non-sugar (inorganic) solute concentration  (mol·l-1) in cells.}
   \item{\code{equilibriumOsmoticConcentration [=  c(leaf = 0.8, sapwood = 0.6)]}: Equilibrium osmotic concentrations (mol·l-1) for leaf and sapwood cells. The difference between leaf and sapwood values helps maintaining phloem transport. The equilibrium sugar concentration is \code{equilibriumOsmoticConcentration - nonSugarConcentration} defaults to \code{[= c(leaf = 0.55, sapwood = 0.35)]}. }
   \item{\code{minimumRelativeSugarForGrowth [= 0.5]}: Minimum concentration of metabolic sugar relative to equilibrium sugar concentration for growth to occur in different tissues. This value (by default 50\%) should always be larger than \code{mortalityRelativeSugarThreshold [= 0.3]}, so that growth stops before starvation, and smaller than 1 (i.e. 100\%), so that growth occurs under normal metabolic status.}
   \item{\code{respirationRates [=   c(leaf = 0.00260274, sapwood = 6.849315e-05, fineroot = 0.002054795)]}: Maintenance respiration rates for different tissues (g gluc · g dry-1 · day-1). Sapwood respiration rates apply to parenchymatic tissue only. Values for sapwood and fine roots are used for all species, whereas the value for leaves is the default used when leaf respiration is not specified via SpParams (\code{RERleaf}).}
   \item{\code{turnoverRates [=   c(sapwood = 0.0001261398, fineroot = 0.001897231)]}: Turnover rates (day-1) for sapwood and fineroots. Default are equivalent to 4.5\% and 50\% annual turnover for sapwood and fine roots, respectively.}
   \item{\code{constructionCosts [=   c(leaf = 1.5, sapwood = 1.47, fineroot = 1.30)]}: Construction costs, including respiration and structural carbon, per dry weight of new tissue (g gluc · g dry -1).}
   \item{\code{maximumRelativeGrowthRates [=   c(leaf = 0.01, sapwood = 0.002, fineroot = 0.1)]}: Maximum relative growth rates for leaves (m2 leaf ·cm-2 sapwood· day-1), sapwood (cm2 sapwood ·cm-2 sapwood· day-1) and fine roots (g dw · g dw -1 · day -1). Values for leaves and fine roots are used for all species, whereas the value for sapwood is the default used when relative growth rate is not specified via SpParams (\code{RGRsapwoodmax}).}
   \item{\code{mortalityMode [= "density/deterministic"]}: String describing how mortality is applied. Current accepted values are combinations of "cohort" vs "density" (for whole-cohort mortality vs reduction of stem density) and "deterministic" vs. "stochastic".}
   \item{\code{mortalityBaselineRate [= 0.01]}: Deterministic proportion or probability specifying the baseline reduction of cohort's density occurring in a year (for \code{mortalityMode = "density/deterministic" or "density/stochastic").}}
   \item{\code{mortalityRelativeSugarThreshold [= 0.3]}: Threshold of stem sugar concentration relative to equilibrium values (by default 30\%), resulting in an increased mortality rate/probability whenever levels are below.}
   \item{\code{mortalityRWCThreshold [= 0.3]}: Threshold of stem relative water content resulting in increased mortality rate/probability whenever levels are below.}
}
\bold{Forest dynamics}:
\itemize{
   \item{\code{ recruitmentMode [= "deterministic"]}: String describing how recruitment is applied. Current accepted values are "deterministic" or "stochastic".}
   \item{\code{ removeDeadCohorts [= TRUE]}: Boolean flag to indicate the removal of cohorts whose density is too low.}
   \item{\code{ minimumCohortDensity [= 1]}: Threshold of density resulting in cohort removal.}
   \item{\code{ seedRain [= NULL]}: Vector of species codes whose seed rain is to be simulated. If \code{NULL} the species identity of seed rain is taken from species currently present in the forest stand and with minimum size (see below).}
   \item{\code{ seedProductionTreeHeight [= 300]}: Default minimum tree height for producing seeds (when species parameter \code{SeedProductionHeight} is missing).}
   \item{\code{ seedProductionShrubHeight [= 30]}: Default minimum shrub height for producing seeds (when species parameter \code{SeedProductionHeight} is missing).}
   \item{\code{ minTempRecr [= 0]}: Default threshold of minimum average temperature of the coldest month necessary for recruiting (when species parameter \code{MinTempRecr} is missing).} 
   \item{\code{ minMoistureRecr [= 0.3]}: Default threshold of minimum moisture index (annual precipitation over annual ETP) necessary for recruiting (when species parameter \code{MinMoistureRecr} is missing).} 
   \item{\code{ minFPARRecr [= 10]}: Default threshold of minimum fraction of PAR (in \%)  reaching the ground necessary for recruiting (when species parameter \code{MinFPARRecr} is missing).} 
   \item{\code{ recrTreeDBH [= 1]}: Default DBH (cm) for recruited trees  (when species parameter \code{RecrTreeDBH} is missing).}
   \item{\code{ recrTreeDensity [= 100]}: Default density (ind·ha-1) for recruited trees  (when species parameter \code{RecrTreeDensity} is missing).}
   \item{\code{ recrTreeHeight [= 100]}: Default height (cm) for recruited trees  (when species parameter \code{RecrTreeHeight} is missing).}
   \item{\code{ recrShrubCover [= 1]}: Default cover (\%) for recruited shrubs  (when species parameter \code{RecrShrubCover} is missing).}
   \item{\code{ recrShrubHeight [= 100]}: Default height (cm) for recruited shrubs  (when species parameter \code{RecrShrubHeight} is missing).}
   \item{\code{ recrTreeZ50 [= 100]}: Default value for Z50 (mm) in recruited trees  (when species parameter \code{RecrZ50} is missing).}
   \item{\code{ recrShrubZ50 [= 50]}: Default value for Z50 (mm) in recruited shrubs  (when species parameter \code{RecrZ50} is missing).}
   \item{\code{ recrTreeZ95 [= 1000]}: Default value for Z95 (mm) in recruited trees  (when species parameter \code{RecrZ50} is missing).}
   \item{\code{ recrShrubZ50 [= 500]}: Default value for Z95 (mm) in recruited shrubs  (when species parameter \code{RecrZ50} is missing).}
}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\seealso{
 \code{\link{spwbInput}}, \code{\link{spwb}}
}
