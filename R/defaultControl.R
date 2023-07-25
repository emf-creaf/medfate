#' Control parameters for simulation models
#'
#' Creates a list control parameters default values for simulations
#' 
#' @param transpirationMode Transpiration model (either 'Granier', 'Sperry' or 'Cochard'). See \code{\link{spwbInput}}.
#' 
#' @details The function returns a list with default parameters. 
#' Users can change those defaults that need to be set to other values and use the list as input for model functions. 
#' The relevant parameters are different for each model function.
#' 
#' @return 
#' A list, with the following options (default values in brackets):
#' \itemize{
#'   \bold{General}:
#'     \itemize{
#'       \item{\code{verbose [= TRUE]}: Boolean flag to indicate console output during calculations. In function \code{fordyn} \code{verbose} is always set to FALSE.}
#'       \item{\code{fillMissingRootParams [= TRUE]}: Boolean flag to indicate that initializing functions should provide estimates for Z50 and Z95 if these are lacking in the forest data. Note that if \code{fillMissingRootParams} is set to \code{FALSE} then simulations may fail if the user does not provide values for Z50 and Z95 in tree or shrub data.}
#'       \item{\code{fillMissingSpParams [= TRUE]}: Boolean flag to indicate that initializing functions should provide estimates for functional parameters if these are lacking in the species parameter table \code{\link{SpParams}}. Note that if \code{fillMissingSpParams} is set to \code{FALSE} then simulations may fail if the user does not provide values for required parameters.}
#'       \item{\code{standResults [= TRUE]}: Boolean flag to keep stand-level results (in a list called 'Stand').}
#'       \item{\code{soilResults [= TRUE]}: Boolean flag to keep stand-level results (in a list called 'Soil').}
#'       \item{\code{plantResults [= TRUE]}: Boolean flag to keep plant-level results (in a list called 'Plants').}
#'       \item{\code{leafResults [= TRUE]}: Boolean flag to keep leaf-level results (in elements called 'SunlitLeaves' and 'ShadeLeaves').}
#'       \item{\code{temperatureResults [= TRUE]}: Boolean flag to keep temperature results (in elements called 'Temperature' and 'TemperatureLayers').}
#'       \item{\code{subdailyResults [= FALSE]}: Boolean flag to force subdaily results to be stored (as a list called 'subdaily' of \code{\link{spwb_day}} objects, one by simulated date) in calls to \code{\link{spwb}}. In function \code{fordyn} \code{subdailyResults} is always set to FALSE.}
#'       \item{\code{fireHazardResults [= FALSE]}: Boolean flag to force calculation of daily fire hazard.}
#'       \item{\code{fireHazardStandardWind [= NA]}: Wind speed (in m/s) for fire-hazard estimation. If missing, actual wind-speed is used.}
#'       \item{\code{fireHazardStandardDFMC [= NA]}: Dead fuel moisture content for fire-hazard estimation. If missing, estimation from current weather is used.}
#'     }
#'   \bold{Water balance} (functions \code{\link{spwb}}, \code{\link{pwb}} or \code{\link{spwb_day}}):
#'     \itemize{
#'       \item{\code{transpirationMode [= "Granier"]}: Transpiration model (either 'Granier', 'Sperry' or 'Cochard'). See \code{\link{spwbInput}}.}
#'       \item{\code{soilFunctions [= "VG"]}: Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten). 
#'                  If \code{transpirationMode} is 'Sperry' or 'Cochard' then soilFunctions is forced to \code{'VG'}.
#'                  Only simulations with 'Granier' are allowed to use Saxton functions.}
#'       \item{\code{defaultWindSpeed [= 2.5]}: Default wind speed value (in m/s) to be used when missing from data. }
#'       \item{\code{defaultCO2 [= 386]}: Default atmospheric (abovecanopy) CO2 concentration (in micromol·mol-1 = ppm). This value will be used whenever CO2 concentration is not specified in the weather input. }
#'       \item{\code{snowpack [= TRUE]}: Boolean flag to indicate the simulation of snow accumulation and melting.}
#'       \item{\code{leafPhenology [= TRUE]}: Boolean flag to indicate the simulation of leaf phenology for winter-deciduous species.}
#'       \item{\code{rockyLayerDrainage [= TRUE]}: Boolean flag to indicate the simulation of drainage from rocky layers (> 95\% of rocks).}
#'       \item{\code{unlimitedSoilWater [= FALSE]}: Boolean flag to indicate the simulation of plant transpiration assuming that soil water is always at field capacity.}
#'       \item{\code{unfoldingDD [= 300]}: Degree-days for complete leaf unfolding after budburst has occurred.}
#'       \item{\code{rhizosphereOverlap [= "total"]}: A string indicating the degree of rhizosphere spatial overlap between plant cohorts:
#'           \itemize{
#'             \item{"none" - no overlap (independent water pools).}
#'             \item{"partial" - partial overlap determined by coarse root volume.}
#'             \item{"total" - total overlap (plants extract from common soil pools).}
#'           }
#'       }
#'       \item{\code{verticalLayerSize [= 100]}: Size of vertical layers (in cm) for the calculation of light extinction (and photosynthesis).}
#'       \item{\code{windMeasurementHeight [= 200]}: Height (in cm) over the canopy corresponding to wind measurements.}
#'       \item{\code{cavitationRefill [= "total"]}: A string indicating how refilling of embolized conduits is done:
#'           \itemize{
#'             \item{"none" - no refilling.}
#'             \item{"annual" - every first day of the year.}
#'             \item{"rate" - following a rate of new sapwood formation.}
#'             \item{"total" - instantaneous complete refilling.}
#'           }
#'       }
#'     }
#'   \bold{Water balance} (functions \code{\link{spwb}}, \code{\link{pwb}} or \code{\link{spwb_day}} when \code{traspirationMode = "Granier"}):
#'     \itemize{
#'       \item{\code{hydraulicRedistributionFraction [= 0.1]}: Fraction of plant transpiration corresponding to hydraulic redistribution.}
#'     }
#'   
#'   \bold{Water balance} (functions \code{\link{spwb}}, \code{\link{pwb}} or \code{\link{spwb_day}} when \code{traspirationMode = "Sperry"} or \code{traspirationMode = "Cochard"}):
#'     \itemize{
#'       \item{\code{ndailysteps [= 24]}: Number of steps into which each day is divided for determination of stomatal conductance, transpiration and photosynthesis (24 equals 1-hour intervals).}
#'       \item{\code{nsubsteps [= 3600]}: Number of substeps into which each step is divided for multi-layer canopy energy balance solving.}
#'       \item{\code{capacitance [= FALSE]}: Whether the effect of plant water compartments is considered in simulations.}
#'       \item{\code{multiLayerBalance [= FALSE]}: Flag to indicate multiple canopy energy balance. If \code{FALSE}, canopy is considered a single layer for energy balance.}
#'       \item{\code{taper [= TRUE]}: Whether taper of xylem conduits is accounted for when calculating aboveground stem conductance from xylem conductivity.}
#'       \item{\code{maximumStemConductance [= 10]}: Maximum allowed value for the stem maximum hydraulic conductance (in mmol·s-1·m-2·MPa-1). Introduced to avoid excessive hydraulic redistribution caused by species with small size (i.e. very large stem conductance).}
#'       \item{\code{klatstem [= 0.01]}: Stem symplastic-apoplastic lateral conductance (in mmol·s-1·m-2·MPa-1). Only relevant when \code{capacitance = TRUE}.}
#'       \item{\code{klatleaf [= 0.01]}: Leaf symplastic-apoplastic lateral conductance (in mmol·s-1·m-2·MPa-1). Only relevant when \code{capacitance = TRUE}.}
#'       \item{\code{numericParams}: A list with the following elements:
#'           \itemize{
#'             \item{\code{maxNsteps [= 400]}: Maximum number of steps in supply function.}
#'             \item{\code{ntrial [= 200]}: Number of iteration trials when finding root of equation system.}
#'             \item{\code{psiTol [= 0.0001]}: Tolerance value for water potential.}
#'             \item{\code{ETol [= 0.0001]}: Tolerance value for flow.}
#'           }
#'       }
#'       \item{\code{thermalCapacityLAI [= 1000000]}: Thermal canopy capacitance per LAI unit.}
#'       \item{\code{fracLeafResistance [= NA]}: Fraction of plant total resistance (leaf+stem+root) that corresponds to leaves. This fraction is used if \code{VCleaf_kmax = NA}.}
#'       \item{\code{fracRootResistance [= 0.40]}: Fraction of plant total resistance (leaf+stem+root) that corresponds to root system.}
#'       \item{\code{averageFracRhizosphereResistance [= 0.15]}: Fraction to total continuum (leaf+stem+root+rhizosphere) resistance that corresponds to rhizosphere (averaged across soil water potential values).}
#'       \item{\code{boundaryLayerSize [= 2000]}: Size of the boundary layer (in cm) over the canopy (relevant for multi-layer canopy energy balance).}
#'       \item{\code{refillMaximumRate [= 0.05]}: Maximum rate of daily refilling of embolized conduits as sapwood area per leaf area (in cm2·m-2·day-1).}
#'     }
#'   \bold{Forest growth} (functions \code{\link{growth}} or \code{\link{growth_day}}):
#'     \itemize{
#'       \item{\code{subdailyCarbonBalance [= FALSE]}: Boolean flag to indicate that labile carbon balance should be conducted at sub-daily steps (applies only to transpirationMode = "Sperry").}
#'       \item{\code{allowDessication [= TRUE]}: Boolean flag to indicate that mortality by dessication is allowed.}
#'       \item{\code{allowStarvation [= TRUE]}: Boolean flag to indicate that mortality by starvation is allowed.}
#'       \item{\code{sinkLimitation [= TRUE]}: Boolean flag to indicate that temperature and turgor limitations to growth are applied.}
#'       \item{\code{shrubDynamics [= TRUE]}: Boolean flag to allow the application of demographic processes to shrubs.}
#'       \item{\code{herbDynamics [= TRUE]}: Boolean flag to allow dynamic herb leaf area as a function of shading due to leaf area of woody cohorts.}
#'       \item{\code{allocationStrategy [= "Al2As"]}: Strategy for allocation (either "Plant_kmax", for constant maximum plant conductance, or "Al2As" for constant Huber value).}
#'       \item{\code{phloemConductanceFactor [= 0.2])}: Factor to transform stem xylem conductance to stem phloem conductance (only for transpirationMode = "Sperry").}
#'       \item{\code{nonSugarConcentration [= 0.25]}: Non-sugar (inorganic) solute concentration  (mol·l-1) in cells.}
#'       \item{\code{equilibriumOsmoticConcentration [=  c(leaf = 0.8, sapwood = 0.6)]}: Equilibrium osmotic concentrations (mol·l-1) for leaf and sapwood cells. The difference between leaf and sapwood values helps maintaining phloem transport. The equilibrium sugar concentration is \code{equilibriumOsmoticConcentration - nonSugarConcentration} defaults to \code{[= c(leaf = 0.55, sapwood = 0.35)]}. }
#'       \item{\code{minimumRelativeStarchForGrowth [= 0.50]}: Default minimum concentration of storage carbon (starch), relative to the maximum storage capacity, for sapwood growth to occur, when not specified via SpParams (\code{RSSG}). }
#'       \item{\code{constructionCosts [= c(leaf = 1.5, sapwood = 1.47, fineroot = 1.30)]}: Default construction costs, including respiration and structural carbon, per dry weight of new tissue (g gluc · g dry -1) when not specified via SpParams (\code{CCleaf}, \code{CCsapwood} and \code{CCfineroot}).}
#'       \item{\code{senescenceRates [= c(sapwood = 0.0001261398, fineroot = 0.001897231)]}: Default senescence rates (day-1) for sapwood and fineroots when not specified via SpParams (\code{SRsapwood} and \code{SRfineroot}). Defaults are equivalent to 9\%, 5\% and 50\% annual turnover for gymnosperm sapwood, angiosperm sapwood and fine roots, respectively.}
#'       \item{\code{maximumRelativeGrowthRates [=   c(leaf = 0.03, cambium = 0.005, sapwood = 0.002, fineroot = 0.1)]}: Default maximum relative growth rates for leaves (m2 leaf ·cm-2 sapwood· day-1), tree sapwood (cm2 sapwood· cm-1 cambium · day-1), shrub sapwood (cm2 sapwood ·cm-2 sapwood· day-1) and fine roots (g dw · g dw -1 · day -1) when not specified via SpParams (\code{RGRleafmax}, \code{RGRcambiummax} , \code{RGRsapwoodmax} and \code{RGRfinerootmax}, respectively).}
#'       \item{\code{mortalityMode [= "density/deterministic"]}: String describing how mortality is applied. Current accepted values are combinations of "cohort" vs "density" (for whole-cohort mortality vs reduction of stem density) and "deterministic" vs. "stochastic".}
#'       \item{\code{mortalityBaselineRate [= 0.0015]}: Default deterministic proportion or probability specifying the baseline reduction of cohort's density occurring in a year (for \code{mortalityMode = "density/deterministic" or "density/stochastic").}}
#'       \item{\code{mortalityRelativeSugarThreshold [= 0.4]}: Threshold of stem sugar concentration relative to the equilibrium sugar concentration, resulting in an increased starvation mortality rate/probability whenever levels are below.}
#'       \item{\code{mortalityRWCThreshold [= 0.4]}: Threshold of stem relative water content resulting in increased mortality rate/probability whenever levels are below.}
#'       \item{\code{recrTreeDBH [= 1]}: Default DBH (cm) for recruited trees  (when species parameter \code{RecrTreeDBH} is missing).}
#'       \item{\code{recrTreeDensity [= 3000]}: Default density (ind·ha-1) for recruited trees  (when species parameter \code{RecrTreeDensity} is missing).}
#'       \item{\code{ingrowthTreeDBH [= 7.5]}: Default DBH (cm) for ingrowth trees  (when species parameter \code{RecrTreeDBH} is missing).}
#'       \item{\code{ingrowthTreeDensity [= 127]}: Default density (ind·ha-1) for ingrowth trees  (when species parameter \code{RecrTreeDensity} is missing).}
#'    }
#'   \bold{Forest dynamics} (function \code{\link{fordyn}}):
#'    \itemize{
#'      \item{\code{allowRecruitment [= TRUE]}: Boolean flag to indicate that recruitment from seeds is allowed.}
#'      \item{\code{recruitmentMode [= "stochastic"]}: String describing how recruitment from seeds is applied. Current accepted values are "deterministic" or "stochastic".}
#'      \item{\code{allowResprouting [= TRUE]}: Boolean flag to indicate that resprouting is allowed.}
#'      \item{\code{removeEmptyCohorts [= TRUE]}: Boolean flag to indicate the removal of cohorts whose density is too low.}
#'      \item{\code{minimumTreeCohortDensity [= 1]}: Threshold of tree density resulting in cohort removal.}
#'      \item{\code{minimumShrubCohortCover [= 0.01]}: Threshold of shrub cover resulting in cohort removal.}
#'      \item{\code{dynamicallyMergeCohorts [= FALSE]}: Boolean flag to indicate that cohorts should be merged when possible. This option speeds up calculations but results in a loss of cohort identity and reinitialization of many state variables.}
#'      \item{\code{seedRain [= NULL]}: Vector of species codes whose seed rain is to be simulated. If \code{NULL} the species identity of seed rain is taken from species currently present in the forest stand and with minimum size (see below).}
#'      \item{\code{seedProductionTreeHeight [= 300]}: Default minimum tree height for producing seeds (when species parameter \code{SeedProductionHeight} is missing).}
#'      \item{\code{seedProductionShrubHeight [= 30]}: Default minimum shrub height for producing seeds (when species parameter \code{SeedProductionHeight} is missing).}
#'      \item{\code{probRecr [= 0.05]}: Default annual probability of seed-recruitment (when species parameter \code{ProbRecr} is missing).}
#'      \item{\code{minTempRecr [= 0]}: Default threshold of minimum average temperature of the coldest month necessary for recruiting from seeds (when species parameter \code{MinTempRecr} is missing).} 
#'      \item{\code{minMoistureRecr [= 0.3]}: Default threshold of minimum moisture index (annual precipitation over annual ETP) necessary for seed-recruiting (when species parameter \code{MinMoistureRecr} is missing).} 
#'      \item{\code{minFPARRecr [= 10]}: Default threshold of minimum fraction of PAR (in \%)  reaching the ground necessary for recruiting (when species parameter \code{MinFPARRecr} is missing).} 
#'      \item{\code{recrTreeHeight [= 620]}: Default height (cm) for recruited trees  (when species parameter \code{RecrTreeHeight} is missing).}
#'      \item{\code{recrShrubCover [= 1]}: Default cover (\%) for shrubs recruited from seed  (when species parameter \code{RecrShrubCover} is missing).}
#'      \item{\code{recrShrubHeight [= 25]}: Default height (cm) for recruited shrubs  (when species parameter \code{RecrShrubHeight} is missing).}
#'      \item{\code{recrTreeZ50 [= 100]}: Default value for Z50 (mm) in seed-recruited trees  (when species parameter \code{RecrZ50} is missing).}
#'      \item{\code{recrShrubZ50 [= 50]}: Default value for Z50 (mm) in seed-recruited shrubs  (when species parameter \code{RecrZ50} is missing).}
#'      \item{\code{recrTreeZ95 [= 1000]}: Default value for Z95 (mm) in seed-recruited trees  (when species parameter \code{RecrZ50} is missing).}
#'      \item{\code{recrShrubZ50 [= 500]}: Default value for Z95 (mm) in seed-recruited shrubs  (when species parameter \code{RecrZ50} is missing).}
#'   }
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{spwbInput}}, \code{\link{spwb}}, \code{\link{growth}}, \code{\link{fordyn}}
#' 
#' @name defaultControl
defaultControl<-function(transpirationMode = "Granier") {
  transpirationMode <- match.arg(transpirationMode, c("Granier", "Sperry", "Cochard"))
  return(list(
    #For all functions
    fillMissingRootParams = TRUE,
    fillMissingSpParams = TRUE,
    verbose = TRUE,
    subdailyResults = FALSE,
    standResults = TRUE,
    soilResults = TRUE,
    plantResults = TRUE,
    leafResults = TRUE,
    temperatureResults = TRUE,
    fireHazardResults = FALSE,
    fireHazardStandardWind = NA,
    fireHazardStandardDFMC = NA,
    
    # For water balance
    transpirationMode = transpirationMode,
    soilFunctions = "VG",
    defaultWindSpeed = 2.5, #m/s
    defaultCO2 = 386, #ppm
    snowpack = TRUE,
    leafPhenology = TRUE,
    rockyLayerDrainage = TRUE,
    unlimitedSoilWater = FALSE,
    rhizosphereOverlap = "total",
    unfoldingDD = 300,
    verticalLayerSize = 100,
    windMeasurementHeight = 200,
    cavitationRefill = "total",
    
    #spwb with granier
    hydraulicRedistributionFraction = 0.1,
    
    #spwb with sperry
    ndailysteps = 24,
    nsubsteps = 3600,
    capacitance = FALSE,
    taper = TRUE,
    multiLayerBalance = FALSE,
    maximumStemConductance = 10,
    klatstem = 0.01, # stem symplastic-apoplastic lateral conductance
    klatleaf = 0.01, # leaf symplastic-apoplastic lateral conductance
    numericParams=list(maxNsteps = 400, ntrial = 200, psiTol = 0.0001, ETol = 0.0000001),
    fracLeafResistance = NA,
    fracRootResistance = 0.4,
    averageFracRhizosphereResistance = 0.15,
    thermalCapacityLAI = 1000000,
    boundaryLayerSize = 2000,
    refillMaximumRate = 0.05,
    
    # growth/mortality
    subdailyCarbonBalance = FALSE,
    allowDessication = TRUE,
    allowStarvation = TRUE,
    sinkLimitation = TRUE,
    shrubDynamics = TRUE,
    herbDynamics = TRUE,
    allocationStrategy = "Al2As",
    phloemConductanceFactor = 0.2, # phloem conductance per leaf area basis (l*m-2*MPa-1*s-1)
    nonSugarConcentration = 0.25, # mol · l-1
    equilibriumOsmoticConcentration = list(leaf = 0.8, sapwood = 0.6),  # (Paljakka et al. 2017)
    minimumRelativeStarchForGrowth = 0.50,
    constructionCosts = list(leaf = 1.5, 
                             sapwood = 1.47, 
                             fineroot = 1.30), #  g gluc · g dw -1
    senescenceRates = list(sapwood = 0.000135, # day-1 Equivalent to annual 4.8% 1-(1-0.048)^(1.0/365)
                      fineroot = 0.001897231), #day-1 Equivalent to annual 50% 1-(1-0.5)^(1.0/365)
    maximumRelativeGrowthRates = list(leaf = 0.03, # m2 leaf ·cm-2 sapwood· day-1
                                      cambium = 0.0025, # cm2 sapwood ·cm-1 cambium· day-1
                                      sapwood = 0.002, # cm2 sapwood ·cm-2 sapwood· day-1
                                      fineroot = 0.1), # g dw · g dw -1 · day -1
    mortalityMode = "density/deterministic",
    mortalityBaselineRate = 0.0015,
    mortalityRelativeSugarThreshold = 0.4,
    mortalityRWCThreshold = 0.4,
    recrTreeDBH = 1,
    recrTreeDensity = 3000,
    ingrowthTreeDBH = 7.5,
    ingrowthTreeDensity = 127,
    
    #dynamics
    allowRecruitment = TRUE,
    recruitmentMode = "stochastic",
    allowResprouting = TRUE,
    removeEmptyCohorts=TRUE,
    minimumTreeCohortDensity = 1,
    minimumShrubCohortCover = 0.01,
    dynamicallyMergeCohorts = FALSE,
    seedRain = NULL,
    seedProductionTreeHeight = 300,
    seedProductionShrubHeight = 30,
    probRecr = 0.05,
    minTempRecr	= 0,
    minMoistureRecr	= 0.3,
    minFPARRecr = 10,
    recrTreeHeight = 620,
    recrShrubCover = 1,
    recrShrubHeight = 25,
    recrTreeZ50 = 100,
    recrShrubZ50 = 50,
    recrTreeZ95 = 1000,
    recrShrubZ95 = 500
    
    
    #For forest dynamics
#     freqZopt = 20,
#     sap2tree=TRUE,
#     initRadius = 5000,
#     setDefaults = TRUE,
#     diamClasses = c(0,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.0, seq(15,200, by=5)),
#     heightClasses = c(0,25,50,75,100,125,150,175, 200,seq(250,6000, by=50)),    
    #For fire
#     sparkLambda = 0, #Average number of sparks/day (0 means no fires)
#     dailySparkingTime = 960.0, #in minutes (960 min = 60 min/h * 16h)
#     sparkingProbabilityVector = NULL, #Relative values for unequal ignition probability
#     brigadeDensity =  0.001, #in number of fire brigades/square km2
#     brigadeDispatchTime = 1.0, #in minutes
#     brigadeFireArrivalTime = 50.0, #in minutes
#     extinctionCapacity = 5000, #in kW/(m*min)
#     fireBrandTimeStep = 1/60, #in minutes
#     numBrandsPerBurningTime = 1, # dimensionless
    #For fire behaviour
   # liveFMCmode = "swb",
   # useModelForLive = FALSE
    # initialMoisture = c(35,35,35,35,35),
    # surfaceToVolumeRatios = c(5600,358,98,6200,8000), #dimensionless
    # heatContent = c(18622, 18622, 18622, 19500, 20000), #kJ/kg
    # deadFuelMoistureExtinction = 30 #Percent
    
    #For dispersal
    # maxDispersal = 500,
    #For management
    # managementUnits = NULL,
    # managementPlan = NULL
))
}