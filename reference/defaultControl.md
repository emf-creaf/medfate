# Control parameters for simulation models

Creates a list control parameters default values for simulations

## Usage

``` r
defaultControl(
  transpirationMode = "Granier",
  soilDomains = "buckets",
  rhizosphereOverlap = "total"
)
```

## Arguments

- transpirationMode:

  String containing transpiration model (either 'Granier', 'Sperry' or
  'Sureau'). See
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- soilDomains:

  String containing soil hydrology model (either 'buckets', 'single' or
  'dual'). See
  [`hydrology_soilWaterBalance`](https://emf-creaf.github.io/medfate/reference/hydrology_soilWaterBalance.md).

- rhizosphereOverlap:

  String indicating the assumption with respect to rhizosphere overlap:

  - `total`: All plants extract water from the same pools.

  - `partial`: Plants partially share water pools, depending on soil
    moisture.

  - `none`: Plants extract water from totally-independent water pools.

## Value

A list, with the following options (default values in brackets):

- **General**:

  - `verbose [= TRUE]`: Boolean flag to indicate console output during
    calculations. In function `fordyn` `verbose` is always set to FALSE.

  - `fillMissingRootParams [= TRUE]`: Boolean flag to indicate that
    initializing functions should provide estimates for Z50 and Z95 if
    these are missing in the forest data. Note that if
    `fillMissingRootParams` is set to `FALSE` then simulations may fail
    if the user does not provide values for Z50 and Z95 in tree or shrub
    data.

  - `fillMissingSpParams [= TRUE]`: Boolean flag to indicate that
    initializing functions should provide estimates for functional
    parameters if these are missing in the species parameter table
    [`SpParams`](https://emf-creaf.github.io/medfate/reference/SpParams.md).
    Note that if `fillMissingSpParams` is set to `FALSE` then
    simulations may fail if the user does not provide values for
    required parameters.

  - `fillMissingWithGenusParams [=TRUE]`: Boolean flag to indicate that
    initializing functions should provide estimates from genus value, if
    species-level values are missing in the species parameter table
    [`SpParams`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
    but genus-level ones are not.

  - `standResults [= TRUE]`: Boolean flag to keep stand-level results
    (in a data frame called 'Stand').

  - `soilResults [= TRUE]`: Boolean flag to keep soil-level results (in
    a list called 'Soil').

  - `soilPoolResults [= FALSE]`: Boolean flag to keep soil pool-level
    results (in a list called 'SoilPools'), if
    `rhizosphereOverlap!="total"`.

  - `snowResults [= TRUE]`: Boolean flag to keep snow results (in a data
    frame called 'Snow').

  - `plantResults [= TRUE]`: Boolean flag to keep plant-level
    water/energy/photosynthesis results (in a list called 'Plants').

  - `labileCarbonBalanceResults [= TRUE]`: Boolean flag to keep
    plant-level labile carbon balance results (in a list called
    'LabileCarbonBalance').

  - `plantStructureResults [= TRUE]`: Boolean flag to keep plant-level
    structure results (in a list called 'PlantStructure').

  - `growthMortalityResults [= TRUE]`: Boolean flag to keep plant-level
    growth and mortality results (in a list called 'GrowthMortality').

  - `leafResults [= TRUE]`: Boolean flag to keep leaf-level results (in
    elements called 'SunlitLeaves' and 'ShadeLeaves').

  - `temperatureResults [= TRUE]`: Boolean flag to keep temperature
    results (in elements called 'Temperature' and 'TemperatureLayers').

  - `subdailyResults [= FALSE]`: Boolean flag to force subdaily results
    to be stored (as a list called 'subdaily' of
    [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
    objects, one by simulated date) in calls to
    [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md). In
    function `fordyn` `subdailyResults` is always set to FALSE.

  - `fireHazardResults [= FALSE]`: Boolean flag to force calculation of
    daily fire hazard.

  - `fireHazardStandardWind [= NA]`: Wind speed (in m/s) for fire-hazard
    estimation. If missing, actual wind-speed is used.

  - `fireHazardStandardDFMC [= NA]`: Dead fuel moisture content for
    fire-hazard estimation. If missing, estimation from current weather
    is used.

  **Water balance** (functions
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
  [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)):

  - `transpirationMode [= "Granier"]`: Transpiration model (either
    'Granier', 'Sperry' or 'Sureau'). See
    [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

  - `soilFunctions [= "VG"]`: Soil water retention curve and
    conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van
    Genuchten). If `transpirationMode` is 'Sperry' or 'Sureau' then
    soilFunctions is forced to `'VG'`. Only simulations with 'Granier'
    are allowed to use Saxton functions.

  - `truncateRootDistribution [= FALSE]`: Boolean flag to indicate that
    fine root distribution is to be truncated (i.e. Z100 is estimated
    from Z95 when not provided).

  - `VG_PTF`: String indicating the pedotransfer functions for van
    Genuchten parameters (either 'Toth' or 'Carsel').

  - `ndailysteps [= 24]`: Number of steps into which each day is divided
    for determination of soil water balance, stomatal conductance,
    transpiration and photosynthesis (24 equals 1-hour intervals).

  - `max_nsubsteps_soil [= 300]`: Maximum number of substeps for soil
    water balance solving.

  - `defaultWindSpeed [= 2.5]`: Default wind speed value (in m/s) to be
    used when missing from data.

  - `defaultCO2 [= 386]`: Default atmospheric (abovecanopy) CO2
    concentration (in micromol·mol-1 = ppm). This value will be used
    whenever CO2 concentration is not specified in the weather input.

  - `defaultRainfallIntensityPerMonth [= c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 5.6, 5.6, 5.6, 5.6, 5.6, 1.5)]`:
    A vector of twelve values indicating the rainfall intensity (mm/h)
    per month. By default synoptic storms (1.5 mm/h) are assumed between
    December and June, and convective storms (5.6 mm/h) are assumed
    between July and November.

  - `leafPhenology [= TRUE]`: Boolean flag to indicate the simulation of
    leaf phenology for winter-deciduous species.

  - `bareSoilEvaporation [= TRUE]`: Boolean flag to indicate the
    simulation of evaporation from bare soil.

  - `unlimitedSoilWater [= FALSE]`: Boolean flag to indicate the
    simulation of plant transpiration assuming that soil water is always
    at field capacity.

  - `unfoldingDD [= 300]`: Degree-days for complete leaf unfolding after
    budburst has occurred.

  - `interceptionMode [= "Gash1995"]`: Infiltration model, either
    "Gash1995" or "Liu2001".

  - `infiltrationMode [= "GreenAmpt1911"]`: Infiltration model, either
    "GreenAmpt1911" or "Boughton1989".

  - `infiltrationCorrection [= 5.0]`: Factor to correct infiltration
    amount in the GreenAmpt1911 model in single-domain simulations.

  - `soilDomains [= "buckets"]`: Either 'buckets' (for multi-bucket
    model), 'single' (for single-domain Richards model) or 'dual' (for
    dual-permeability model). See
    [`hydrology_soilWaterBalance`](https://emf-creaf.github.io/medfate/reference/hydrology_soilWaterBalance.md).

  - `rhizosphereOverlap [= "total"]`: A string indicating the degree of
    rhizosphere spatial overlap between plant cohorts:

    - "none" - no overlap (independent water pools).

    - "partial" - partial overlap determined by coarse root volume.

    - "total" - total overlap (plants extract from common soil pools).

  - `fullRhizosphereOverlapConductivity [= 0.01]`: The minimum soil
    hydraulic conductivity (in cm/day) allowing a full connectivity of
    water pools, when `rhizosphereOverlap = "partial"`.

  - `verticalLayerSize [= 100]`: Size of vertical layers (in cm) for the
    calculation of light extinction (and photosynthesis).

  - `windMeasurementHeight [= 200]`: Height (in cm) over the canopy
    corresponding to wind measurements.

  - `segmentedXylemVulnerability [= TRUE/FALSE]`: If `FALSE` leaf and
    root vulnerability curves will be equal to those of stem. By
    default, `segmentedXylemVulnerability = TRUE` for
    `transpirationMode = "Sperry"` and
    `segmentedXylemVulnerability = FALSE` for
    `transpirationMode = "Sureau"`.

  - `leafCavitationEffects, stemCavitationEffects [= FALSE/TRUE]`: A
    flag indicating whether cavitation effects on conductance of leaves
    and stem are applied. Only relevant for
    `transpirationMode = "Sperry"`.

  - `leafCavitationRecovery, stemCavitationRecovery [= "annual"]`: A
    string indicating how recovery of previous cavitation leaf/stem
    xylem is done (only relevant for functions
    [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) and
    [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)):

    - "none" - no recovery.

    - "annual" - every first day of the year.

    - "rate" - following a rate of new leaf or sapwood formation.

    - "total" - instantaneous complete recovery.

  - `cavitationRecoveryMaximumRate [= 0.05]`: Maximum rate of daily
    refilling of embolized conduits as sapwood area per leaf area (in
    cm2·m-2·day-1).

  - `lfmcComponent [= "fine"]`: Plant component used to estimate LFMC,
    either "leaf" or "fine" (for fine fuel).

  **Water balance** (functions
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
  [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
  when `traspirationMode = "Granier"` only):

  - `hydraulicRedistributionFraction [= 0.1]`: Fraction of plant
    transpiration corresponding to hydraulic redistribution.

  **Water balance** (functions
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
  [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
  when `traspirationMode = "Sperry"` or `traspirationMode = "Sureau"`):

  - `nsubsteps_canopy [= 3600]`: Number of substeps into which each step
    is divided for multi-layer canopy energy balance solving.

  - `multiLayerBalance [= FALSE]`: Flag to indicate multiple canopy
    energy balance. If `FALSE`, canopy is considered a single layer for
    energy balance.

  - `sapFluidityVariation [= TRUE]`: Flag to indicate that temperature
    affects sap fluidity (and indirectly plant conductance).

  - `TPhase_gmin [= 37.5]`: Temperature for transition phase of gmin.

  - `Q10_1_gmin [= 1.2]`: Temperature dependance of gmin when T less
    than or equal to TPhase.

  - `Q10_2_gmin [= 4.8]`: Temperature dependance of gmin when T greater
    than TPhase.

  - `taper [= TRUE]`: Whether taper of xylem conduits is accounted for
    when calculating aboveground stem conductance from xylem
    conductivity.

  - `thermalCapacityLAI [= 1000000]`: Thermal canopy capacitance per LAI
    unit.

  - `rootRadialConductance [= 4]`: Radial conductance in roots
    (mmol·s-1·m-2·MPa-1).

  - `averageFracRhizosphereResistance [= 0.15]`: Fraction to total
    continuum (leaf+stem+root+rhizosphere) resistance that corresponds
    to rhizosphere (averaged across soil water potential values).

  - `boundaryLayerSize [= 2000]`: Size of the boundary layer (in cm)
    over the canopy (relevant for multi-layer canopy energy balance).

  **Water balance** (functions
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
  [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
  when `traspirationMode = "Sperry"` only):

  - `numericParams`: A list with the following elements:

    - `maxNsteps [= 400]`: Maximum number of steps in supply function.

    - `ntrial [= 200]`: Number of iteration trials when finding root of
      equation system.

    - `psiTol [= 0.0001]`: Tolerance value for water potential.

    - `ETol [= 0.0001]`: Tolerance value for flow.

  **Water balance** (functions
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
  [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
  when `traspirationMode = "Sureau"` only):

  - `plantCapacitance [= TRUE]`: Whether the effect of (symplasmic)
    plant water compartments is considered in simulations.

  - `cavitationFlux [= TRUE]`: Whether the effect of water flux
    generated by cavitation of apoplasmic tissues is considered in
    simulations.

  - `leafCuticularTranspiration [= TRUE]`: Whether the effect of leaf
    cuticular transpiration is considered in simulations.

  - `stemCuticularTranspiration [= FALSE]`: Whether the effect of stem
    cuticular transpiration is considered in simulations.

  - `C_SApoInit [= 2.0e-5]`: Maximum capacitance of the stem apoplasm
    (mmol·m-2).

  - `C_LApoInit [= 1.0e-5]`: Maximum capacitance of the leaf apoplasm
    (mmol·m-2).

  - `k_SSym [= 0.26]`: Conductance from stem apoplasm to stem symplasm
    (mmol·s-1·m-2·MPa-1).

  - `fractionLeafSymplasm [= 0.5]`: Fraction of the leaf resistance from
    leaf apoplasm to leaf symplasm (\[0-1\]).

  - `gs_NightFrac [= 0.05]`: Stomatal conductance at night as fraction
    of maximum stomatal conductance (\[0-1\]).

  - `stomatalSubmodel [= "Baldocchi"]`: Stomatal regulation sub-model,
    either "Jarvis" or "Baldocchi".

  - `JarvisPAR [= 0.003]`: Parameter regulating the response of stomatal
    conductance to light (PAR) in the Jarvis model.

  - `fTRBToLeaf [= 0.8]`: Fraction of surface of bark exposed to air per
    leaf area.

  **Forest growth** (functions
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md) or
  [`growth_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)):

  - `subdailyCarbonBalance [= FALSE]`: Boolean flag to indicate that
    labile carbon balance should be conducted at sub-daily steps
    (applies only to transpirationMode = "Sperry").

  - `allowDessication [= TRUE]`: Boolean flag to indicate that mortality
    by dessication is allowed.

  - `allowStarvation [= TRUE]`: Boolean flag to indicate that mortality
    by starvation is allowed.

  - `sinkLimitation [= TRUE]`: Boolean flag to indicate that temperature
    and turgor limitations to growth are applied.

  - `shrubDynamics [= TRUE]`: Boolean flag to allow the application of
    demographic processes to shrubs.

  - `herbDynamics [= TRUE]`: Boolean flag to allow dynamic herb leaf
    area as a function of shading due to leaf area of woody cohorts.

  - `allocationStrategy [= "Al2As"]`: Strategy for allocation (either
    "Plant_kmax", for constant maximum plant conductance, or "Al2As" for
    constant Huber value).

  - `phloemConductanceFactor [= 0.2])`: Factor to transform stem xylem
    conductance to stem phloem conductance (only for transpirationMode =
    "Sperry").

  - `nonSugarConcentration [= 0.25]`: Non-sugar (inorganic) solute
    concentration (mol·l-1) in cells.

  - `equilibriumOsmoticConcentration [= c(leaf = 0.8, sapwood = 0.6)]`:
    Equilibrium osmotic concentrations (mol·l-1) for leaf and sapwood
    cells. The difference between leaf and sapwood values helps
    maintaining phloem transport. The equilibrium sugar concentration is
    `equilibriumOsmoticConcentration - nonSugarConcentration` defaults
    to `[= c(leaf = 0.55, sapwood = 0.35)]`.

  - `minimumRelativeStarchForGrowth [= 0.50]`: Default minimum
    concentration of storage carbon (starch), relative to the maximum
    storage capacity, for sapwood growth to occur, when not specified
    via SpParams (`RSSG`).

  - `constructionCosts [= c(leaf = 1.5, sapwood = 1.47, fineroot = 1.30)]`:
    Default construction costs, including respiration and structural
    carbon, per dry weight of new tissue (g gluc · g dry -1) when not
    specified via SpParams (`CCleaf`, `CCsapwood` and `CCfineroot`).

  - `senescenceRates [= c(sapwood = 0.0001261398, fineroot = 0.001897231)]`:
    Default senescence rates (day-1) for sapwood and fineroots when not
    specified via SpParams (`SRsapwood` and `SRfineroot`). Defaults are
    equivalent to 9%, 5% and 50% annual turnover for gymnosperm sapwood,
    angiosperm sapwood and fine roots, respectively.

  - `maximumRelativeGrowthRates [= c(leaf = 0.09, cambium = 0.005, sapwood = 0.002, fineroot = 0.1)]`:
    Default maximum relative growth rates for leaves (m2 leaf ·cm-2
    sapwood· day-1), tree sapwood (cm2 sapwood· cm-1 cambium · day-1),
    shrub sapwood (cm2 sapwood ·cm-2 sapwood· day-1) and fine roots (g
    dw · g dw -1 · day -1) when not specified via SpParams
    (`RGRleafmax`, `RGRcambiummax` , `RGRsapwoodmax` and
    `RGRfinerootmax`, respectively).

  - `mortalityMode [= "density/deterministic"]`: String describing how
    mortality is applied. Current accepted values are combinations of
    "cohort" vs "density" (for whole-cohort mortality vs reduction of
    stem density) and "deterministic" vs. "stochastic".

  - `mortalityBaselineRate [= 0.0015]`: Default deterministic proportion
    or probability specifying the baseline reduction of cohort's density
    occurring in a year (for
    `mortalityMode = "density/deterministic" or "density/stochastic").`

  - `mortalityRelativeSugarThreshold [= 0.4]`: Threshold of stem sugar
    concentration relative to the equilibrium sugar concentration,
    resulting in an increased starvation mortality rate/probability
    whenever levels are below.

  - `mortalityRWCThreshold [= 0.4]`: Threshold of stem relative water
    content resulting in increased mortality rate/probability whenever
    levels are below.

  - `recrTreeDBH [= 1]`: Default DBH (cm) for recruited trees (when
    species parameter `RecrTreeDBH` is missing).

  - `recrTreeDensity [= 3000]`: Default density (ind·ha-1) for recruited
    trees (when species parameter `RecrTreeDensity` is missing).

  - `ingrowthTreeDBH [= 7.5]`: Default DBH (cm) for ingrowth trees (when
    species parameter `IngrowthTreeDBH` is missing).

  - `ingrowthTreeDensity [= 127]`: Default density (ind·ha-1) for
    ingrowth trees (when species parameter `IngrowthTreeDensity` is
    missing).

  **Forest dynamics** (function
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md)):

  - `allowSeedBankDynamics [= TRUE]`: Boolean flag to indicate that seed
    production and seed bank dynamics is simulated.

  - `allowRecruitment [= TRUE]`: Boolean flag to indicate that
    recruitment from seeds is allowed.

  - `allowResprouting [= TRUE]`: Boolean flag to indicate that
    resprouting is allowed.

  - `recruitmentMode [= "stochastic"]`: String describing how
    recruitment from seeds is applied. Current accepted values are
    "deterministic" or "stochastic".

  - `removeEmptyCohorts [= TRUE]`: Boolean flag to indicate the removal
    of cohorts whose density is too low.

  - `minimumTreeCohortDensity [= 1]`: Threshold of tree density
    resulting in cohort removal.

  - `minimumShrubCohortCover [= 0.01]`: Threshold of shrub cover
    resulting in cohort removal.

  - `dynamicallyMergeCohorts [= TRUE]`: Boolean flag to indicate that
    cohorts should be merged when possible. This option speeds up
    calculations but results in a loss of cohort identity and
    reinitialization of many state variables.

  - `keepCohortsWithID [= TRUE]`: Boolean flag to indicate that cohorts
    having a non-missing value in a column `"ID"` (if present) should
    not be merged or removed.

  - `seedRain [= NULL]`: Vector of species names whose seed rain is to
    be added to seed bank, regardless of local seed production.

  - `seedProductionTreeHeight [= 300]`: Default minimum tree height for
    producing seeds (when species parameter `SeedProductionHeight` is
    missing).

  - `seedProductionShrubHeight [= 30]`: Default minimum shrub height for
    producing seeds (when species parameter `SeedProductionHeight` is
    missing).

  - `probRecr [= 0.05]`: Default annual probability of seed-recruitment
    (when species parameter `ProbRecr` is missing).

  - `minTempRecr [= 0]`: Default threshold of minimum average
    temperature of the coldest month necessary for recruiting from seeds
    (when species parameter `MinTempRecr` is missing).

  - `minMoistureRecr [= 0.3]`: Default threshold of minimum moisture
    index (annual precipitation over annual ETP) necessary for
    seed-recruiting (when species parameter `MinMoistureRecr` is
    missing).

  - `minFPARRecr [= 10]`: Default threshold of minimum fraction of PAR
    (in %) reaching the ground necessary for recruiting (when species
    parameter `MinFPARRecr` is missing).

  - `recrTreeHeight [= 620]`: Default height (cm) for recruited trees
    (when species parameter `RecrTreeHeight` is missing).

  - `recrShrubCover [= 1]`: Default cover (%) for shrubs recruited from
    seed (when species parameter `RecrShrubCover` is missing).

  - `recrShrubHeight [= 25]`: Default height (cm) for recruited shrubs
    (when species parameter `RecrShrubHeight` is missing).

  - `recrTreeZ50 [= 100]`: Default value for Z50 (mm) in seed-recruited
    trees (when species parameter `RecrZ50` is missing).

  - `recrShrubZ50 [= 50]`: Default value for Z50 (mm) in seed-recruited
    shrubs (when species parameter `RecrZ50` is missing).

  - `recrTreeZ95 [= 1000]`: Default value for Z95 (mm) in seed-recruited
    trees (when species parameter `RecrZ50` is missing).

  - `recrShrubZ50 [= 500]`: Default value for Z95 (mm) in seed-recruited
    shrubs (when species parameter `RecrZ50` is missing).

## Details

The function returns a list with default parameters. Users can change
those defaults that need to be set to other values and use the list as
input for model functions. The relevant parameters are different for
each model function.

## See also

[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
[`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md)

## Author

Miquel De Cáceres Ainsa, CREAF
