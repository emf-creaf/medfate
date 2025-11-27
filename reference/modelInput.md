# Input for simulation models

Functions `spwbInput()` and `growthInput()` take an object of class
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md) and
a soil data input to create input objects for simulation functions
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) (or
[`pwb`](https://emf-creaf.github.io/medfate/reference/pwb.md)) and
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
respectively.

## Usage

``` r
spwbInput(x, soil, SpParams, control)

growthInput(x, soil, SpParams, control)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- soil:

  An object of class
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
  containing soil parameters per soil layer.

- SpParams:

  A data frame with species parameters (see
  [`SpParamsDefinition`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  and
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- control:

  A list with default control parameters (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

## Value

Function `spwbInput()` returns a list of class `spwbInput` with the
following elements (rows of data frames are identified as specified by
function
[`plant_ID`](https://emf-creaf.github.io/medfate/reference/plant_values.md)):

- `control`: List with control parameters (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- `soil`: A data frame with initialized soil parameters (see
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)).

- `snowpack`: The amount of snow (in mm) in the snow pack over the soil.

- `canopy`: A list of stand-level state variables.

- `cohorts`: A data frame with cohort information, with columns `SP` and
  `Name`.

- `above`: A data frame with columns such as `H`, `CR` and `LAI_live`
  (see internal function `forest2aboveground`).

- `below`: A data frame with columns `Z50`, `Z95` and `Z100`.

- `belowLayers`: A list. If `control$transpirationMode = "Granier"` it
  contains elements:

  - `V`: A matrix with the proportion of fine roots of each cohort (in
    rows) in each soil layer (in columns).

  - `L`: A matrix with the length of coarse roots of each cohort (in
    rows) in each soil layer (in columns).

  - `Wpool`: A matrix with the soil moisture relative to field capacity
    around the rhizosphere of each cohort (in rows) in each soil layer
    (in columns).

  If `control$transpirationMode = "Sperry"` or
  `control$transpirationMode = "Sureau"` there are the following
  additional elements:

  - `VGrhizo_kmax`: A matrix with maximum rhizosphere conductance values
    of each cohort (in rows) in each soil layer (in columns).

  - `VGroot_kmax`: A matrix with maximum root xylem conductance values
    of each cohort (in rows) in each soil layer (in columns).

  - `RhizoPsi`: A matrix with the water potential around the rhizosphere
    of each cohort (in rows) in each soil layer (in columns).

- `paramsPhenology`: A data frame with leaf phenology parameters:

  - `PhenologyType`: Leaf phenology type.

  - `LeafDuration`: Leaf duration (in years).

  - `Sgdd`: Degree days needed for leaf budburst (for winter decideous
    species).

  - `Tbgdd`: Base temperature for the calculation of degree days to leaf
    budburst.

  - `Ssen`: Degree days corresponding to leaf senescence.

  - `Phsen`: Photoperiod corresponding to start counting senescence
    degree-days.

  - `Tbsen`: Base temperature for the calculation of degree days to leaf
    senescence.

- `paramsAnatomy`: A data frame with plant anatomy parameters for each
  cohort:

  - `Hmax`: Maximum plant height (cm).

  - `Hmed`: Median plant height (cm).

  - `Al2As`: Leaf area to sapwood area ratio (in m2·m-2).

  - `Ar2Al`: Fine root area to leaf area ratio (in m2·m-2).

  - `SLA`: Specific leaf area (mm2/mg = m2/kg).

  - `LeafWidth`: Leaf width (in cm).

  - `LeafDensity`: Density of leaf tissue (dry weight over volume).

  - `WoodDensity`: Density of wood tissue (dry weight over volume).

  - `FineRootDensity`: Density of fine root tissue (dry weight over
    volume).

  - `SRL`: Specific Root length (cm·g-1).

  - `RLD`: Root length density (cm·cm-3).

  - `r635`: Ratio between the weight of leaves plus branches and the
    weight of leaves alone for branches of 6.35 mm.

- `paramsInterception`: A data frame with rain interception and light
  extinction parameters for each cohort:

  - `kPAR`: PAR extinction coefficient.

  - `g`: Canopy water retention capacity per LAI unit (mm/LAI).

  If `control$transpirationMode = "Sperry"` or
  `control$transpirationMode = "Sureau"` additional columns are:

  - `gammaSWR`: Reflectance (albedo) coefficient for SWR .

  - `alphaSWR`: Absorbance coefficient for SWR .

- `paramsTranspiration`: A data frame with parameters for transpiration
  and photosynthesis. If `control$transpirationMode = "Granier"`,
  columns are:

  - `Gswmin`: Minimum stomatal conductance to water vapor (in mol
    H2O·m-2·s-1).

  - `Tmax_LAI`: Coefficient relating LAI with the ratio of maximum
    transpiration over potential evapotranspiration.

  - `Tmax_LAIsq`: Coefficient relating squared LAI with the ratio of
    maximum transpiration over potential evapotranspiration.

  - `Psi_Extract`: Water potential corresponding to 50% relative
    transpiration (in MPa).

  - `Exp_Extract`: Parameter of the Weibull function regulating
    transpiration reduction.

  - `VCstem_c`, `VCstem_d`: Parameters of the stem xylem vulnerability
    curve (Weibull).

  - `WUE`: Daily water use efficiency (gross photosynthesis over
    transpiration) under no light, water or CO2 limitations and VPD =
    1kPa (g C/mm water).

  - `WUE_par`: Coefficient regulating the influence of % PAR on gross
    photosynthesis.

  - `WUE_co2`: Coefficient regulating the influence of atmospheric CO2
    concentration on gross photosynthesis.

  - `WUE_vpd`: Coefficient regulating the influence of vapor pressure
    deficit (VPD) on gross photosynthesis.

  If `control$transpirationMode = "Sperry"` columns are:

  - `Gswmin`: Minimum stomatal conductance to water vapor (in mol
    H2O·m-2·s-1).

  - `Gswmax`: Maximum stomatal conductance to water vapor (in mol
    H2O·m-2·s-1).

  - `Vmax298`: Maximum Rubisco carboxilation rate at 25ºC (in micromol
    CO2·s-1·m-2).

  - `Jmax298`: Maximum rate of electron transport at 25ºC (in micromol
    photons·s-1·m-2).

  - `Kmax_stemxylem`: Sapwood-specific hydraulic conductivity of stem
    xylem (in kg H2O·s-1·m-1·MPa-1).

  - `Kmax_rootxylem`: Sapwood-specific hydraulic conductivity of root
    xylem (in kg H2O·s-1·m-1·MPa-1).

  - `VCleaf_kmax`: Maximum leaf hydraulic conductance (in mmol
    H2O·s-1·m-2·MPa-1).

  - `VCleaf_c`, `VCleaf_d`: Parameters of the leaf vulnerability curve
    (Weibull).

  - `VCstem_kmax`: Maximum stem xylem conductance (in mmol
    H2O·s-1·m-2·MPa-1).

  - `VCstem_c`, `VCstem_d`: Parameters of the stem xylem vulnerability
    curve (Weibull).

  - `VCroot_c`, `VCroot_d`: Parameters of the root xylem vulnerability
    curve (Weibull).

  - `Plant_kmax`: Maximum whole-plant conductance (in mmol
    H2O·s-1·m-2·MPa-1).

  - `FR_leaf`, `FR_stem`, `FR_root`: Fraction of whole-plant resistance
    corresponding to each segment.

  If `control$transpirationMode = "Sureau"` columns are:

  - `Gswmin`: Minimum stomatal conductance to water vapor (in mol
    H2O·m-2·s-1).

  - `Gswmax`: Maximum stomatal conductance to water vapor (in mol
    H2O·m-2·s-1).

  - `Gsw_AC_slope`: Slope of the Gsw vs Ac/Cs relationship (see
    [`photo_photosynthesisBaldocchi`](https://emf-creaf.github.io/medfate/reference/photo.md)).

  - `Gs_P50`, `Gs_slope`: Parameters of the curve describing the
    decrease in stomatal conductance as a function of leaf water
    potential (sigmoid).

  - `Vmax298`: Maximum Rubisco carboxylation rate at 25ºC (in micromol
    CO2·s-1·m-2).

  - `Jmax298`: Maximum rate of electron transport at 25ºC (in micromol
    photons·s-1·m-2).

  - `Kmax_stemxylem`: Sapwood-specific hydraulic conductivity of stem
    xylem (in kg H2O·s-1·m-1·MPa-1).

  - `Kmax_rootxylem`: Sapwood-specific hydraulic conductivity of root
    xylem (in kg H2O·s-1·m-1·MPa-1).

  - `VCleaf_kmax`: Maximum leaf hydraulic conductance (in mmol
    H2O·s-1·m-2·MPa-1).

  - `VCleaf_c`, `VCleaf_d`: Parameters of the leaf vulnerability curve
    (Weibull).

  - `VCleaf_P50`, `VCleaf_slope`: Parameters of the leaf vulnerability
    curve (sigmoid).

  - `VCstem_kmax`: Maximum stem xylem conductance (in mmol
    H2O·s-1·m-2·MPa-1).

  - `VCstem_c`, `VCstem_d`: Parameters of the stem xylem vulnerability
    curve (Weibull).

  - `VCstem_P50`, `VCstem_slope`: Parameters of the stem xylem
    vulnerability curve (sigmoid).

  - `VCroot_c`, `VCroot_d`: Parameters of the root xylem vulnerability
    curve (Weibull).

  - `VCroot_P50`, `VCroot_slope`: Parameters of the root xylem
    vulnerability curve (sigmoid).

  - `Plant_kmax`: Maximum whole-plant conductance (in mmol
    H2O·s-1·m-2·MPa-1).

  - `FR_leaf`, `FR_stem`, `FR_root`: Fraction of whole-plant resistance
    corresponding to each segment.

- `paramsWaterStorage`: A data frame with plant water storage parameters
  for each cohort:

  - `LeafPI0`: Osmotic potential at full turgor of leaves (MPa).

  - `LeafEPS`: Modulus of elasticity (capacity of the cell wall to
    resist changes in volume in response to changes in turgor) of leaves
    (MPa).

  - `LeafAF`: Apoplastic fraction (proportion of water outside the
    living cells) in leaves.

  - `Vleaf`: Storage water capacity in leaves, per leaf area (L/m2).

  - `StemPI0`: Osmotic potential at full turgor of symplastic xylem
    tissue (MPa).

  - `StemEPS`: Modulus of elasticity (capacity of the cell wall to
    resist changes in volume in response to changes in turgor) of
    symplastic xylem tissue (Mpa).

  - `StemAF`: Apoplastic fraction (proportion of water outside the
    living cells) in stem xylem.

  - `Vstem`: Storage water capacity in sapwood, per leaf area (L/m2).

- `internalPhenology` and `internalWater`: data frames to store internal
  state variables.

- `internalFCCS`: A data frame with fuel characteristics, according to
  [`fuel_FCCS`](https://emf-creaf.github.io/medfate/reference/fuel_properties.md)
  (only if `fireHazardResults = TRUE`, in the control list).

Function `growthInput()` returns a list of class `growthInput` with the
same elements as `spwbInput`, but with additional information.

- Element `above` includes the following additional columns:

  - `LA_live`: Live leaf area per individual (m2/ind).

  - `LA_dead`: Dead leaf area per individual (m2/ind).

  - `SA`: Live sapwood area per individual (cm2/ind).

- `paramsGrowth`: A data frame with growth parameters for each cohort:

  - `RERleaf`: Maintenance respiration rates (at 20ºC) for leaves (in g
    gluc·g dry-1·day-1).

  - `RERsapwood`: Maintenance respiration rates (at 20ºC) for sapwood
    (in g gluc·g dry-1·day-1).

  - `RERfineroot`: Maintenance respiration rates (at 20ºC) for fine
    roots (in g gluc·g dry-1·day-1).

  - `CCleaf`: Leaf construction costs (in g gluc·g dry-1).

  - `CCsapwood`: Sapwood construction costs (in g gluc·g dry-1).

  - `CCfineroot`: Fine root construction costs (in g gluc·g dry-1).

  - `RGRleafmax`: Maximum leaf relative growth rate (in m2·cm-2·day-1).

  - `RGRsapwoodmax`: Maximum sapwood relative growth rate (in
    cm2·cm-2·day-1).

  - `RGRfinerootmax`: Maximum fine root relative growth rate (in g dry·g
    dry-1·day-1).

  - `SRsapwood`: Sapwood daily senescence rate (in day-1).

  - `SRfineroot`: Fine root daily senescence rate (in day-1).

  - `RSSG`: Minimum relative starch for sapwood growth (proportion).

  - `fHDmin`: Minimum value of the height-to-diameter ratio
    (dimensionless).

  - `fHDmax`: Maximum value of the height-to-diameter ratio
    (dimensionless).

  - `WoodC`: Wood carbon content per dry weight (g C /g dry).

- `paramsMortalityRegeneration`: A data frame with
  mortality/regeneration parameters for each cohort:

  - `MortalityBaselineRate`: Deterministic proportion or probability
    specifying the baseline reduction of cohort's density occurring in a
    year.

  - `SurvivalModelStep`: Time step in years of the empirical survival
    model depending on stand basal area (e.g. 10).

  - `SurvivalB0`: Intercept of the logistic baseline survival model
    depending on stand basal area.

  - `SurvivalB1`: Slope of the logistic baseline survival model
    depending on stand basal area.

  - `RecrTreeDensity`: Density of tree recruits from seeds.

  - `IngrowthTreeDensity`: Density of trees reaching ingrowth DBH.

  - `RecrTreeDBH`: DBH for tree recruits from seeds or resprouting (e.g.
    1 cm).

  - `IngrowthTreeDBH`: Ingrowth DBH for trees (e.g. 7.5 cm).

- `paramsAllometry`: A data frame with allometric parameters for each
  cohort:

  - `Aash`: Regression coefficient relating the square of shrub height
    with shrub area.

  - `Absh`, `Bbsh`: Allometric coefficients relating phytovolume with
    dry weight of shrub individuals.

  - `Acr`, `B1cr`, `B2cr`, `B3cr`, `C1cr`, `C2cr`: Regression
    coefficients used to calculate crown ratio of trees.

  - `Acw`, `Bcw`: Regression coefficients used to calculated crown width
    of trees.

- `internalAllocation`: A data frame with internal allocation variables
  for each cohort:

  - `allocationTarget`: Value of the allocation target variable.

  - `leafAreaTarget`: Target leaf area (m2) per individual.

  - `sapwoodAreaTarget`: Target sapwood area (cm2) per individual.

  - `fineRootBiomassTarget`: Target fine root biomass (g dry) per
    individual.

  - `crownBudPercent`: Percentage of the crown with buds.

- `internalCarbon`: A data frame with the concentration (mol·gluc·l-1)
  of metabolic and storage carbon compartments for leaves and sapwood.

- `internalMortality`: A data frame to store the cumulative mortality
  (density for trees and cover for shrubs) predicted during the
  simulation, also distinguishing mortality due to starvation or
  dessication.

## Details

Functions `spwbInput()` and `growthInput()` initialize inputs
differently depending on control parameters.

*IMPORTANT NOTE*: Older function names
[`forest2spwbInput`](https://emf-creaf.github.io/medfate/reference/forest2aboveground.md)
and
[`forest2growthInput`](https://emf-creaf.github.io/medfate/reference/forest2aboveground.md)
are now deprecated, but they can still be used for back-compatibility.

## See also

[`resetInputs`](https://emf-creaf.github.io/medfate/reference/resetInputs.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md),
[`defaultSoilParams`](https://emf-creaf.github.io/medfate/reference/defaultSoilParams.md),
[`plant_ID`](https://emf-creaf.github.io/medfate/reference/plant_values.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example plot plant data
data(exampleforest)

# Example of aboveground parameters taken from a forest
# described using LAI and crown ratio
data(exampleforest2)

#Default species parameterization
data(SpParamsMED)


# Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

# Initialize control parameters using 'Granier' transpiration mode
control <- defaultControl("Granier")

# Prepare spwb input
spwbInput(exampleforest, examplesoil, SpParamsMED, control)
#> $control
#> $control$fillMissingRootParams
#> [1] TRUE
#> 
#> $control$fillMissingSpParams
#> [1] TRUE
#> 
#> $control$fillMissingWithGenusParams
#> [1] TRUE
#> 
#> $control$verbose
#> [1] TRUE
#> 
#> $control$subdailyResults
#> [1] FALSE
#> 
#> $control$standResults
#> [1] TRUE
#> 
#> $control$soilResults
#> [1] TRUE
#> 
#> $control$soilPoolResults
#> [1] FALSE
#> 
#> $control$snowResults
#> [1] TRUE
#> 
#> $control$plantResults
#> [1] TRUE
#> 
#> $control$labileCarbonBalanceResults
#> [1] TRUE
#> 
#> $control$plantStructureResults
#> [1] TRUE
#> 
#> $control$growthMortalityResults
#> [1] TRUE
#> 
#> $control$leafResults
#> [1] TRUE
#> 
#> $control$temperatureResults
#> [1] TRUE
#> 
#> $control$fireHazardResults
#> [1] FALSE
#> 
#> $control$fireHazardStandardWind
#> [1] NA
#> 
#> $control$fireHazardStandardDFMC
#> [1] NA
#> 
#> $control$transpirationMode
#> [1] "Granier"
#> 
#> $control$soilDomains
#> [1] "buckets"
#> 
#> $control$rhizosphereOverlap
#> [1] "total"
#> 
#> $control$truncateRootDistribution
#> [1] FALSE
#> 
#> $control$fullRhizosphereOverlapConductivity
#> [1] 0.01
#> 
#> $control$soilFunctions
#> [1] "VG"
#> 
#> $control$VG_PTF
#> [1] "Toth"
#> 
#> $control$ndailysteps
#> [1] 24
#> 
#> $control$max_nsubsteps_soil
#> [1] 300
#> 
#> $control$defaultWindSpeed
#> [1] 2.5
#> 
#> $control$defaultCO2
#> [1] 386
#> 
#> $control$defaultRainfallIntensityPerMonth
#>  [1] 1.5 1.5 1.5 1.5 1.5 1.5 5.6 5.6 5.6 5.6 5.6 1.5
#> 
#> $control$leafPhenology
#> [1] TRUE
#> 
#> $control$bareSoilEvaporation
#> [1] TRUE
#> 
#> $control$unlimitedSoilWater
#> [1] FALSE
#> 
#> $control$interceptionMode
#> [1] "Gash1995"
#> 
#> $control$infiltrationMode
#> [1] "GreenAmpt1911"
#> 
#> $control$infiltrationCorrection
#> [1] 5
#> 
#> $control$unfoldingDD
#> [1] 300
#> 
#> $control$verticalLayerSize
#> [1] 100
#> 
#> $control$windMeasurementHeight
#> [1] 200
#> 
#> $control$segmentedXylemVulnerability
#> [1] TRUE
#> 
#> $control$stemCavitationRecovery
#> [1] "annual"
#> 
#> $control$leafCavitationRecovery
#> [1] "annual"
#> 
#> $control$lfmcComponent
#> [1] "fine"
#> 
#> $control$hydraulicRedistributionFraction
#> [1] 0.1
#> 
#> $control$nsubsteps_canopy
#> [1] 3600
#> 
#> $control$taper
#> [1] TRUE
#> 
#> $control$multiLayerBalance
#> [1] FALSE
#> 
#> $control$sapFluidityVariation
#> [1] TRUE
#> 
#> $control$TPhase_gmin
#> [1] 37.5
#> 
#> $control$Q10_1_gmin
#> [1] 1.2
#> 
#> $control$Q10_2_gmin
#> [1] 4.8
#> 
#> $control$rootRadialConductance
#> [1] 4
#> 
#> $control$averageFracRhizosphereResistance
#> [1] 0.15
#> 
#> $control$thermalCapacityLAI
#> [1] 1e+06
#> 
#> $control$boundaryLayerSize
#> [1] 2000
#> 
#> $control$cavitationRecoveryMaximumRate
#> [1] 0.05
#> 
#> $control$sunlitShade
#> [1] TRUE
#> 
#> $control$numericParams
#> $control$numericParams$maxNsteps
#> [1] 400
#> 
#> $control$numericParams$ntrial
#> [1] 200
#> 
#> $control$numericParams$psiTol
#> [1] 1e-04
#> 
#> $control$numericParams$ETol
#> [1] 1e-07
#> 
#> 
#> $control$leafCavitationEffects
#> [1] FALSE
#> 
#> $control$stemCavitationEffects
#> [1] TRUE
#> 
#> $control$stomatalSubmodel
#> [1] "Baldocchi"
#> 
#> $control$plantCapacitance
#> [1] TRUE
#> 
#> $control$cavitationFlux
#> [1] TRUE
#> 
#> $control$leafCuticularTranspiration
#> [1] TRUE
#> 
#> $control$stemCuticularTranspiration
#> [1] FALSE
#> 
#> $control$C_SApoInit
#> [1] 2e-05
#> 
#> $control$C_LApoInit
#> [1] 1e-05
#> 
#> $control$k_SSym
#> [1] 0.26
#> 
#> $control$fractionLeafSymplasm
#> [1] 0.5
#> 
#> $control$gs_NightFrac
#> [1] 0.05
#> 
#> $control$JarvisPAR
#> [1] 0.003
#> 
#> $control$fTRBToLeaf
#> [1] 0.8
#> 
#> $control$subdailyCarbonBalance
#> [1] FALSE
#> 
#> $control$allowDessication
#> [1] TRUE
#> 
#> $control$allowStarvation
#> [1] TRUE
#> 
#> $control$sinkLimitation
#> [1] TRUE
#> 
#> $control$shrubDynamics
#> [1] TRUE
#> 
#> $control$herbDynamics
#> [1] TRUE
#> 
#> $control$allocationStrategy
#> [1] "Al2As"
#> 
#> $control$phloemConductanceFactor
#> [1] 0.2
#> 
#> $control$nonSugarConcentration
#> [1] 0.25
#> 
#> $control$equilibriumOsmoticConcentration
#> $control$equilibriumOsmoticConcentration$leaf
#> [1] 0.8
#> 
#> $control$equilibriumOsmoticConcentration$sapwood
#> [1] 0.6
#> 
#> 
#> $control$minimumRelativeStarchForGrowth
#> [1] 0.5
#> 
#> $control$constructionCosts
#> $control$constructionCosts$leaf
#> [1] 1.5
#> 
#> $control$constructionCosts$sapwood
#> [1] 1.47
#> 
#> $control$constructionCosts$fineroot
#> [1] 1.3
#> 
#> 
#> $control$senescenceRates
#> $control$senescenceRates$sapwood
#> [1] 0.000135
#> 
#> $control$senescenceRates$fineroot
#> [1] 0.001897231
#> 
#> 
#> $control$maximumRelativeGrowthRates
#> $control$maximumRelativeGrowthRates$leaf
#> [1] 0.09
#> 
#> $control$maximumRelativeGrowthRates$cambium
#> [1] 0.0025
#> 
#> $control$maximumRelativeGrowthRates$sapwood
#> [1] 0.002
#> 
#> $control$maximumRelativeGrowthRates$fineroot
#> [1] 0.1
#> 
#> 
#> $control$mortalityMode
#> [1] "density/deterministic"
#> 
#> $control$mortalityBaselineRate
#> [1] 0.0015
#> 
#> $control$mortalityRelativeSugarThreshold
#> [1] 0.4
#> 
#> $control$mortalityRWCThreshold
#> [1] 0.4
#> 
#> $control$recrTreeDBH
#> [1] 1
#> 
#> $control$recrTreeDensity
#> [1] 3000
#> 
#> $control$ingrowthTreeDBH
#> [1] 7.5
#> 
#> $control$ingrowthTreeDensity
#> [1] 127
#> 
#> $control$allowSeedBankDynamics
#> [1] TRUE
#> 
#> $control$allowRecruitment
#> [1] TRUE
#> 
#> $control$allowResprouting
#> [1] TRUE
#> 
#> $control$recruitmentMode
#> [1] "stochastic"
#> 
#> $control$removeEmptyCohorts
#> [1] TRUE
#> 
#> $control$minimumTreeCohortDensity
#> [1] 1
#> 
#> $control$minimumShrubCohortCover
#> [1] 0.01
#> 
#> $control$dynamicallyMergeCohorts
#> [1] TRUE
#> 
#> $control$keepCohortsWithObsID
#> [1] FALSE
#> 
#> $control$seedRain
#> NULL
#> 
#> $control$seedProductionTreeHeight
#> [1] 300
#> 
#> $control$seedProductionShrubHeight
#> [1] 30
#> 
#> $control$probRecr
#> [1] 0.05
#> 
#> $control$minTempRecr
#> [1] 0
#> 
#> $control$minMoistureRecr
#> [1] 0.3
#> 
#> $control$minFPARRecr
#> [1] 10
#> 
#> $control$recrTreeHeight
#> [1] 620
#> 
#> $control$recrShrubCover
#> [1] 1
#> 
#> $control$recrShrubHeight
#> [1] 25
#> 
#> $control$recrTreeZ50
#> [1] 100
#> 
#> $control$recrShrubZ50
#> [1] 50
#> 
#> $control$recrTreeZ95
#> [1] 1000
#> 
#> $control$recrShrubZ95
#> [1] 500
#> 
#> 
#> $soil
#>   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha
#> 1    300   25   25 Silt loam NA       NA 1.5  25 0.0485 5401.471 89.16112
#> 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112
#> 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112
#> 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112
#>       VG_n VG_theta_res VG_theta_sat W Temp
#> 1 1.303861        0.041     0.423715 1   NA
#> 2 1.303861        0.041     0.423715 1   NA
#> 3 1.303861        0.041     0.423715 1   NA
#> 4 1.303861        0.041     0.423715 1   NA
#> 
#> $snowpack
#> [1] 0
#> 
#> $canopy
#>    zlow zmid  zup LAIlive LAIexpanded LAIdead Tair Cair VPair
#> 1     0   50  100      NA          NA      NA   NA   NA    NA
#> 2   100  150  200      NA          NA      NA   NA   NA    NA
#> 3   200  250  300      NA          NA      NA   NA   NA    NA
#> 4   300  350  400      NA          NA      NA   NA   NA    NA
#> 5   400  450  500      NA          NA      NA   NA   NA    NA
#> 6   500  550  600      NA          NA      NA   NA   NA    NA
#> 7   600  650  700      NA          NA      NA   NA   NA    NA
#> 8   700  750  800      NA          NA      NA   NA   NA    NA
#> 9   800  850  900      NA          NA      NA   NA   NA    NA
#> 10  900  950 1000      NA          NA      NA   NA   NA    NA
#> 11 1000 1050 1100      NA          NA      NA   NA   NA    NA
#> 12 1100 1150 1200      NA          NA      NA   NA   NA    NA
#> 13 1200 1250 1300      NA          NA      NA   NA   NA    NA
#> 14 1300 1350 1400      NA          NA      NA   NA   NA    NA
#> 15 1400 1450 1500      NA          NA      NA   NA   NA    NA
#> 16 1500 1550 1600      NA          NA      NA   NA   NA    NA
#> 17 1600 1650 1700      NA          NA      NA   NA   NA    NA
#> 18 1700 1750 1800      NA          NA      NA   NA   NA    NA
#> 19 1800 1850 1900      NA          NA      NA   NA   NA    NA
#> 20 1900 1950 2000      NA          NA      NA   NA   NA    NA
#> 21 2000 2050 2100      NA          NA      NA   NA   NA    NA
#> 22 2100 2150 2200      NA          NA      NA   NA   NA    NA
#> 23 2200 2250 2300      NA          NA      NA   NA   NA    NA
#> 24 2300 2350 2400      NA          NA      NA   NA   NA    NA
#> 25 2400 2450 2500      NA          NA      NA   NA   NA    NA
#> 26 2500 2550 2600      NA          NA      NA   NA   NA    NA
#> 27 2600 2650 2700      NA          NA      NA   NA   NA    NA
#> 28 2700 2750 2800      NA          NA      NA   NA   NA    NA
#> 
#> $herbLAI
#> [1] 0.1736369
#> 
#> $herbLAImax
#> [1] 0.252
#> 
#> $cohorts
#>         SP              Name
#> T1_148 148  Pinus halepensis
#> T2_168 168      Quercus ilex
#> S1_165 165 Quercus coccifera
#> 
#> $above
#>          H        CR   LAI_live LAI_expanded LAI_dead Age ObsID
#> T1_148 800 0.6605196 0.84874773   0.84874773        0  NA  <NA>
#> T2_168 660 0.6055642 0.70557382   0.70557382        0  NA  <NA>
#> S1_165  80 0.8032817 0.03062604   0.03062604        0  NA  <NA>
#> 
#> $below
#>        Z50  Z95 Z100
#> T1_148 100  600   NA
#> T2_168 300 1000   NA
#> S1_165 200 1000   NA
#> 
#> $belowLayers
#> $belowLayers$V
#>                1         2          3           4
#> T1_148 0.8604899 0.1194556 0.01511005 0.004944476
#> T2_168 0.5008953 0.4505941 0.04064831 0.007862284
#> S1_165 0.6799879 0.2737911 0.03567632 0.010544678
#> 
#> $belowLayers$L
#>          1   2    3    4
#> T1_148 150 650 1500 3000
#> T2_168 150 650 1500 3000
#> S1_165 150 650 1500 3000
#> 
#> $belowLayers$Wpool
#>        1 2 3 4
#> T1_148 1 1 1 1
#> T2_168 1 1 1 1
#> S1_165 1 1 1 1
#> 
#> 
#> $paramsPhenology
#>             PhenologyType LeafDuration t0gdd  Sgdd Tbgdd  Ssen Phsen Tbsen xsen
#> T1_148 oneflush-evergreen     2.536875  50.0 200.0  0.00  8268  12.5  28.5    2
#> T2_168 oneflush-evergreen     2.183837  54.5 240.7  4.34 10178  12.5  28.5    2
#> S1_165 oneflush-evergreen     1.250000  54.5 240.7  4.34 10178  12.5  28.5    2
#>        ysen
#> T1_148    2
#> T2_168    2
#> S1_165    2
#> 
#> $paramsAnatomy
#>           Al2As Ar2Al      SLA LeafDensity WoodDensity FineRootDensity      SRL
#> T1_148 1317.523     1 5.140523   0.2982842   0.6077016       0.2982842 3172.572
#> T2_168 3908.823     1 6.340000   0.4893392   0.9008264       0.4893392 4398.812
#> S1_165 4189.325     1 4.980084   0.3709679   0.4389106       0.3709679 4398.812
#>        RLD     r635
#> T1_148  10 1.964226
#> T2_168  10 1.805872
#> S1_165  10 2.289452
#> 
#> $paramsInterception
#>        kPAR      kSWR   g
#> T1_148 0.50 0.3703704 1.0
#> T2_168 0.55 0.4074074 0.5
#> S1_165 0.55 0.4074074 0.5
#> 
#> $paramsTranspiration
#>             Gswmin  Tmax_LAI   Tmax_LAIsq Psi_Extract Exp_Extract  VCleaf_c
#> T1_148 0.003086667 0.1869849 -0.008372458  -0.9218219    1.504542 11.137050
#> T2_168 0.004473333 0.1251027 -0.005601615  -1.9726871    1.149052  1.339370
#> S1_165 0.010455247 0.1340000 -0.006000000  -2.1210726    1.300000  2.254991
#>         VCleaf_d  VCstem_c  VCstem_d      WUE   WUE_par     WUE_co2    WUE_vpd
#> T1_148 -2.380849 12.709999 -5.290000 8.525550 0.5239136 0.002586327 -0.2647169
#> T2_168 -2.582279  3.560000 -7.720000 8.968208 0.1412266 0.002413091 -0.5664879
#> S1_165 -3.133381  3.095442 -7.857378 7.900000 0.3643000 0.002757000 -0.4636000
#> 
#> $paramsWaterStorage
#>           maxFMC maxMCleaf maxMCstem   LeafPI0   LeafEPS LeafAF     Vleaf
#> T1_148 126.03063  151.9063  99.19498 -1.591429  8.918571 0.3525 0.5258525
#> T2_168  93.15304  131.4346  45.64970 -1.483333 19.260000 0.1700 0.2199087
#> S1_165  96.53441   47.3984 134.64052 -2.370000 17.230000 0.2400 0.4108968
#>          StemPI0   StemEPS    StemAF  Vsapwood
#> T1_148 -2.008039 13.256355 0.9236406 4.1638559
#> T2_168 -3.227438 46.420610 0.6238125 0.8135590
#> S1_165 -1.305868  6.297155 0.6238125 0.3177724
#> 
#> $internalPhenology
#>        gdd sen budFormation leafUnfolding leafSenescence leafDormancy phi
#> T1_148   0   0        FALSE         FALSE          FALSE        FALSE   0
#> T2_168   0   0        FALSE         FALSE          FALSE        FALSE   0
#> S1_165   0   0        FALSE         FALSE          FALSE        FALSE   0
#> 
#> $internalWater
#>        PlantPsi LeafPLC StemPLC
#> T1_148   -0.033       0       0
#> T2_168   -0.033       0       0
#> S1_165   -0.033       0       0
#> 
#> $internalLAIDistribution
#> $internalLAIDistribution$PrevLAIexpanded
#> [1] NA NA NA
#> 
#> $internalLAIDistribution$PrevLAIdead
#> [1] NA NA NA
#> 
#> $internalLAIDistribution$PARcohort
#> [1] 0 0 0
#> 
#> $internalLAIDistribution$live
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> $internalLAIDistribution$expanded
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> $internalLAIDistribution$dead
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> 
#> $internalFCCS
#> data frame with 0 columns and 0 rows
#> 
#> $version
#> [1] "4.8.5"
#> 
#> attr(,"class")
#> [1] "spwbInput" "list"     
                
# Prepare input for 'Sperry' transpiration mode
control <- defaultControl("Sperry")
spwbInput(exampleforest,examplesoil,SpParamsMED, control)
#> $control
#> $control$fillMissingRootParams
#> [1] TRUE
#> 
#> $control$fillMissingSpParams
#> [1] TRUE
#> 
#> $control$fillMissingWithGenusParams
#> [1] TRUE
#> 
#> $control$verbose
#> [1] TRUE
#> 
#> $control$subdailyResults
#> [1] FALSE
#> 
#> $control$standResults
#> [1] TRUE
#> 
#> $control$soilResults
#> [1] TRUE
#> 
#> $control$soilPoolResults
#> [1] FALSE
#> 
#> $control$snowResults
#> [1] TRUE
#> 
#> $control$plantResults
#> [1] TRUE
#> 
#> $control$labileCarbonBalanceResults
#> [1] TRUE
#> 
#> $control$plantStructureResults
#> [1] TRUE
#> 
#> $control$growthMortalityResults
#> [1] TRUE
#> 
#> $control$leafResults
#> [1] TRUE
#> 
#> $control$temperatureResults
#> [1] TRUE
#> 
#> $control$fireHazardResults
#> [1] FALSE
#> 
#> $control$fireHazardStandardWind
#> [1] NA
#> 
#> $control$fireHazardStandardDFMC
#> [1] NA
#> 
#> $control$transpirationMode
#> [1] "Sperry"
#> 
#> $control$soilDomains
#> [1] "buckets"
#> 
#> $control$rhizosphereOverlap
#> [1] "total"
#> 
#> $control$truncateRootDistribution
#> [1] FALSE
#> 
#> $control$fullRhizosphereOverlapConductivity
#> [1] 0.01
#> 
#> $control$soilFunctions
#> [1] "VG"
#> 
#> $control$VG_PTF
#> [1] "Toth"
#> 
#> $control$ndailysteps
#> [1] 24
#> 
#> $control$max_nsubsteps_soil
#> [1] 300
#> 
#> $control$defaultWindSpeed
#> [1] 2.5
#> 
#> $control$defaultCO2
#> [1] 386
#> 
#> $control$defaultRainfallIntensityPerMonth
#>  [1] 1.5 1.5 1.5 1.5 1.5 1.5 5.6 5.6 5.6 5.6 5.6 1.5
#> 
#> $control$leafPhenology
#> [1] TRUE
#> 
#> $control$bareSoilEvaporation
#> [1] TRUE
#> 
#> $control$unlimitedSoilWater
#> [1] FALSE
#> 
#> $control$interceptionMode
#> [1] "Gash1995"
#> 
#> $control$infiltrationMode
#> [1] "GreenAmpt1911"
#> 
#> $control$infiltrationCorrection
#> [1] 5
#> 
#> $control$unfoldingDD
#> [1] 300
#> 
#> $control$verticalLayerSize
#> [1] 100
#> 
#> $control$windMeasurementHeight
#> [1] 200
#> 
#> $control$segmentedXylemVulnerability
#> [1] TRUE
#> 
#> $control$stemCavitationRecovery
#> [1] "annual"
#> 
#> $control$leafCavitationRecovery
#> [1] "annual"
#> 
#> $control$lfmcComponent
#> [1] "fine"
#> 
#> $control$hydraulicRedistributionFraction
#> [1] 0.1
#> 
#> $control$nsubsteps_canopy
#> [1] 3600
#> 
#> $control$taper
#> [1] TRUE
#> 
#> $control$multiLayerBalance
#> [1] FALSE
#> 
#> $control$sapFluidityVariation
#> [1] TRUE
#> 
#> $control$TPhase_gmin
#> [1] 37.5
#> 
#> $control$Q10_1_gmin
#> [1] 1.2
#> 
#> $control$Q10_2_gmin
#> [1] 4.8
#> 
#> $control$rootRadialConductance
#> [1] 4
#> 
#> $control$averageFracRhizosphereResistance
#> [1] 0.15
#> 
#> $control$thermalCapacityLAI
#> [1] 1e+06
#> 
#> $control$boundaryLayerSize
#> [1] 2000
#> 
#> $control$cavitationRecoveryMaximumRate
#> [1] 0.05
#> 
#> $control$sunlitShade
#> [1] TRUE
#> 
#> $control$numericParams
#> $control$numericParams$maxNsteps
#> [1] 400
#> 
#> $control$numericParams$ntrial
#> [1] 200
#> 
#> $control$numericParams$psiTol
#> [1] 1e-04
#> 
#> $control$numericParams$ETol
#> [1] 1e-07
#> 
#> 
#> $control$leafCavitationEffects
#> [1] FALSE
#> 
#> $control$stemCavitationEffects
#> [1] TRUE
#> 
#> $control$stomatalSubmodel
#> [1] "Baldocchi"
#> 
#> $control$plantCapacitance
#> [1] TRUE
#> 
#> $control$cavitationFlux
#> [1] TRUE
#> 
#> $control$leafCuticularTranspiration
#> [1] TRUE
#> 
#> $control$stemCuticularTranspiration
#> [1] FALSE
#> 
#> $control$C_SApoInit
#> [1] 2e-05
#> 
#> $control$C_LApoInit
#> [1] 1e-05
#> 
#> $control$k_SSym
#> [1] 0.26
#> 
#> $control$fractionLeafSymplasm
#> [1] 0.5
#> 
#> $control$gs_NightFrac
#> [1] 0.05
#> 
#> $control$JarvisPAR
#> [1] 0.003
#> 
#> $control$fTRBToLeaf
#> [1] 0.8
#> 
#> $control$subdailyCarbonBalance
#> [1] FALSE
#> 
#> $control$allowDessication
#> [1] TRUE
#> 
#> $control$allowStarvation
#> [1] TRUE
#> 
#> $control$sinkLimitation
#> [1] TRUE
#> 
#> $control$shrubDynamics
#> [1] TRUE
#> 
#> $control$herbDynamics
#> [1] TRUE
#> 
#> $control$allocationStrategy
#> [1] "Al2As"
#> 
#> $control$phloemConductanceFactor
#> [1] 0.2
#> 
#> $control$nonSugarConcentration
#> [1] 0.25
#> 
#> $control$equilibriumOsmoticConcentration
#> $control$equilibriumOsmoticConcentration$leaf
#> [1] 0.8
#> 
#> $control$equilibriumOsmoticConcentration$sapwood
#> [1] 0.6
#> 
#> 
#> $control$minimumRelativeStarchForGrowth
#> [1] 0.5
#> 
#> $control$constructionCosts
#> $control$constructionCosts$leaf
#> [1] 1.5
#> 
#> $control$constructionCosts$sapwood
#> [1] 1.47
#> 
#> $control$constructionCosts$fineroot
#> [1] 1.3
#> 
#> 
#> $control$senescenceRates
#> $control$senescenceRates$sapwood
#> [1] 0.000135
#> 
#> $control$senescenceRates$fineroot
#> [1] 0.001897231
#> 
#> 
#> $control$maximumRelativeGrowthRates
#> $control$maximumRelativeGrowthRates$leaf
#> [1] 0.09
#> 
#> $control$maximumRelativeGrowthRates$cambium
#> [1] 0.0025
#> 
#> $control$maximumRelativeGrowthRates$sapwood
#> [1] 0.002
#> 
#> $control$maximumRelativeGrowthRates$fineroot
#> [1] 0.1
#> 
#> 
#> $control$mortalityMode
#> [1] "density/deterministic"
#> 
#> $control$mortalityBaselineRate
#> [1] 0.0015
#> 
#> $control$mortalityRelativeSugarThreshold
#> [1] 0.4
#> 
#> $control$mortalityRWCThreshold
#> [1] 0.4
#> 
#> $control$recrTreeDBH
#> [1] 1
#> 
#> $control$recrTreeDensity
#> [1] 3000
#> 
#> $control$ingrowthTreeDBH
#> [1] 7.5
#> 
#> $control$ingrowthTreeDensity
#> [1] 127
#> 
#> $control$allowSeedBankDynamics
#> [1] TRUE
#> 
#> $control$allowRecruitment
#> [1] TRUE
#> 
#> $control$allowResprouting
#> [1] TRUE
#> 
#> $control$recruitmentMode
#> [1] "stochastic"
#> 
#> $control$removeEmptyCohorts
#> [1] TRUE
#> 
#> $control$minimumTreeCohortDensity
#> [1] 1
#> 
#> $control$minimumShrubCohortCover
#> [1] 0.01
#> 
#> $control$dynamicallyMergeCohorts
#> [1] TRUE
#> 
#> $control$keepCohortsWithObsID
#> [1] FALSE
#> 
#> $control$seedRain
#> NULL
#> 
#> $control$seedProductionTreeHeight
#> [1] 300
#> 
#> $control$seedProductionShrubHeight
#> [1] 30
#> 
#> $control$probRecr
#> [1] 0.05
#> 
#> $control$minTempRecr
#> [1] 0
#> 
#> $control$minMoistureRecr
#> [1] 0.3
#> 
#> $control$minFPARRecr
#> [1] 10
#> 
#> $control$recrTreeHeight
#> [1] 620
#> 
#> $control$recrShrubCover
#> [1] 1
#> 
#> $control$recrShrubHeight
#> [1] 25
#> 
#> $control$recrTreeZ50
#> [1] 100
#> 
#> $control$recrShrubZ50
#> [1] 50
#> 
#> $control$recrTreeZ95
#> [1] 1000
#> 
#> $control$recrShrubZ95
#> [1] 500
#> 
#> 
#> $soil
#>   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha
#> 1    300   25   25 Silt loam NA       NA 1.5  25 0.0485 5401.471 89.16112
#> 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112
#> 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112
#> 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112
#>       VG_n VG_theta_res VG_theta_sat W Temp
#> 1 1.303861        0.041     0.423715 1   NA
#> 2 1.303861        0.041     0.423715 1   NA
#> 3 1.303861        0.041     0.423715 1   NA
#> 4 1.303861        0.041     0.423715 1   NA
#> 
#> $snowpack
#> [1] 0
#> 
#> $canopy
#>    zlow zmid  zup LAIlive LAIexpanded LAIdead Tair Cair VPair
#> 1     0   50  100      NA          NA      NA   NA   NA    NA
#> 2   100  150  200      NA          NA      NA   NA   NA    NA
#> 3   200  250  300      NA          NA      NA   NA   NA    NA
#> 4   300  350  400      NA          NA      NA   NA   NA    NA
#> 5   400  450  500      NA          NA      NA   NA   NA    NA
#> 6   500  550  600      NA          NA      NA   NA   NA    NA
#> 7   600  650  700      NA          NA      NA   NA   NA    NA
#> 8   700  750  800      NA          NA      NA   NA   NA    NA
#> 9   800  850  900      NA          NA      NA   NA   NA    NA
#> 10  900  950 1000      NA          NA      NA   NA   NA    NA
#> 11 1000 1050 1100      NA          NA      NA   NA   NA    NA
#> 12 1100 1150 1200      NA          NA      NA   NA   NA    NA
#> 13 1200 1250 1300      NA          NA      NA   NA   NA    NA
#> 14 1300 1350 1400      NA          NA      NA   NA   NA    NA
#> 15 1400 1450 1500      NA          NA      NA   NA   NA    NA
#> 16 1500 1550 1600      NA          NA      NA   NA   NA    NA
#> 17 1600 1650 1700      NA          NA      NA   NA   NA    NA
#> 18 1700 1750 1800      NA          NA      NA   NA   NA    NA
#> 19 1800 1850 1900      NA          NA      NA   NA   NA    NA
#> 20 1900 1950 2000      NA          NA      NA   NA   NA    NA
#> 21 2000 2050 2100      NA          NA      NA   NA   NA    NA
#> 22 2100 2150 2200      NA          NA      NA   NA   NA    NA
#> 23 2200 2250 2300      NA          NA      NA   NA   NA    NA
#> 24 2300 2350 2400      NA          NA      NA   NA   NA    NA
#> 25 2400 2450 2500      NA          NA      NA   NA   NA    NA
#> 26 2500 2550 2600      NA          NA      NA   NA   NA    NA
#> 27 2600 2650 2700      NA          NA      NA   NA   NA    NA
#> 28 2700 2750 2800      NA          NA      NA   NA   NA    NA
#> 
#> $herbLAI
#> [1] 0.1736369
#> 
#> $herbLAImax
#> [1] 0.252
#> 
#> $cohorts
#>         SP              Name
#> T1_148 148  Pinus halepensis
#> T2_168 168      Quercus ilex
#> S1_165 165 Quercus coccifera
#> 
#> $above
#>          H        CR   LAI_live LAI_expanded LAI_dead Age ObsID
#> T1_148 800 0.6605196 0.84874773   0.84874773        0  NA  <NA>
#> T2_168 660 0.6055642 0.70557382   0.70557382        0  NA  <NA>
#> S1_165  80 0.8032817 0.03062604   0.03062604        0  NA  <NA>
#> 
#> $below
#>        Z50  Z95 Z100
#> T1_148 100  600   NA
#> T2_168 300 1000   NA
#> S1_165 200 1000   NA
#> 
#> $belowLayers
#> $belowLayers$V
#>                1         2          3           4
#> T1_148 0.8604899 0.1194556 0.01511005 0.004944476
#> T2_168 0.5008953 0.4505941 0.04064831 0.007862284
#> S1_165 0.6799879 0.2737911 0.03567632 0.010544678
#> 
#> $belowLayers$L
#>               1        2        3        4
#> T1_148 2289.062 1566.552 2250.052 4226.166
#> T2_168 1817.571 2100.346 2410.127 4285.194
#> S1_165 1085.030 1380.808 2170.587 4146.637
#> 
#> $belowLayers$VGrhizo_kmax
#>               1        2         3         4
#> T1_148 11478618  1593494  201562.5  65957.49
#> T2_168 36121247 32493859 2931286.8 566975.78
#> S1_165 10941459  4405482  574055.7 169670.90
#> 
#> $belowLayers$VCroot_kmax
#>               1         2          3           4
#> T1_148 2.382795 0.4833484 0.04256689 0.007416044
#> T2_168 1.568929 1.2213562 0.09601747 0.010445417
#> S1_165 2.407779 0.7618041 0.06314806 0.009770000
#> 
#> $belowLayers$Wpool
#>        1 2 3 4
#> T1_148 1 1 1 1
#> T2_168 1 1 1 1
#> S1_165 1 1 1 1
#> 
#> $belowLayers$RhizoPsi
#>             1      2      3      4
#> T1_148 -0.033 -0.033 -0.033 -0.033
#> T2_168 -0.033 -0.033 -0.033 -0.033
#> S1_165 -0.033 -0.033 -0.033 -0.033
#> 
#> 
#> $paramsPhenology
#>             PhenologyType LeafDuration t0gdd  Sgdd Tbgdd  Ssen Phsen Tbsen xsen
#> T1_148 oneflush-evergreen     2.536875  50.0 200.0  0.00  8268  12.5  28.5    2
#> T2_168 oneflush-evergreen     2.183837  54.5 240.7  4.34 10178  12.5  28.5    2
#> S1_165 oneflush-evergreen     1.250000  54.5 240.7  4.34 10178  12.5  28.5    2
#>        ysen
#> T1_148    2
#> T2_168    2
#> S1_165    2
#> 
#> $paramsAnatomy
#>        Hmed    Al2As      SLA LeafWidth LeafDensity WoodDensity FineRootDensity
#> T1_148  850 1317.523 5.140523 0.1384772   0.2982842   0.6077016       0.2982842
#> T2_168  500 3908.823 6.340000 1.7674359   0.4893392   0.9008264       0.4893392
#> S1_165   80 4189.325 4.980084 1.3761085   0.3709679   0.4389106       0.3709679
#>        conduit2sapwood      SRL RLD     r635
#> T1_148       0.9236406 3172.572  10 1.964226
#> T2_168       0.6238125 4398.812  10 1.805872
#> S1_165       0.6238125 4398.812  10 2.289452
#> 
#> $paramsInterception
#>        LeafAngle LeafAngleSD   Beta_p   Beta_q ClumpingIndex kPAR alphaSWR
#> T1_148      53.7       21.55 1.907817 1.289641          0.75 0.50      0.7
#> T2_168      53.7       21.55 1.907817 1.289641          0.75 0.55      0.7
#> S1_165      53.7       21.55 1.907817 1.289641          0.75 0.55      0.7
#>        gammaSWR   g
#> T1_148     0.14 1.0
#> T2_168     0.18 0.5
#> S1_165     0.18 0.5
#> 
#> $paramsTranspiration
#>             Gswmin    Gswmax  Vmax298  Jmax298 Kmax_stemxylem Kmax_rootxylem
#> T1_148 0.003086667 0.2850000 72.19617 124.1687           0.15           0.60
#> T2_168 0.004473333 0.2007222 68.51600 118.7863           0.40           1.60
#> S1_165 0.010455247 0.2830167 62.78100 118.4486           0.29           1.16
#>        VCleaf_kmax VCleafapo_kmax VCleaf_slope VCleaf_P50  VCleaf_c  VCleaf_d
#> T1_148    4.000000        8.00000    133.86620  -2.303772 11.137050 -2.380849
#> T2_168    4.000000        8.00000     19.14428  -1.964085  1.339370 -2.582279
#> S1_165    9.579077       19.15815     25.47382  -2.663333  2.254991 -3.133381
#>        kleaf_symp VCstem_kmax VCstem_slope VCstem_P50  VCstem_c  VCstem_d
#> T1_148    8.00000    1.339563     68.30291  -5.139633 12.709999 -5.290000
#> T2_168    8.00000    1.620936     14.60786  -6.964747  3.560000 -7.720000
#> S1_165   19.15815    4.599269     12.31134  -6.980000  3.095442 -7.857378
#>        VCroot_kmax VCroot_slope VCroot_P50  VCroot_c  VCroot_d VGrhizo_kmax
#> T1_148    2.916127    103.96607  -2.966325 11.137050 -3.065569     13339632
#> T2_168    2.896748     22.32794  -1.684034  1.339370 -2.214081     72113368
#> S1_165    3.242501     31.37100  -1.173000  1.402489 -1.523324     16090667
#>        Plant_kmax   FR_leaf   FR_stem   FR_root
#> T1_148  0.7465846 0.1866462 0.5573346 0.2560193
#> T2_168  0.8249857 0.2062464 0.5089563 0.2847972
#> S1_165  1.5867376 0.1656462 0.3449978 0.4893561
#> 
#> $paramsWaterStorage
#>           maxFMC maxMCleaf maxMCstem   LeafPI0   LeafEPS LeafAF     Vleaf
#> T1_148 126.03063  151.9063  99.19498 -1.591429  8.918571 0.3525 0.5258525
#> T2_168  93.15304  131.4346  45.64970 -1.483333 19.260000 0.1700 0.2199087
#> S1_165  96.53441   47.3984 134.64052 -2.370000 17.230000 0.2400 0.4108968
#>          StemPI0   StemEPS    StemAF Vsapwood
#> T1_148 -2.008039 13.256355 0.9236406 6.174277
#> T2_168 -3.227438 46.420610 0.6238125 1.278142
#> S1_165 -1.305868  6.297155 0.6238125 1.064511
#> 
#> $internalPhenology
#>        gdd sen budFormation leafUnfolding leafSenescence leafDormancy phi
#> T1_148   0   0        FALSE         FALSE          FALSE        FALSE   0
#> T2_168   0   0        FALSE         FALSE          FALSE        FALSE   0
#> S1_165   0   0        FALSE         FALSE          FALSE        FALSE   0
#> 
#> $internalWater
#>        Einst RootCrownPsi LeafPsi StemPsi LeafSympPsi StemSympPsi LeafPLC
#> T1_148     0       -0.033  -0.033  -0.033      -0.033      -0.033       0
#> T2_168     0       -0.033  -0.033  -0.033      -0.033      -0.033       0
#> S1_165     0       -0.033  -0.033  -0.033      -0.033      -0.033       0
#>        StemPLC
#> T1_148       0
#> T2_168       0
#> S1_165       0
#> 
#> $internalLAIDistribution
#> $internalLAIDistribution$PrevLAIexpanded
#> [1] NA NA NA
#> 
#> $internalLAIDistribution$PrevLAIdead
#> [1] NA NA NA
#> 
#> $internalLAIDistribution$PARcohort
#> [1] 0 0 0
#> 
#> $internalLAIDistribution$live
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> $internalLAIDistribution$expanded
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> $internalLAIDistribution$dead
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> 
#> $internalFCCS
#> data frame with 0 columns and 0 rows
#> 
#> $version
#> [1] "4.8.5"
#> 
#> attr(,"class")
#> [1] "spwbInput" "list"     

# Prepare input for 'Sureau' transpiration mode
control <- defaultControl("Sureau")
spwbInput(exampleforest,examplesoil,SpParamsMED, control)
#> $control
#> $control$fillMissingRootParams
#> [1] TRUE
#> 
#> $control$fillMissingSpParams
#> [1] TRUE
#> 
#> $control$fillMissingWithGenusParams
#> [1] TRUE
#> 
#> $control$verbose
#> [1] TRUE
#> 
#> $control$subdailyResults
#> [1] FALSE
#> 
#> $control$standResults
#> [1] TRUE
#> 
#> $control$soilResults
#> [1] TRUE
#> 
#> $control$soilPoolResults
#> [1] FALSE
#> 
#> $control$snowResults
#> [1] TRUE
#> 
#> $control$plantResults
#> [1] TRUE
#> 
#> $control$labileCarbonBalanceResults
#> [1] TRUE
#> 
#> $control$plantStructureResults
#> [1] TRUE
#> 
#> $control$growthMortalityResults
#> [1] TRUE
#> 
#> $control$leafResults
#> [1] TRUE
#> 
#> $control$temperatureResults
#> [1] TRUE
#> 
#> $control$fireHazardResults
#> [1] FALSE
#> 
#> $control$fireHazardStandardWind
#> [1] NA
#> 
#> $control$fireHazardStandardDFMC
#> [1] NA
#> 
#> $control$transpirationMode
#> [1] "Sureau"
#> 
#> $control$soilDomains
#> [1] "buckets"
#> 
#> $control$rhizosphereOverlap
#> [1] "total"
#> 
#> $control$truncateRootDistribution
#> [1] FALSE
#> 
#> $control$fullRhizosphereOverlapConductivity
#> [1] 0.01
#> 
#> $control$soilFunctions
#> [1] "VG"
#> 
#> $control$VG_PTF
#> [1] "Toth"
#> 
#> $control$ndailysteps
#> [1] 24
#> 
#> $control$max_nsubsteps_soil
#> [1] 300
#> 
#> $control$defaultWindSpeed
#> [1] 2.5
#> 
#> $control$defaultCO2
#> [1] 386
#> 
#> $control$defaultRainfallIntensityPerMonth
#>  [1] 1.5 1.5 1.5 1.5 1.5 1.5 5.6 5.6 5.6 5.6 5.6 1.5
#> 
#> $control$leafPhenology
#> [1] TRUE
#> 
#> $control$bareSoilEvaporation
#> [1] TRUE
#> 
#> $control$unlimitedSoilWater
#> [1] FALSE
#> 
#> $control$interceptionMode
#> [1] "Gash1995"
#> 
#> $control$infiltrationMode
#> [1] "GreenAmpt1911"
#> 
#> $control$infiltrationCorrection
#> [1] 5
#> 
#> $control$unfoldingDD
#> [1] 300
#> 
#> $control$verticalLayerSize
#> [1] 100
#> 
#> $control$windMeasurementHeight
#> [1] 200
#> 
#> $control$segmentedXylemVulnerability
#> [1] FALSE
#> 
#> $control$stemCavitationRecovery
#> [1] "annual"
#> 
#> $control$leafCavitationRecovery
#> [1] "annual"
#> 
#> $control$lfmcComponent
#> [1] "fine"
#> 
#> $control$hydraulicRedistributionFraction
#> [1] 0.1
#> 
#> $control$nsubsteps_canopy
#> [1] 3600
#> 
#> $control$taper
#> [1] TRUE
#> 
#> $control$multiLayerBalance
#> [1] FALSE
#> 
#> $control$sapFluidityVariation
#> [1] TRUE
#> 
#> $control$TPhase_gmin
#> [1] 37.5
#> 
#> $control$Q10_1_gmin
#> [1] 1.2
#> 
#> $control$Q10_2_gmin
#> [1] 4.8
#> 
#> $control$rootRadialConductance
#> [1] 4
#> 
#> $control$averageFracRhizosphereResistance
#> [1] 0.15
#> 
#> $control$thermalCapacityLAI
#> [1] 1e+06
#> 
#> $control$boundaryLayerSize
#> [1] 2000
#> 
#> $control$cavitationRecoveryMaximumRate
#> [1] 0.05
#> 
#> $control$sunlitShade
#> [1] TRUE
#> 
#> $control$numericParams
#> $control$numericParams$maxNsteps
#> [1] 400
#> 
#> $control$numericParams$ntrial
#> [1] 200
#> 
#> $control$numericParams$psiTol
#> [1] 1e-04
#> 
#> $control$numericParams$ETol
#> [1] 1e-07
#> 
#> 
#> $control$leafCavitationEffects
#> [1] FALSE
#> 
#> $control$stemCavitationEffects
#> [1] TRUE
#> 
#> $control$stomatalSubmodel
#> [1] "Baldocchi"
#> 
#> $control$plantCapacitance
#> [1] TRUE
#> 
#> $control$cavitationFlux
#> [1] TRUE
#> 
#> $control$leafCuticularTranspiration
#> [1] TRUE
#> 
#> $control$stemCuticularTranspiration
#> [1] FALSE
#> 
#> $control$C_SApoInit
#> [1] 2e-05
#> 
#> $control$C_LApoInit
#> [1] 1e-05
#> 
#> $control$k_SSym
#> [1] 0.26
#> 
#> $control$fractionLeafSymplasm
#> [1] 0.5
#> 
#> $control$gs_NightFrac
#> [1] 0.05
#> 
#> $control$JarvisPAR
#> [1] 0.003
#> 
#> $control$fTRBToLeaf
#> [1] 0.8
#> 
#> $control$subdailyCarbonBalance
#> [1] FALSE
#> 
#> $control$allowDessication
#> [1] TRUE
#> 
#> $control$allowStarvation
#> [1] TRUE
#> 
#> $control$sinkLimitation
#> [1] TRUE
#> 
#> $control$shrubDynamics
#> [1] TRUE
#> 
#> $control$herbDynamics
#> [1] TRUE
#> 
#> $control$allocationStrategy
#> [1] "Al2As"
#> 
#> $control$phloemConductanceFactor
#> [1] 0.2
#> 
#> $control$nonSugarConcentration
#> [1] 0.25
#> 
#> $control$equilibriumOsmoticConcentration
#> $control$equilibriumOsmoticConcentration$leaf
#> [1] 0.8
#> 
#> $control$equilibriumOsmoticConcentration$sapwood
#> [1] 0.6
#> 
#> 
#> $control$minimumRelativeStarchForGrowth
#> [1] 0.5
#> 
#> $control$constructionCosts
#> $control$constructionCosts$leaf
#> [1] 1.5
#> 
#> $control$constructionCosts$sapwood
#> [1] 1.47
#> 
#> $control$constructionCosts$fineroot
#> [1] 1.3
#> 
#> 
#> $control$senescenceRates
#> $control$senescenceRates$sapwood
#> [1] 0.000135
#> 
#> $control$senescenceRates$fineroot
#> [1] 0.001897231
#> 
#> 
#> $control$maximumRelativeGrowthRates
#> $control$maximumRelativeGrowthRates$leaf
#> [1] 0.09
#> 
#> $control$maximumRelativeGrowthRates$cambium
#> [1] 0.0025
#> 
#> $control$maximumRelativeGrowthRates$sapwood
#> [1] 0.002
#> 
#> $control$maximumRelativeGrowthRates$fineroot
#> [1] 0.1
#> 
#> 
#> $control$mortalityMode
#> [1] "density/deterministic"
#> 
#> $control$mortalityBaselineRate
#> [1] 0.0015
#> 
#> $control$mortalityRelativeSugarThreshold
#> [1] 0.4
#> 
#> $control$mortalityRWCThreshold
#> [1] 0.4
#> 
#> $control$recrTreeDBH
#> [1] 1
#> 
#> $control$recrTreeDensity
#> [1] 3000
#> 
#> $control$ingrowthTreeDBH
#> [1] 7.5
#> 
#> $control$ingrowthTreeDensity
#> [1] 127
#> 
#> $control$allowSeedBankDynamics
#> [1] TRUE
#> 
#> $control$allowRecruitment
#> [1] TRUE
#> 
#> $control$allowResprouting
#> [1] TRUE
#> 
#> $control$recruitmentMode
#> [1] "stochastic"
#> 
#> $control$removeEmptyCohorts
#> [1] TRUE
#> 
#> $control$minimumTreeCohortDensity
#> [1] 1
#> 
#> $control$minimumShrubCohortCover
#> [1] 0.01
#> 
#> $control$dynamicallyMergeCohorts
#> [1] TRUE
#> 
#> $control$keepCohortsWithObsID
#> [1] FALSE
#> 
#> $control$seedRain
#> NULL
#> 
#> $control$seedProductionTreeHeight
#> [1] 300
#> 
#> $control$seedProductionShrubHeight
#> [1] 30
#> 
#> $control$probRecr
#> [1] 0.05
#> 
#> $control$minTempRecr
#> [1] 0
#> 
#> $control$minMoistureRecr
#> [1] 0.3
#> 
#> $control$minFPARRecr
#> [1] 10
#> 
#> $control$recrTreeHeight
#> [1] 620
#> 
#> $control$recrShrubCover
#> [1] 1
#> 
#> $control$recrShrubHeight
#> [1] 25
#> 
#> $control$recrTreeZ50
#> [1] 100
#> 
#> $control$recrShrubZ50
#> [1] 50
#> 
#> $control$recrTreeZ95
#> [1] 1000
#> 
#> $control$recrShrubZ95
#> [1] 500
#> 
#> 
#> $soil
#>   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha
#> 1    300   25   25 Silt loam NA       NA 1.5  25 0.0485 5401.471 89.16112
#> 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112
#> 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112
#> 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112
#>       VG_n VG_theta_res VG_theta_sat W Temp
#> 1 1.303861        0.041     0.423715 1   NA
#> 2 1.303861        0.041     0.423715 1   NA
#> 3 1.303861        0.041     0.423715 1   NA
#> 4 1.303861        0.041     0.423715 1   NA
#> 
#> $snowpack
#> [1] 0
#> 
#> $canopy
#>    zlow zmid  zup LAIlive LAIexpanded LAIdead Tair Cair VPair
#> 1     0   50  100      NA          NA      NA   NA   NA    NA
#> 2   100  150  200      NA          NA      NA   NA   NA    NA
#> 3   200  250  300      NA          NA      NA   NA   NA    NA
#> 4   300  350  400      NA          NA      NA   NA   NA    NA
#> 5   400  450  500      NA          NA      NA   NA   NA    NA
#> 6   500  550  600      NA          NA      NA   NA   NA    NA
#> 7   600  650  700      NA          NA      NA   NA   NA    NA
#> 8   700  750  800      NA          NA      NA   NA   NA    NA
#> 9   800  850  900      NA          NA      NA   NA   NA    NA
#> 10  900  950 1000      NA          NA      NA   NA   NA    NA
#> 11 1000 1050 1100      NA          NA      NA   NA   NA    NA
#> 12 1100 1150 1200      NA          NA      NA   NA   NA    NA
#> 13 1200 1250 1300      NA          NA      NA   NA   NA    NA
#> 14 1300 1350 1400      NA          NA      NA   NA   NA    NA
#> 15 1400 1450 1500      NA          NA      NA   NA   NA    NA
#> 16 1500 1550 1600      NA          NA      NA   NA   NA    NA
#> 17 1600 1650 1700      NA          NA      NA   NA   NA    NA
#> 18 1700 1750 1800      NA          NA      NA   NA   NA    NA
#> 19 1800 1850 1900      NA          NA      NA   NA   NA    NA
#> 20 1900 1950 2000      NA          NA      NA   NA   NA    NA
#> 21 2000 2050 2100      NA          NA      NA   NA   NA    NA
#> 22 2100 2150 2200      NA          NA      NA   NA   NA    NA
#> 23 2200 2250 2300      NA          NA      NA   NA   NA    NA
#> 24 2300 2350 2400      NA          NA      NA   NA   NA    NA
#> 25 2400 2450 2500      NA          NA      NA   NA   NA    NA
#> 26 2500 2550 2600      NA          NA      NA   NA   NA    NA
#> 27 2600 2650 2700      NA          NA      NA   NA   NA    NA
#> 28 2700 2750 2800      NA          NA      NA   NA   NA    NA
#> 
#> $herbLAI
#> [1] 0.1736369
#> 
#> $herbLAImax
#> [1] 0.252
#> 
#> $cohorts
#>         SP              Name
#> T1_148 148  Pinus halepensis
#> T2_168 168      Quercus ilex
#> S1_165 165 Quercus coccifera
#> 
#> $above
#>          H        CR   LAI_live LAI_expanded LAI_dead Age ObsID
#> T1_148 800 0.6605196 0.84874773   0.84874773        0  NA  <NA>
#> T2_168 660 0.6055642 0.70557382   0.70557382        0  NA  <NA>
#> S1_165  80 0.8032817 0.03062604   0.03062604        0  NA  <NA>
#> 
#> $below
#>        Z50  Z95 Z100
#> T1_148 100  600   NA
#> T2_168 300 1000   NA
#> S1_165 200 1000   NA
#> 
#> $belowLayers
#> $belowLayers$V
#>                1         2          3           4
#> T1_148 0.8604899 0.1194556 0.01511005 0.004944476
#> T2_168 0.5008953 0.4505941 0.04064831 0.007862284
#> S1_165 0.6799879 0.2737911 0.03567632 0.010544678
#> 
#> $belowLayers$L
#>               1        2        3        4
#> T1_148 2289.062 1566.552 2250.052 4226.166
#> T2_168 1817.571 2100.346 2410.127 4285.194
#> S1_165 1085.030 1380.808 2170.587 4146.637
#> 
#> $belowLayers$VGrhizo_kmax
#>                 1         2        3        4
#> T1_148  278101462  38606817  4883412  1598004
#> T2_168  406976825 366107176 33026706  6388096
#> S1_165 1116059122 449371399 58555276 17306902
#> 
#> $belowLayers$VCroot_kmax
#>               1         2          3           4
#> T1_148 2.382795 0.4833484 0.04256689 0.007416044
#> T2_168 1.568929 1.2213562 0.09601747 0.010445417
#> S1_165 2.407779 0.7618041 0.06314806 0.009770000
#> 
#> $belowLayers$Wpool
#>        1 2 3 4
#> T1_148 1 1 1 1
#> T2_168 1 1 1 1
#> S1_165 1 1 1 1
#> 
#> $belowLayers$RhizoPsi
#>             1      2      3      4
#> T1_148 -0.033 -0.033 -0.033 -0.033
#> T2_168 -0.033 -0.033 -0.033 -0.033
#> S1_165 -0.033 -0.033 -0.033 -0.033
#> 
#> 
#> $paramsPhenology
#>             PhenologyType LeafDuration t0gdd  Sgdd Tbgdd  Ssen Phsen Tbsen xsen
#> T1_148 oneflush-evergreen     2.536875  50.0 200.0  0.00  8268  12.5  28.5    2
#> T2_168 oneflush-evergreen     2.183837  54.5 240.7  4.34 10178  12.5  28.5    2
#> S1_165 oneflush-evergreen     1.250000  54.5 240.7  4.34 10178  12.5  28.5    2
#>        ysen
#> T1_148    2
#> T2_168    2
#> S1_165    2
#> 
#> $paramsAnatomy
#>        Hmed    Al2As      SLA LeafWidth LeafDensity WoodDensity FineRootDensity
#> T1_148  850 1317.523 5.140523 0.1384772   0.2982842   0.6077016       0.2982842
#> T2_168  500 3908.823 6.340000 1.7674359   0.4893392   0.9008264       0.4893392
#> S1_165   80 4189.325 4.980084 1.3761085   0.3709679   0.4389106       0.3709679
#>        conduit2sapwood      SRL RLD     r635
#> T1_148       0.9236406 3172.572  10 1.964226
#> T2_168       0.6238125 4398.812  10 1.805872
#> S1_165       0.6238125 4398.812  10 2.289452
#> 
#> $paramsInterception
#>        LeafAngle LeafAngleSD   Beta_p   Beta_q ClumpingIndex kPAR alphaSWR
#> T1_148      53.7       21.55 1.907817 1.289641          0.75 0.50      0.7
#> T2_168      53.7       21.55 1.907817 1.289641          0.75 0.55      0.7
#> S1_165      53.7       21.55 1.907817 1.289641          0.75 0.55      0.7
#>        gammaSWR   g
#> T1_148     0.14 1.0
#> T2_168     0.18 0.5
#> S1_165     0.18 0.5
#> 
#> $paramsTranspiration
#>             Gswmin    Gswmax Gsw_AC_slope    Gs_P50 Gs_slope  Vmax298  Jmax298
#> T1_148 0.003086667 0.2850000     6.238912 -2.303772       30 72.19617 124.1687
#> T2_168 0.004473333 0.2007222     4.957957 -1.964085       30 68.51600 118.7863
#> S1_165 0.010455247 0.2830167     6.590920 -2.663333       30 62.78100 118.4486
#>        Kmax_stemxylem Kmax_rootxylem VCleaf_kmax VCleafapo_kmax VCleaf_slope
#> T1_148           0.15           0.60    4.000000        8.00000     68.30291
#> T2_168           0.40           1.60    4.000000        8.00000     14.60786
#> S1_165           0.29           1.16    9.579077       19.15815     12.52202
#>        VCleaf_P50  VCleaf_c  VCleaf_d kleaf_symp VCstem_kmax VCstem_slope
#> T1_148  -5.139633 12.709999 -5.290000    8.00000    1.339563     68.30291
#> T2_168  -6.964747  3.560000 -7.720000    8.00000    1.620936     14.60786
#> S1_165  -6.980000  3.095442 -7.857378   19.15815    4.599269     12.52202
#>        VCstem_P50  VCstem_c  VCstem_d kstem_symp VCroot_kmax VCroot_slope
#> T1_148  -5.139633 12.709999 -5.290000       0.26    2.916127     68.30291
#> T2_168  -6.964747  3.560000 -7.720000       0.26    2.896748     14.60786
#> S1_165  -6.980000  3.095442 -7.857378       0.26    3.242501     12.52202
#>        VCroot_P50  VCroot_c  VCroot_d VGrhizo_kmax Plant_kmax   FR_leaf
#> T1_148  -5.139633 12.709999 -5.290000    323189695  0.7465846 0.1866462
#> T2_168  -6.964747  3.560000 -7.720000    812498804  0.8249857 0.2062464
#> S1_165  -6.980000  3.095442 -7.857378   1641292700  1.5867376 0.1656462
#>          FR_stem   FR_root
#> T1_148 0.5573346 0.2560193
#> T2_168 0.5089563 0.2847972
#> S1_165 0.3449978 0.4893561
#> 
#> $paramsWaterStorage
#>           maxFMC maxMCleaf maxMCstem   LeafPI0   LeafEPS LeafAF     Vleaf
#> T1_148 126.03063  151.9063  99.19498 -1.591429  8.918571 0.3525 0.5258525
#> T2_168  93.15304  131.4346  45.64970 -1.483333 19.260000 0.1700 0.2199087
#> S1_165  96.53441   47.3984 134.64052 -2.370000 17.230000 0.2400 0.4108968
#>          StemPI0   StemEPS    StemAF Vsapwood
#> T1_148 -2.008039 13.256355 0.9236406 6.174277
#> T2_168 -3.227438 46.420610 0.6238125 1.278142
#> S1_165 -1.305868  6.297155 0.6238125 1.064511
#> 
#> $internalPhenology
#>        gdd sen budFormation leafUnfolding leafSenescence leafDormancy phi
#> T1_148   0   0        FALSE         FALSE          FALSE        FALSE   0
#> T2_168   0   0        FALSE         FALSE          FALSE        FALSE   0
#> S1_165   0   0        FALSE         FALSE          FALSE        FALSE   0
#> 
#> $internalWater
#>        Einst Elim Emin_L Emin_S RootCrownPsi LeafPsi StemPsi LeafSympPsi
#> T1_148     0    0      0      0       -0.033  -0.033  -0.033      -0.033
#> T2_168     0    0      0      0       -0.033  -0.033  -0.033      -0.033
#> S1_165     0    0      0      0       -0.033  -0.033  -0.033      -0.033
#>        StemSympPsi LeafPLC StemPLC
#> T1_148      -0.033       0       0
#> T2_168      -0.033       0       0
#> S1_165      -0.033       0       0
#> 
#> $internalLAIDistribution
#> $internalLAIDistribution$PrevLAIexpanded
#> [1] NA NA NA
#> 
#> $internalLAIDistribution$PrevLAIdead
#> [1] NA NA NA
#> 
#> $internalLAIDistribution$PARcohort
#> [1] 0 0 0
#> 
#> $internalLAIDistribution$live
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> $internalLAIDistribution$expanded
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> $internalLAIDistribution$dead
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> 
#> $internalFCCS
#> data frame with 0 columns and 0 rows
#> 
#> $version
#> [1] "4.8.5"
#> 
#> attr(,"class")
#> [1] "spwbInput" "list"     

# Example of initialization from a forest 
# described using LAI and crown ratio
control <- defaultControl("Granier")
spwbInput(exampleforest2, examplesoil, SpParamsMED, control)
#> $control
#> $control$fillMissingRootParams
#> [1] TRUE
#> 
#> $control$fillMissingSpParams
#> [1] TRUE
#> 
#> $control$fillMissingWithGenusParams
#> [1] TRUE
#> 
#> $control$verbose
#> [1] TRUE
#> 
#> $control$subdailyResults
#> [1] FALSE
#> 
#> $control$standResults
#> [1] TRUE
#> 
#> $control$soilResults
#> [1] TRUE
#> 
#> $control$soilPoolResults
#> [1] FALSE
#> 
#> $control$snowResults
#> [1] TRUE
#> 
#> $control$plantResults
#> [1] TRUE
#> 
#> $control$labileCarbonBalanceResults
#> [1] TRUE
#> 
#> $control$plantStructureResults
#> [1] TRUE
#> 
#> $control$growthMortalityResults
#> [1] TRUE
#> 
#> $control$leafResults
#> [1] TRUE
#> 
#> $control$temperatureResults
#> [1] TRUE
#> 
#> $control$fireHazardResults
#> [1] FALSE
#> 
#> $control$fireHazardStandardWind
#> [1] NA
#> 
#> $control$fireHazardStandardDFMC
#> [1] NA
#> 
#> $control$transpirationMode
#> [1] "Granier"
#> 
#> $control$soilDomains
#> [1] "buckets"
#> 
#> $control$rhizosphereOverlap
#> [1] "total"
#> 
#> $control$truncateRootDistribution
#> [1] FALSE
#> 
#> $control$fullRhizosphereOverlapConductivity
#> [1] 0.01
#> 
#> $control$soilFunctions
#> [1] "VG"
#> 
#> $control$VG_PTF
#> [1] "Toth"
#> 
#> $control$ndailysteps
#> [1] 24
#> 
#> $control$max_nsubsteps_soil
#> [1] 300
#> 
#> $control$defaultWindSpeed
#> [1] 2.5
#> 
#> $control$defaultCO2
#> [1] 386
#> 
#> $control$defaultRainfallIntensityPerMonth
#>  [1] 1.5 1.5 1.5 1.5 1.5 1.5 5.6 5.6 5.6 5.6 5.6 1.5
#> 
#> $control$leafPhenology
#> [1] TRUE
#> 
#> $control$bareSoilEvaporation
#> [1] TRUE
#> 
#> $control$unlimitedSoilWater
#> [1] FALSE
#> 
#> $control$interceptionMode
#> [1] "Gash1995"
#> 
#> $control$infiltrationMode
#> [1] "GreenAmpt1911"
#> 
#> $control$infiltrationCorrection
#> [1] 5
#> 
#> $control$unfoldingDD
#> [1] 300
#> 
#> $control$verticalLayerSize
#> [1] 100
#> 
#> $control$windMeasurementHeight
#> [1] 200
#> 
#> $control$segmentedXylemVulnerability
#> [1] TRUE
#> 
#> $control$stemCavitationRecovery
#> [1] "annual"
#> 
#> $control$leafCavitationRecovery
#> [1] "annual"
#> 
#> $control$lfmcComponent
#> [1] "fine"
#> 
#> $control$hydraulicRedistributionFraction
#> [1] 0.1
#> 
#> $control$nsubsteps_canopy
#> [1] 3600
#> 
#> $control$taper
#> [1] TRUE
#> 
#> $control$multiLayerBalance
#> [1] FALSE
#> 
#> $control$sapFluidityVariation
#> [1] TRUE
#> 
#> $control$TPhase_gmin
#> [1] 37.5
#> 
#> $control$Q10_1_gmin
#> [1] 1.2
#> 
#> $control$Q10_2_gmin
#> [1] 4.8
#> 
#> $control$rootRadialConductance
#> [1] 4
#> 
#> $control$averageFracRhizosphereResistance
#> [1] 0.15
#> 
#> $control$thermalCapacityLAI
#> [1] 1e+06
#> 
#> $control$boundaryLayerSize
#> [1] 2000
#> 
#> $control$cavitationRecoveryMaximumRate
#> [1] 0.05
#> 
#> $control$sunlitShade
#> [1] TRUE
#> 
#> $control$numericParams
#> $control$numericParams$maxNsteps
#> [1] 400
#> 
#> $control$numericParams$ntrial
#> [1] 200
#> 
#> $control$numericParams$psiTol
#> [1] 1e-04
#> 
#> $control$numericParams$ETol
#> [1] 1e-07
#> 
#> 
#> $control$leafCavitationEffects
#> [1] FALSE
#> 
#> $control$stemCavitationEffects
#> [1] TRUE
#> 
#> $control$stomatalSubmodel
#> [1] "Baldocchi"
#> 
#> $control$plantCapacitance
#> [1] TRUE
#> 
#> $control$cavitationFlux
#> [1] TRUE
#> 
#> $control$leafCuticularTranspiration
#> [1] TRUE
#> 
#> $control$stemCuticularTranspiration
#> [1] FALSE
#> 
#> $control$C_SApoInit
#> [1] 2e-05
#> 
#> $control$C_LApoInit
#> [1] 1e-05
#> 
#> $control$k_SSym
#> [1] 0.26
#> 
#> $control$fractionLeafSymplasm
#> [1] 0.5
#> 
#> $control$gs_NightFrac
#> [1] 0.05
#> 
#> $control$JarvisPAR
#> [1] 0.003
#> 
#> $control$fTRBToLeaf
#> [1] 0.8
#> 
#> $control$subdailyCarbonBalance
#> [1] FALSE
#> 
#> $control$allowDessication
#> [1] TRUE
#> 
#> $control$allowStarvation
#> [1] TRUE
#> 
#> $control$sinkLimitation
#> [1] TRUE
#> 
#> $control$shrubDynamics
#> [1] TRUE
#> 
#> $control$herbDynamics
#> [1] TRUE
#> 
#> $control$allocationStrategy
#> [1] "Al2As"
#> 
#> $control$phloemConductanceFactor
#> [1] 0.2
#> 
#> $control$nonSugarConcentration
#> [1] 0.25
#> 
#> $control$equilibriumOsmoticConcentration
#> $control$equilibriumOsmoticConcentration$leaf
#> [1] 0.8
#> 
#> $control$equilibriumOsmoticConcentration$sapwood
#> [1] 0.6
#> 
#> 
#> $control$minimumRelativeStarchForGrowth
#> [1] 0.5
#> 
#> $control$constructionCosts
#> $control$constructionCosts$leaf
#> [1] 1.5
#> 
#> $control$constructionCosts$sapwood
#> [1] 1.47
#> 
#> $control$constructionCosts$fineroot
#> [1] 1.3
#> 
#> 
#> $control$senescenceRates
#> $control$senescenceRates$sapwood
#> [1] 0.000135
#> 
#> $control$senescenceRates$fineroot
#> [1] 0.001897231
#> 
#> 
#> $control$maximumRelativeGrowthRates
#> $control$maximumRelativeGrowthRates$leaf
#> [1] 0.09
#> 
#> $control$maximumRelativeGrowthRates$cambium
#> [1] 0.0025
#> 
#> $control$maximumRelativeGrowthRates$sapwood
#> [1] 0.002
#> 
#> $control$maximumRelativeGrowthRates$fineroot
#> [1] 0.1
#> 
#> 
#> $control$mortalityMode
#> [1] "density/deterministic"
#> 
#> $control$mortalityBaselineRate
#> [1] 0.0015
#> 
#> $control$mortalityRelativeSugarThreshold
#> [1] 0.4
#> 
#> $control$mortalityRWCThreshold
#> [1] 0.4
#> 
#> $control$recrTreeDBH
#> [1] 1
#> 
#> $control$recrTreeDensity
#> [1] 3000
#> 
#> $control$ingrowthTreeDBH
#> [1] 7.5
#> 
#> $control$ingrowthTreeDensity
#> [1] 127
#> 
#> $control$allowSeedBankDynamics
#> [1] TRUE
#> 
#> $control$allowRecruitment
#> [1] TRUE
#> 
#> $control$allowResprouting
#> [1] TRUE
#> 
#> $control$recruitmentMode
#> [1] "stochastic"
#> 
#> $control$removeEmptyCohorts
#> [1] TRUE
#> 
#> $control$minimumTreeCohortDensity
#> [1] 1
#> 
#> $control$minimumShrubCohortCover
#> [1] 0.01
#> 
#> $control$dynamicallyMergeCohorts
#> [1] TRUE
#> 
#> $control$keepCohortsWithObsID
#> [1] FALSE
#> 
#> $control$seedRain
#> NULL
#> 
#> $control$seedProductionTreeHeight
#> [1] 300
#> 
#> $control$seedProductionShrubHeight
#> [1] 30
#> 
#> $control$probRecr
#> [1] 0.05
#> 
#> $control$minTempRecr
#> [1] 0
#> 
#> $control$minMoistureRecr
#> [1] 0.3
#> 
#> $control$minFPARRecr
#> [1] 10
#> 
#> $control$recrTreeHeight
#> [1] 620
#> 
#> $control$recrShrubCover
#> [1] 1
#> 
#> $control$recrShrubHeight
#> [1] 25
#> 
#> $control$recrTreeZ50
#> [1] 100
#> 
#> $control$recrShrubZ50
#> [1] 50
#> 
#> $control$recrTreeZ95
#> [1] 1000
#> 
#> $control$recrShrubZ95
#> [1] 500
#> 
#> 
#> $soil
#>   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha
#> 1    300   25   25 Silt loam NA       NA 1.5  25 0.0485 5401.471 89.16112
#> 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112
#> 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112
#> 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112
#>       VG_n VG_theta_res VG_theta_sat W Temp
#> 1 1.303861        0.041     0.423715 1   NA
#> 2 1.303861        0.041     0.423715 1   NA
#> 3 1.303861        0.041     0.423715 1   NA
#> 4 1.303861        0.041     0.423715 1   NA
#> 
#> $snowpack
#> [1] 0
#> 
#> $canopy
#>    zlow zmid  zup LAIlive LAIexpanded LAIdead Tair Cair VPair
#> 1     0   50  100      NA          NA      NA   NA   NA    NA
#> 2   100  150  200      NA          NA      NA   NA   NA    NA
#> 3   200  250  300      NA          NA      NA   NA   NA    NA
#> 4   300  350  400      NA          NA      NA   NA   NA    NA
#> 5   400  450  500      NA          NA      NA   NA   NA    NA
#> 6   500  550  600      NA          NA      NA   NA   NA    NA
#> 7   600  650  700      NA          NA      NA   NA   NA    NA
#> 8   700  750  800      NA          NA      NA   NA   NA    NA
#> 9   800  850  900      NA          NA      NA   NA   NA    NA
#> 10  900  950 1000      NA          NA      NA   NA   NA    NA
#> 11 1000 1050 1100      NA          NA      NA   NA   NA    NA
#> 12 1100 1150 1200      NA          NA      NA   NA   NA    NA
#> 13 1200 1250 1300      NA          NA      NA   NA   NA    NA
#> 14 1300 1350 1400      NA          NA      NA   NA   NA    NA
#> 15 1400 1450 1500      NA          NA      NA   NA   NA    NA
#> 16 1500 1550 1600      NA          NA      NA   NA   NA    NA
#> 17 1600 1650 1700      NA          NA      NA   NA   NA    NA
#> 18 1700 1750 1800      NA          NA      NA   NA   NA    NA
#> 19 1800 1850 1900      NA          NA      NA   NA   NA    NA
#> 20 1900 1950 2000      NA          NA      NA   NA   NA    NA
#> 21 2000 2050 2100      NA          NA      NA   NA   NA    NA
#> 22 2100 2150 2200      NA          NA      NA   NA   NA    NA
#> 23 2200 2250 2300      NA          NA      NA   NA   NA    NA
#> 24 2300 2350 2400      NA          NA      NA   NA   NA    NA
#> 25 2400 2450 2500      NA          NA      NA   NA   NA    NA
#> 26 2500 2550 2600      NA          NA      NA   NA   NA    NA
#> 27 2600 2650 2700      NA          NA      NA   NA   NA    NA
#> 28 2700 2750 2800      NA          NA      NA   NA   NA    NA
#> 
#> $herbLAI
#> [1] 0
#> 
#> $herbLAImax
#> [1] 0
#> 
#> $cohorts
#>         SP              Name
#> T1_148 148  Pinus halepensis
#> T2_168 168      Quercus ilex
#> S1_165 165 Quercus coccifera
#> 
#> $above
#>          H   CR LAI_live LAI_expanded LAI_dead Age ObsID
#> T1_148 800 0.66     0.80         0.80        0  NA  <NA>
#> T2_168 660 0.60     0.50         0.50        0  NA  <NA>
#> S1_165  80 0.80     0.03         0.03        0  NA  <NA>
#> 
#> $below
#>        Z50  Z95 Z100
#> T1_148 100  600   NA
#> T2_168 300 1000   NA
#> S1_165 200 1000   NA
#> 
#> $belowLayers
#> $belowLayers$V
#>                1         2          3           4
#> T1_148 0.8604899 0.1194556 0.01511005 0.004944476
#> T2_168 0.5008953 0.4505941 0.04064831 0.007862284
#> S1_165 0.6799879 0.2737911 0.03567632 0.010544678
#> 
#> $belowLayers$L
#>          1   2    3    4
#> T1_148 150 650 1500 3000
#> T2_168 150 650 1500 3000
#> S1_165 150 650 1500 3000
#> 
#> $belowLayers$Wpool
#>        1 2 3 4
#> T1_148 1 1 1 1
#> T2_168 1 1 1 1
#> S1_165 1 1 1 1
#> 
#> 
#> $paramsPhenology
#>             PhenologyType LeafDuration t0gdd  Sgdd Tbgdd  Ssen Phsen Tbsen xsen
#> T1_148 oneflush-evergreen     2.536875  50.0 200.0  0.00  8268  12.5  28.5    2
#> T2_168 oneflush-evergreen     2.183837  54.5 240.7  4.34 10178  12.5  28.5    2
#> S1_165 oneflush-evergreen     1.250000  54.5 240.7  4.34 10178  12.5  28.5    2
#>        ysen
#> T1_148    2
#> T2_168    2
#> S1_165    2
#> 
#> $paramsAnatomy
#>           Al2As Ar2Al      SLA LeafDensity WoodDensity FineRootDensity      SRL
#> T1_148 1317.523     1 5.140523   0.2982842   0.6077016       0.2982842 3172.572
#> T2_168 3908.823     1 6.340000   0.4893392   0.9008264       0.4893392 4398.812
#> S1_165 4189.325     1 4.980084   0.3709679   0.4389106       0.3709679 4398.812
#>        RLD     r635
#> T1_148  10 1.964226
#> T2_168  10 1.805872
#> S1_165  10 2.289452
#> 
#> $paramsInterception
#>        kPAR      kSWR   g
#> T1_148 0.50 0.3703704 1.0
#> T2_168 0.55 0.4074074 0.5
#> S1_165 0.55 0.4074074 0.5
#> 
#> $paramsTranspiration
#>             Gswmin  Tmax_LAI   Tmax_LAIsq Psi_Extract Exp_Extract  VCleaf_c
#> T1_148 0.003086667 0.1869849 -0.008372458  -0.9218219    1.504542 11.137050
#> T2_168 0.004473333 0.1251027 -0.005601615  -1.9726871    1.149052  1.339370
#> S1_165 0.010455247 0.1340000 -0.006000000  -2.1210726    1.300000  2.254991
#>         VCleaf_d  VCstem_c  VCstem_d      WUE   WUE_par     WUE_co2    WUE_vpd
#> T1_148 -2.380849 12.709999 -5.290000 8.525550 0.5239136 0.002586327 -0.2647169
#> T2_168 -2.582279  3.560000 -7.720000 8.968208 0.1412266 0.002413091 -0.5664879
#> S1_165 -3.133381  3.095442 -7.857378 7.900000 0.3643000 0.002757000 -0.4636000
#> 
#> $paramsWaterStorage
#>           maxFMC maxMCleaf maxMCstem   LeafPI0   LeafEPS LeafAF     Vleaf
#> T1_148 126.03063  151.9063  99.19498 -1.591429  8.918571 0.3525 0.5258525
#> T2_168  93.15304  131.4346  45.64970 -1.483333 19.260000 0.1700 0.2199087
#> S1_165  96.53441   47.3984 134.64052 -2.370000 17.230000 0.2400 0.4108968
#>          StemPI0   StemEPS    StemAF  Vsapwood
#> T1_148 -2.008039 13.256355 0.9236406 4.1638559
#> T2_168 -3.227438 46.420610 0.6238125 0.8135590
#> S1_165 -1.305868  6.297155 0.6238125 0.3177724
#> 
#> $internalPhenology
#>        gdd sen budFormation leafUnfolding leafSenescence leafDormancy phi
#> T1_148   0   0        FALSE         FALSE          FALSE        FALSE   0
#> T2_168   0   0        FALSE         FALSE          FALSE        FALSE   0
#> S1_165   0   0        FALSE         FALSE          FALSE        FALSE   0
#> 
#> $internalWater
#>        PlantPsi LeafPLC StemPLC
#> T1_148   -0.033       0       0
#> T2_168   -0.033       0       0
#> S1_165   -0.033       0       0
#> 
#> $internalLAIDistribution
#> $internalLAIDistribution$PrevLAIexpanded
#> [1] NA NA NA
#> 
#> $internalLAIDistribution$PrevLAIdead
#> [1] NA NA NA
#> 
#> $internalLAIDistribution$PARcohort
#> [1] 0 0 0
#> 
#> $internalLAIDistribution$live
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> $internalLAIDistribution$expanded
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> $internalLAIDistribution$dead
#>    T1_148 T2_168 S1_165
#> 1       0      0      0
#> 2       0      0      0
#> 3       0      0      0
#> 4       0      0      0
#> 5       0      0      0
#> 6       0      0      0
#> 7       0      0      0
#> 8       0      0      0
#> 9       0      0      0
#> 10      0      0      0
#> 11      0      0      0
#> 12      0      0      0
#> 13      0      0      0
#> 14      0      0      0
#> 15      0      0      0
#> 16      0      0      0
#> 17      0      0      0
#> 18      0      0      0
#> 19      0      0      0
#> 20      0      0      0
#> 21      0      0      0
#> 22      0      0      0
#> 23      0      0      0
#> 24      0      0      0
#> 25      0      0      0
#> 26      0      0      0
#> 27      0      0      0
#> 28      0      0      0
#> 
#> 
#> $internalFCCS
#> data frame with 0 columns and 0 rows
#> 
#> $version
#> [1] "4.8.5"
#> 
#> attr(,"class")
#> [1] "spwbInput" "list"     
```
