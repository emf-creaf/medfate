# Transpiration modes

High-level sub-models representing transpiration, plant hydraulics,
photosynthesis and water relations within plants.

## Usage

``` r
transp_transpirationSperry(
  x,
  meteo,
  day,
  latitude,
  elevation,
  slope,
  aspect,
  canopyEvaporation = 0,
  snowMelt = 0,
  soilEvaporation = 0,
  herbTranspiration = 0,
  stepFunctions = NA_integer_,
  modifyInput = TRUE
)

transp_transpirationSureau(
  x,
  meteo,
  day,
  latitude,
  elevation,
  slope,
  aspect,
  canopyEvaporation = 0,
  snowMelt = 0,
  soilEvaporation = 0,
  herbTranspiration = 0,
  modifyInput = TRUE
)

transp_transpirationGranier(
  x,
  meteo,
  day,
  latitude,
  elevation,
  slope,
  aspect,
  modifyInput = TRUE
)
```

## Arguments

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  or
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
  built using the 'Granier', 'Sperry' or 'Sureau' transpiration modes.

- meteo:

  A data frame with daily meteorological data series (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- day:

  An integer to identify a day (row) within the `meteo` data frame.

- latitude:

  Latitude (in degrees).

- elevation, slope, aspect:

  Elevation above sea level (in m), slope (in degrees) and aspect (in
  degrees from North).

- canopyEvaporation:

  Canopy evaporation (from interception) for `day` (mm).

- snowMelt:

  Snow melt values for `day` (mm).

- soilEvaporation:

  Bare soil evaporation for `day` (mm).

- herbTranspiration:

  Transpiration of herbaceous plants for `day` (mm).

- stepFunctions:

  An integer to indicate a simulation step for which photosynthesis and
  profit maximization functions are desired.

- modifyInput:

  Boolean flag to indicate that the input `x` object is allowed to be
  modified during the simulation.

## Value

A list with the following elements:

- `"cohorts"`: A data frame with cohort information, copied from
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- `"Stand"`: A vector of stand-level variables.

- `"Plants"`: A data frame of results for each plant cohort. When using
  `transp_transpirationGranier`, element `"Plants"` includes:

  - `"LAI"`: Leaf area index of the plant cohort.

  - `"LAIlive"`: Leaf area index of the plant cohort, assuming all
    leaves are unfolded.

  - `"AbsorbedSWRFraction"`: Fraction of SWR absorbed by each cohort.

  - `"Transpiration"`: Transpirated water (in mm) corresponding to each
    cohort.

  - `"GrossPhotosynthesis"`: Gross photosynthesis (in gC/m2)
    corresponding to each cohort.

  - `"psi"`: Water potential (in MPa) of the plant cohort (average over
    soil layers).

  - `"DDS"`: Daily drought stress \[0-1\] (relative whole-plant
    conductance).

  When using `transp_transpirationSperry` or
  `transp_transpirationSureau`, element `"Plants"` includes:

  - `"LAI"`: Leaf area index of the plant cohort.

  - `"LAIlive"`: Leaf area index of the plant cohort, assuming all
    leaves are unfolded.

  - `"Extraction"`: Water extracted from the soil (in mm) for each
    cohort.

  - `"Transpiration"`: Transpirated water (in mm) corresponding to each
    cohort.

  - `"GrossPhotosynthesis"`: Gross photosynthesis (in gC/m2)
    corresponding to each cohort.

  - `"NetPhotosynthesis"`: Net photosynthesis (in gC/m2) corresponding
    to each cohort.

  - `"RootPsi"`: Minimum water potential (in MPa) at the root collar.

  - `"StemPsi"`: Minimum water potential (in MPa) at the stem.

  - `"StemPLC"`: Proportion of conductance loss in stem.

  - `"LeafPsiMin"`: Minimum (predawn) water potential (in MPa) at the
    leaf (representing an average leaf).

  - `"LeafPsiMax"`: Maximum (midday) water potential (in MPa) at the
    leaf (representing an average leaf).

  - `"LeafPsiMin_SL"`: Minimum (predawn) water potential (in MPa) at
    sunlit leaves.

  - `"LeafPsiMax_SL"`: Maximum (midday) water potential (in MPa) at
    sunlit leaves.

  - `"LeafPsiMin_SH"`: Minimum (predawn) water potential (in MPa) at
    shade leaves.

  - `"LeafPsiMax_SH"`: Maximum (midday) water potential (in MPa) at
    shade leaves.

  - `"dEdP"`: Overall soil-plant conductance (derivative of the supply
    function).

  - `"DDS"`: Daily drought stress \[0-1\] (relative whole-plant
    conductance).

  - `"StemRWC"`: Relative water content of stem tissue (including
    symplasm and apoplasm).

  - `"LeafRWC"`: Relative water content of leaf tissue (including
    symplasm and apoplasm).

  - `"LFMC"`: Live fuel moisture content (in percent of dry weight).

  - `"WaterBalance"`: Plant water balance (extraction - transpiration).

- `"Extraction"`: A data frame with mm of water extracted from each soil
  layer (in columns) by each cohort (in rows). The sum of a given row is
  equal to the total extraction of the corresponding plant cohort.

- `"ExtractionPools"`: A named list with as many elements as plant
  cohorts, where each element is a matrix data with mm of water
  extracted from each layer (in columns) of the water pool of each
  cohort (in rows). The sum of a given matrix is equal to the total
  extraction of the corresponding plant cohort.

  The remaining items are only given by `transp_transpirationSperry` or
  `transp_transpirationSureau`:

- `"EnergyBalance"`: A list with the following elements:

  - `"Temperature"`: A data frame with the temperature of the atmosphere
    ('Tatm'), canopy ('Tcan') and soil ('Tsoil.1', 'Tsoil.2', ...) for
    each time step.

  - `"CanopyEnergyBalance"`: A data frame with the components of the
    canopy energy balance (in W/m2) for each time step.

  - `"SoilEnergyBalance"`: A data frame with the components of the soil
    energy balance (in W/m2) for each time step.

- `"RhizoPsi"`: Minimum water potential (in MPa) inside roots, after
  crossing rhizosphere, per cohort and soil layer.

- `"Sunlitleaves"` and `"ShadeLeaves"`: Data frames for sunlit leaves
  and shade leaves and the following columns per cohort:

  - `"LAI"`: Cumulative leaf area index of sunlit/shade leaves.

  - `"Vmax298"`: Average maximum carboxilation rate for sunlit/shade
    leaves.

  - `"Jmax298"`: Average maximum electron transport rate for
    sunlit/shade leaves.

- `"ExtractionInst"`: Water extracted by each plant cohort during each
  time step.

- `"PlantsInst"`: A list with instantaneous (per time step) results for
  each plant cohort:

  - `"E"`: A data frame with the cumulative transpiration (mm) for each
    plant cohort during each time step.

  - `"Ag"`: A data frame with the cumulative gross photosynthesis
    (gC/m2) for each plant cohort during each time step.

  - `"An"`: A data frame with the cumulative net photosynthesis (gC/m2)
    for each plant cohort during each time step.

  - `"Sunlitleaves"` and `"ShadeLeaves"`: Lists with instantaneous (for
    each time step) results for sunlit leaves and shade leaves and the
    following items:

    - `"Abs_SWR"`: A data frame with instantaneous absorbed short-wave
      radiation (SWR).

    - `"Net_LWR"`: A data frame with instantaneous net long-wave
      radiation (LWR).

    - `"An"`: A data frame with instantaneous net photosynthesis (in
      micromol/m2/s).

    - `"Ci"`: A data frame with instantaneous intercellular CO2
      concentration (in ppm).

    - `"GW"`: A data frame with instantaneous stomatal conductance (in
      mol/m2/s).

    - `"VPD"`: A data frame with instantaneous vapour pressure deficit
      (in kPa).

    - `"Temp"`: A data frame with leaf temperature (in degrees Celsius).

    - `"Psi"`: A data frame with leaf water potential (in MPa).

  - `"dEdP"`: A data frame with the slope of the plant supply function
    (an estimation of whole-plant conductance).

  - `"RootPsi"`: A data frame with root crown water potential (in MPa)
    for each plant cohort during each time step.

  - `"StemPsi"`: A data frame with stem water potential (in MPa) for
    each plant cohort during each time step.

  - `"LeafPsi"`: A data frame with leaf (average) water potential (in
    MPa) for each plant cohort during each time step.

  - `"StemPLC"`: A data frame with the proportion loss of conductance
    \[0-1\] for each plant cohort during each time step.

  - `"StemRWC"`: A data frame with the (average) relative water content
    of stem tissue \[0-1\] for each plant cohort during each time step.

  - `"LeafRWC"`: A data frame with the relative water content of leaf
    tissue \[0-1\] for each plant cohort during each time step.

  - `"StemSympRWC"`: A data frame with the (average) relative water
    content of symplastic stem tissue \[0-1\] for each plant cohort
    during each time step.

  - `"LeafSympRWC"`: A data frame with the relative water content of
    symplastic leaf tissue \[0-1\] for each plant cohort during each
    time step.

  - `"PWB"`: A data frame with plant water balance (extraction -
    transpiration).

- `"LightExtinction"`: A list of information regarding radiation balance
  through the canopy, as returned by function
  [`light_instantaneousLightExtinctionAbsortion`](https://emf-creaf.github.io/medfate/reference/light_advanced.md).

- `"CanopyTurbulence"`: Canopy turbulence (see
  [`wind_canopyTurbulence`](https://emf-creaf.github.io/medfate/reference/wind.md)).

- `"SupplyFunctions"`: If `stepFunctions` is not missing, a list of
  supply functions, photosynthesis functions and profit maximization
  functions.

## Details

Three sub-models are available:

- Sub-model in function `transp_transpirationGranier` was described in
  De Cáceres et al. (2015), and implements an approach originally
  described in Granier et al. (1999).

- Sub-model in function `transp_transpirationSperry` was described in De
  Cáceres et al. (2021), and implements a modelling approach originally
  described in Sperry et al. (2017).

- Sub-model in function `transp_transpirationSureau` was described for
  SurEau-Ecos v2.0 model in Ruffault et al. (2022).

## References

De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
forest inventory data to predict drought stress: the role of forest
structural changes vs. climate changes. Agricultural and Forest
Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).

De Cáceres M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L,
Poyatos R, Cabon A, Granda V, Forner A, Valladares F, Martínez-Vilalta J
(2021) Unravelling the effect of species mixing on water use and drought
stress in holm oak forests: a modelling approach. Agricultural and
Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).

Granier A, Bréda N, Biron P, Villette S (1999) A lumped water balance
model to evaluate duration and intensity of drought constraints in
forest stands. Ecol Modell 116:269–283.
https://doi.org/10.1016/S0304-3800(98)00205-1.

Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022)
SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations
of plant water status and drought-induced mortality at the ecosystem
level. Geoscientific Model Development 15, 5593-5626
(doi:10.5194/gmd-15-5593-2022).

Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S.
Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to
the environment from the optimization of photosynthetic gain and
hydraulic cost. Plant Cell and Environment 40, 816-830 (doi:
10.1111/pce.12852).

## See also

[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
[`plot.spwb_day`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md)

## Author

- Miquel De Cáceres Ainsa, CREAF

- Nicolas Martin-StPaul, URFM-INRAE

## Examples

``` r
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

#Initialize control parameters
control <- defaultControl("Granier")

#Initialize input
x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

# Transpiration according to Granier's model, plant water potential 
# and plant stress for a given day
t1 <- transp_transpirationGranier(x1, examplemeteo, 1, 
                                 latitude = 41.82592, elevation = 100, slope = 0, aspect = 0, 
                                 modifyInput = FALSE)

#Switch to 'Sperry' transpiration mode
control <- defaultControl("Sperry")

#Initialize input
x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

# Transpiration according to Sperry's model
t2 <- transp_transpirationSperry(x2, examplemeteo, 1, 
                                latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
                                modifyInput = FALSE)
                                
#Switch to 'Sureau' transpiration mode
control <- defaultControl("Sureau")

#Initialize input
x3 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

# Transpiration according to Sureau model
t3 <- transp_transpirationSureau(x3, examplemeteo, 1, 
                                  latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
                                  modifyInput = FALSE)
                                
```
