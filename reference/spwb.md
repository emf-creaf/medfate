# Soil-plant water balance

Function `spwb()` is a water balance model that determines changes in
soil moisture, soil water potentials, plant transpiration and drought
stress at daily steps for a given forest stand during a period specified
in the input climatic data. Plant transpiration and photosynthesis
processes are conducted with different level of detail depending on the
transpiration mode.

## Usage

``` r
spwb(
  x,
  meteo,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  CO2ByYear = numeric(0),
  waterTableDepth = NA_real_
)
```

## Arguments

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- meteo:

  A data frame with daily meteorological data series. Row names of the
  data frame should correspond to date strings with format "yyyy-mm-dd"
  (see [`Date`](https://rdrr.io/r/base/Dates.html)). Alternatively, a
  column called `"dates"` or `"Dates"` can contain
  [`Date`](https://rdrr.io/r/base/Dates.html) or
  [`POSIXct`](https://rdrr.io/r/base/DateTimeClasses.html) classes. The
  following columns are required and cannot have missing values:

  - `MinTemperature`: Minimum temperature (in degrees Celsius).

  - `MaxTemperature`: Maximum temperature (in degrees Celsius).

  - `Precipitation`: Precipitation (in mm).

  The following columns are required but can contain missing values
  (NOTE: missing values will raise warnings):

  - `MinRelativeHumidity`: Minimum relative humidity (in percent).

  - `MaxRelativeHumidity`: Maximum relative humidity (in percent).

  - `Radiation`: Solar radiation (in MJ/m2/day).

  The following columns are optional:

  - `WindSpeed`: Above-canopy wind speed (in m/s). This column may not
    exist, or can be left with `NA` values. In both cases simulations
    will assume a constant value specified in
    [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md).

  - `CO2`: Atmospheric (above-canopy) CO2 concentration (in ppm). This
    column may not exist, or can be left with `NA` values. In both cases
    simulations will assume a constant value specified in
    [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md).

  - `Patm`: Atmospheric pressure (in kPa). This column may not exist, or
    can be left with `NA` values. In both cases, a value is estimated
    from elevation.

- latitude:

  Latitude (in degrees).

- elevation, slope, aspect:

  Elevation above sea level (in m), slope (in degrees) and aspect (in
  degrees from North).

- CO2ByYear:

  A named numeric vector with years as names and atmospheric CO2
  concentration (in ppm) as values. Used to specify annual changes in
  CO2 concentration along the simulation (as an alternative to
  specifying daily values in `meteo`).

- waterTableDepth:

  Water table depth (in mm). When not missing, capillarity rise will be
  allowed if lower than total soil depth.

## Value

Function `spwb` returns a list of class 'spwb'. Since lists are
difficult to handle, we recommend using function
[`extract`](https://emf-creaf.github.io/medfate/reference/extract.md) to
reshape simulation results (including their units) from those objects.

List elements are as follows:

- `"latitude"`: Latitude (in degrees) given as input.

- `"topography"`: Vector with elevation, slope and aspect given as
  input.

- `"weather"`: A copy of the input weather data frame.

- `"spwbInput"`: An copy of the object `x` of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  given as input.

- `"spwbOutput"`: An copy of the final state of the object `x` of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- `"WaterBalance"`: A data frame where different variables (in columns)
  are given for each simulated day (in rows):

  - `"PET"`: Potential evapotranspiration (in mm).

  - `"Precipitation"`: Input precipitation (in mm).

  - `"Rain"`: Precipitation as rainfall (in mm).

  - `"Snow"`: Precipitation as snowfall (in mm).

  - `"NetRain"`: Net rain, after accounting for interception (in mm).

  - `"Infiltration"`: The amount of water infiltrating into the soil (in
    mm).

  - `"InfiltrationExcess"`: Excess infiltration in the topmost layer
    leading to an increase in runoff (in mm).

  - `"SaturationExcess"`: Excess saturation in the topmost layer leading
    to an increase in runoff (in mm).

  - `"CapillarityRise"`: Water entering the soil via capillarity
    rise (mm) from the water table, if `waterTableDepth` is supplied.

  - `"Runoff"`: The amount of water exported via surface runoff (in mm).

  - `"DeepDrainage"`: The amount of water exported via deep drainage (in
    mm).

  - `"Evapotranspiration"`: Evapotranspiration (in mm).

  - `"SoilEvaporation"`: Bare soil evaporation (in mm).

  - `"HerbTranspiration"`: Transpiration due to the herbaceous layer (in
    mm).

  - `"PlantExtraction"`: Amount of water extracted from soil by woody
    plants (in mm).

  - `"Transpiration"`: Woody plant transpiration (in mm).

  - `"HydraulicRedistribution"`: Water redistributed among soil layers,
    transported through the plant hydraulic network.

- `"EnergyBalance"`: A data frame with the daily values of energy
  balance components for the soil and the canopy (only for
  `transpirationMode = "Sperry"` or `transpirationMode = "Sureau"`).

- `"Temperature"`: A data frame with the daily values of
  minimum/mean/maximum temperatures for the atmosphere (input), canopy
  and soil (only for `transpirationMode = "Sperry"` or
  `transpirationMode = "Sureau"`).

- `"Soil"`: A list with the following subelements:

  - `"SWC"`: Soil water content (percent of soil volume) in each soil
    layer (and overall).

  - `"RWC"`: Relative soil moisture content (relative to field capacity)
    in each soil layer (and overall).

  - `"REW"`: Relative extractable water (min. psi = -5 MPa) in each soil
    layer (and overall).

  - `"ML"`: Soil water volume in each soil layer (in L/m2) (and
    overall).

  - `"Psi"`: Soil water potential in each soil layer (in MPa) (and
    overall).

  - `"PlantExt"`: Plant extraction from each soil layer (in mm) (and
    overall).

  - `"HydraulicInput"`: Water that entered the layer coming from other
    layers and transported via the plant hydraulic network (in mm) (and
    overall).

- `"Snow"`: A data frame where the following variable (in columns) is
  given for each simulated day (in rows):

  - `"SWE"`: Snow water equivalent (mm) of the snow pack.

- `"Stand"`: A data frame where different variables (in columns) are
  given for each simulated day (in rows):

  - `"LAI"`: LAI of the stand (including the herbaceous layer and live +
    dead leaves of woody plants) (in m2/m2).

  - `"LAIherb"`: LAI of the herbaceous layer (in m2/m2).

  - `"LAIlive"`: LAI of the woody plants assuming all leaves are
    unfolded (in m2/m2).

  - `"LAIexpanded"`: LAI of the woody plants with leaves actually
    unfolded (in m2/m2).

  - `"LAIdead"`: LAI of the woody plants corresponding to dead leaves
    (in m2/m2).

  - `"Cm"`: Water retention capacity of the canopy (in mm) (accounting
    for leaf phenology).

  - `"LgroundPAR"`: The percentage of PAR that reaches the ground
    (accounting for leaf phenology).

  - `"LgroundSWR"`: The percentage of SWR that reaches the ground
    (accounting for leaf phenology).

- `"Plants"`: A list of daily results for plant cohorts (see below).

- `"SunlitLeaves"`: A list of daily results for sunlit leaves of plant
  cohorts (only for `transpirationMode = "Sperry"` or
  `transpirationMode = "Sureau"`; see below).

- `"ShadeLeaves"`: A list of daily results for sunlit leaves of plant
  cohorts (only for `transpirationMode = "Sperry"` or
  `transpirationMode = "Sureau"`; see below).

- `"subdaily"`: A list of objects of class
  [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
  one per day simulated (only if required in `control` parameters, see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

When `transpirationMode = "Granier"`, element `"Plants"` is a list with
the following subelements:

- `"LAI"`: A data frame with the daily leaf area index for each plant
  cohort.

- `"LAIlive"`: A data frame with the daily leaf area index for each
  plant cohort, assuming all leaves are unfolded (in m2/m2).

- `"FPAR"`: A data frame with the fraction of PAR at the canopy level of
  each plant cohort.

- `"AbsorbedSWRFraction"`: A data frame with the fraction of SWR
  absorbed by each plant cohort.

- `"Transpiration"`: A data frame with the amount of daily transpiration
  (in mm) for each plant cohort.

- `"GrossPhotosynthesis"`: A data frame with the amount of daily gross
  photosynthesis (in g C·m-2) for each plant cohort.

- `"PlantPsi"`: A data frame with the average daily water potential of
  each plant (in MPa).

- `"LeafPLC"`: A data frame with the average daily proportion of leaf
  conductance loss of each plant (\[0-1\]).

- `"StemPLC"`: A data frame with the average daily proportion of stem
  conductance loss of each plant (\[0-1\]).

- `"PlantWaterBalance"`: A data frame with the daily balance between
  transpiration and soil water extraction for each plant cohort.

- `"LeafRWC"`: A data frame with the average daily leaf relative water
  content of each plant (in percent).

- `"StemRWC"`: A data frame with the average daily stem relative water
  content of each plant (in percent).

- `"LFMC"`: A data frame with the daily live fuel moisture content (in
  percent of dry weight).

- `"PlantStress"`: A data frame with the amount of daily stress \[0-1\]
  suffered by each plant cohort (relative whole-plant conductance).

If `transpirationMode="Sperry"` or `transpirationMode="Sureau"`, element
`"Plants"` is a list with the following subelements:

- `"LAI"`: A data frame with the daily leaf area index for each plant
  cohort.

- `"LAIlive"`: A data frame with the daily leaf area index for each
  plant cohort, assuming all leaves are unfolded (in m2/m2).

- `"FPAR"`: A data frame with the fraction of PAR at the canopy level of
  each plant cohort.

- `"AbsorbedSWR"`: A data frame with the daily SWR absorbed by each
  plant cohort.

- `"NetLWR"`: A data frame with the daily net LWR by each plant cohort.

- `"Transpiration"`: A data frame with the amount of daily transpiration
  (in mm) for each plant cohorts.

- `"GrossPhotosynthesis"`: A data frame with the amount of daily gross
  photosynthesis (in g C·m-2) for each plant cohort.

- `"NetPhotosynthesis"`: A data frame with the amount of daily net
  photosynthesis (in g C·m-2) for each plant cohort.

- `"dEdP"`: A data frame with mean daily values of soil-plant
  conductance (derivative of the supply function) for each plant cohort.

- `"PlantWaterBalance"`: A data frame with the daily balance between
  transpiration and soil water extraction for each plant cohort.

- `"LeafPsiMin"`: A data frame with the minimum (midday) daily (average)
  leaf water potential of each plant (in MPa).

- `"LeafPsiMax"`: A data frame with the maximum (predawn) daily
  (average) leaf water potential of each plant (in MPa).

- `"LeafRWC"`: A data frame with the average daily leaf relative water
  content of each plant (in percent).

- `"StemRWC"`: A data frame with the average daily stem relative water
  content of each plant (in percent).

- `"LFMC"`: A data frame with the daily live fuel moisture content (in
  percent of dry weight).

- `"StemPsi"`: A data frame with the minimum daily stem water potential
  of each plant (in MPa).

- `"LeafPLC"`: A data frame with the average daily proportion of leaf
  conductance loss of each plant (\[0-1\]).

- `"StemPLC"`: A data frame with the average daily proportion of stem
  conductance loss of each plant (\[0-1\]).

- `"RootPsi"`: A data frame with the minimum daily root water potential
  of each plant (in MPa).

- `"RhizoPsi"`: A list of data frames (one per plant cohort) with the
  minimum daily root water potential of each plant (in MPa).

- `"LFMC"`: A data frame with the daily live fuel moisture content (in
  percent of dry weight).

- `"PlantStress"`: A data frame with the amount of daily stress \[0-1\]
  suffered by each plant cohort (relative whole-plant conductance).

If `transpirationMode="Sperry"` or `transpirationMode="Sureau"`, then
elements `"SunlitLeaves"` and `"ShadeLeaves"` are list with daily
results for sunlit and shade leaves:

- `"PsiMin"`: A data frame with the minimum (midday) daily sunlit or
  shade leaf water potential (in MPa).

- `"PsiMax"`: A data frame with the maximum (predawn) daily sunlit or
  shade leaf water potential (in MPa).

- `"TempMin"`: A data frame with the minimum daily sunlit or shade leaf
  temperature (in Celsius).

- `"TempMax"`: A data frame with the maximum daily sunlit or shade leaf
  temperature (in Celsius).

- `"GSWMin"`: A data frame with the minimum daily sunlit or shade leaf
  stomatal conductance (in mol·m-2·s-1).

- `"GSWMax"`: A data frame with the maximum daily sunlit or shade leaf
  stomatal conductance (in mol·m-2·s-1).

## Details

The simulation functions allow using three different sub-models of
transpiration and photosynthesis:

- The sub-model corresponding to 'Granier' transpiration mode is
  illustrated by function
  [`transp_transpirationGranier`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)
  and was described in De Caceres et al. (2015), and implements an
  approach originally described in Granier et al. (1999).

- The sub-model corresponding to 'Sperry' transpiration mode is
  illustrated by function
  [`transp_transpirationSperry`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)
  and was described in De Caceres et al. (2021), and implements a
  modelling approach originally described in Sperry et al. (2017).

- The sub-model corresponding to 'Sureau' transpiration mode is
  illustrated by function
  [`transp_transpirationSureau`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)
  and was described for model SurEau-Ecos v2.0 in Ruffault et al.
  (2022).

Simulations using the 'Sperry' or 'Sureau' transpiration mode are
computationally much more expensive than 'Granier' because they include
explicit plant hydraulics and canopy/soil energy balance at subdaily
time steps.

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

[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
[`plot.spwb`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md),
[`extract`](https://emf-creaf.github.io/medfate/reference/extract.md),
[`summary.spwb`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md),
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`pwb`](https://emf-creaf.github.io/medfate/reference/pwb.md),
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
[`aspwb`](https://emf-creaf.github.io/medfate/reference/aspwb.md)

## Author

- Miquel De Cáceres Ainsa, CREAF

- Nicolas Martin-StPaul, URFM-INRAE

## Examples

``` r
# \donttest{
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

#Call simulation function
S1 <- spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant water content (mm): 4.73001
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 4.72839
#> Final soil water content (mm): 274.93
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00162134
#> Plant water balance result (mm): -0.00163359
#> Change in soil water content (mm): -15.9454
#> Soil water balance result (mm): -15.9454
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 0
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 370
#>   Infiltration (mm) 401 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 25  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 247
#>   Plant extraction from soil (mm) 247  Plant water balance (mm) -0 Hydraulic redistribution (mm) 3
#>   Runoff (mm) 21 Deep drainage (mm) 131

#Switch to 'Sperry' transpiration mode
control <- defaultControl("Sperry")

#Initialize input
x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
S2 <- spwb(x2, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant water content (mm): 6.78662
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 6.78165
#> Final soil water content (mm): 273.985
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00497469
#> Plant water balance result (mm): 2.04696e-16
#> Change in soil water content (mm): -16.8895
#> Soil water balance result (mm): -16.8895
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 370
#>   Infiltration (mm) 402 Infiltration excess (mm) 19 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 22  Herbaceous transpiration (mm) 13 Woody plant transpiration (mm) 242
#>   Plant extraction from soil (mm) 242  Plant water balance (mm) 0 Hydraulic redistribution (mm) 4
#>   Runoff (mm) 19 Deep drainage (mm) 142

#Switch to 'Sureau' transpiration mode
control <- defaultControl("Sureau")

#Initialize input
x3 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
S3 <- spwb(x3, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant water content (mm): 6.78662
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 6.76076
#> Final soil water content (mm): 278.244
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.0258564
#> Plant water balance result (mm): -0.111014
#> Change in soil water content (mm): -12.631
#> Soil water balance result (mm): -12.631
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 371
#>   Infiltration (mm) 401 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 30  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 187
#>   Plant extraction from soil (mm) 187  Plant water balance (mm) -0 Hydraulic redistribution (mm) 0
#>   Runoff (mm) 21 Deep drainage (mm) 183
# }
                
```
