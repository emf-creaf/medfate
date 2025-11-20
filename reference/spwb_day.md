# Single-day simulation

Function `spwb_day` performs water balance for a single day and
`growth_day` performs water and carbon balance for a single day.

## Usage

``` r
growth_day(
  x,
  date,
  meteovec,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  runon = 0,
  lateralFlows = NULL,
  waterTableDepth = NA_real_,
  modifyInput = TRUE
)

spwb_day(
  x,
  date,
  meteovec,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  runon = 0,
  lateralFlows = NULL,
  waterTableDepth = NA_real_,
  modifyInput = TRUE
)
```

## Arguments

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  or
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- date:

  Date as string "yyyy-mm-dd".

- meteovec:

  A named numerical vector with weather data. See variable names in
  parameter `meteo` of
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md).

- latitude:

  Latitude (in degrees).

- elevation, slope, aspect:

  Elevation above sea level (in m), slope (in degrees) and aspect (in
  degrees from North).

- runon:

  Surface water amount running on the target area from upslope (in mm).

- lateralFlows:

  Lateral source/sink terms for each soil layer (interflow/to from
  adjacent locations) as mm/day.

- waterTableDepth:

  Water table depth (in mm). When not missing, capillarity rise will be
  allowed if lower than total soil depth.

- modifyInput:

  Boolean flag to indicate that the input `x` object is allowed to be
  modified during the simulation.

## Value

Function `spwb_day()` returns a list of class `spwb_day` with the
following elements:

- `"cohorts"`: A data frame with cohort information, copied from
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- `"topography"`: Vector with elevation, slope and aspect given as
  input.

- `"weather"`: A vector with the input weather.

- `"WaterBalance"`: A vector of water balance components (rain, snow,
  net rain, infiltration, ...) for the simulated day, equivalent to one
  row of 'WaterBalance' object given in
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md).

- `"Soil"`: A data frame with results for each soil layer:

  - `"Psi"`: Soil water potential (in MPa) at the end of the day.

  - `"HerbTranspiration"`: Water extracted by herbaceous plants from
    each soil layer (in mm).

  - `"HydraulicInput"`: Water entering each soil layer from other
    layers, transported via plant roots (in mm).

  - `"HydraulicOutput"`: Water leaving each soil layer (going to other
    layers or the transpiration stream) (in mm).

  - `"PlantExtraction"`: Water extracted by woody plants from each soil
    layer (in mm).

- `"Stand"`: A named vector with with stand values for the simulated
  day, equivalent to one row of 'Stand' object returned by
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md).

- `"Plants"`: A data frame of results for each plant cohort (see
  [`transp_transpirationGranier`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)
  or
  [`transp_transpirationSperry`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)).

The following items are only returned when
`transpirationMode = "Sperry"` or `transpirationMode = "Sureau"`:

- `"EnergyBalance"`: Energy balance of the stand (see
  [`transp_transpirationSperry`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)).

- `"RhizoPsi"`: Minimum water potential (in MPa) inside roots, after
  crossing rhizosphere, per cohort and soil layer.

- `"SunlitLeaves"` and `"ShadeLeaves"`: For each leaf type, a data frame
  with values of LAI, Vmax298 and Jmax298 for leaves of this type in
  each plant cohort.

- `"ExtractionInst"`: Water extracted by each plant cohort during each
  time step.

- `"PlantsInst"`: A list with instantaneous (per time step) results for
  each plant cohort (see
  [`transp_transpirationSperry`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)).

- `"LightExtinction"`: A list of information regarding radiation balance
  through the canopy, as returned by function
  [`light_instantaneousLightExtinctionAbsortion`](https://emf-creaf.github.io/medfate/reference/light_advanced.md).

- `"CanopyTurbulence"`: Canopy turbulence (see
  [`wind_canopyTurbulence`](https://emf-creaf.github.io/medfate/reference/wind.md)).

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
computationally much more expensive than 'Granier'.

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
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`plot.spwb_day`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md),
[`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
[`plot.growth_day`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md)

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

#Define soil parameters
examplesoil <- defaultSoilParams(4)

# Day to be simulated
d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

#Simulate water balance one day only (Granier mode)
control <- defaultControl("Granier")
x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
sd1 <- spwb_day(x1, date, meteovec,  
                latitude = 41.82592, elevation = 100, slope=0, aspect=0) 

#Simulate water balance for one day only (Sperry mode)
control <- defaultControl("Sperry")
x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
sd2 <-spwb_day(x2, date, meteovec,
              latitude = 41.82592, elevation = 100, slope=0, aspect=0)

#Simulate water balance for one day only (Sureau mode)
control <- defaultControl("Sureau")
x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
sd3 <-spwb_day(x3, date, meteovec,
              latitude = 41.82592, elevation = 100, slope=0, aspect=0)


#Simulate water and carbon balance for one day only (Granier mode)
control <- defaultControl("Granier")
x4  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
sd4 <- growth_day(x4, date, meteovec,
                latitude = 41.82592, elevation = 100, slope=0, aspect=0)

#Simulate water and carbon balance for one day only (Sperry mode)
control <- defaultControl("Sperry")
x5  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
sd5 <- growth_day(x5, date, meteovec,
                latitude = 41.82592, elevation = 100, slope=0, aspect=0)

#Simulate water and carbon balance for one day only (Sureau mode)
control <- defaultControl("Sureau")
x6  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
sd6 <- growth_day(x6, date, meteovec,
                latitude = 41.82592, elevation = 100, slope=0, aspect=0)
```
