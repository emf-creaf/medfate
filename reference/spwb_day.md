# Single-day soil-plant water balance

Function `spwb_day` performs water balance for a single day.

## Usage

``` r
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

- `"EnergyBalance"`: Energy balance of the stand (only returned when
  `transpirationMode = "Sperry"` or `transpirationMode = "Sureau"`; see
  [`transp_transpirationSperry`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)).

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

The following additional items are only returned when
`transpirationMode = "Sperry"` or `transpirationMode = "Sureau"`:

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

Detailed model description is available in the medfate book.

Soil-plant water balance simulations allow using different sub-models
for bulk soil water flows and different sub-models of transpiration and
photosynthesis:

\(1\) Sub-models of transpiration and photosynthesis (control parameter
`transpirationMode`):

- The sub-model corresponding to 'Granier' transpiration mode is
  illustrated by internal function
  [`transp_transpirationGranier`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)
  and was described in De Caceres et al. (2015), and implements an
  approach originally described in Granier et al. (1999).

- The sub-model corresponding to 'Sperry' transpiration mode is
  illustrated by internal function
  [`transp_transpirationSperry`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)
  and was described in De Caceres et al. (2021), and implements a
  modelling approach originally described in Sperry et al. (2017).

- The sub-model corresponding to 'Sureau' transpiration mode is
  illustrated by internal function
  [`transp_transpirationSureau`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)
  and was described for model SurEau-Ecos v2.0 in Ruffault et al.
  (2022).

Simulations using the 'Sperry' or 'Sureau' transpiration mode are
computationally much more expensive than 'Granier' because they include
explicit plant hydraulics and canopy/soil energy balance at subdaily
time steps.

\(2\) Sub-models of bulk soil water flows (control parameter
`soilDomains`; see internal function
[`hydrology_soilWaterBalance`](https://emf-creaf.github.io/medfate/reference/hydrology_soilWaterBalance.md)):

- The submodel corresponding to 'buckets' is a multi-bucket model
  adds/substracts water to each layer and if content is above field
  capacity the excess percolates to the layer below.

- The submodel corresponding to 'single' is a single-domain solving 1D
  Richards equation (Bonan et al. 2019).

- The submodel corresponding to 'dual' is a dual-permeability model from
  MACRO 5.0 (Jarvis et al. 1991; Larsbo et al. 2005)

## References

Bonan, G. (2019). Climate change and terrestrial ecosystem modeling.
Cambridge University Press, Cambridge, UK.

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

Jarvis, N.J., Jansson, P‐E., Dik, P.E. & Messing, I. (1991). Modelling
water and solute transport in macroporous soil. I. Model description and
sensitivity analysis. Journal of Soil Science, 42, 59–70.

Larsbo, M., Roulier, S., Stenemo, F., Kasteel, R. & Jarvis, N. (2005).
An Improved Dual‐Permeability Model of Water Flow and Solute Transport
in the Vadose Zone. Vadose Zone Journal, 4, 398–406.

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
[`growth_day`](https://emf-creaf.github.io/medfate/reference/growth_day.md)

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

```
