# Fire hazard

Estimates potential fire behaviour at each daily step of a simulation

## Usage

``` r
fireHazard(
  x,
  SpParams,
  forest = NULL,
  standardConditions = FALSE,
  freq = "days",
  fun = "max"
)
```

## Arguments

- x:

  An object of class
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/pwb.md),
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
  [`growth_day`](https://emf-creaf.github.io/medfate/reference/growth_day.md)
  or
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

- SpParams:

  A data frame with species parameters (see
  [`SpParamsDefinition`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  and
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- forest:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
  (needed if `x` is not of class
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md)).

- standardConditions:

  A logical flag to indicate that standard fire weather conditions are
  to be used (instead of deriving fuel moisture and windspeed from `x`).

- freq:

  Frequency of summary statistics (see
  [`cut.Date`](https://rdrr.io/r/base/cut.POSIXt.html)).

- fun:

  Summary function (by default, maximum values).

## Value

A matrix with fire behaviour variables (columns) for each simulated day
(rows) or coarser time steps if summaries are requested.

## Details

Live fuel moisture of shrub and canopy layers is estimated from plant
water status. Dead fuel moisture is estimated following Resco-de-Dios et
al. (2015).

## References

Resco de Dios, V., A. W. Fellows, R. H. Nolan, M. M. Boer, R. A.
Bradstock, F. Domingo, and M. L. Goulden. 2015. A semi-mechanistic model
for predicting the moisture content of fine litter. Agricultural and
Forest Meteorology 203:64–73.

Ruffault J, Limousin JM, Pimont F, Dupuy JL, De Cáceres M, Cochard H,
Mouillot F, Blackman C, Torres-Ruiz JM, Parsons R, Moreno M, Delzon S,
Jansen S, Olioso A, Choat B, Martin-StPaul N. 2023. Plant hydraulic
modelling of leaf and canopy fuel moisture content reveals increasing
vulnerability of a Mediterranean forest to wildfires under extreme
drought. New Phytologist. (10.1111/nph.18614).

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`fuel_FCCS`](https://emf-creaf.github.io/medfate/reference/fuel_properties.md),
[`fire_FCCS`](https://emf-creaf.github.io/medfate/reference/fire_behaviour.md)

## Author

Miquel De Cáceres Ainsa, CREAF

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
#> Initial plant water content (mm): 4.69853
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 4.69659
#> Final soil water content (mm): 275.04
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00193896
#> Plant water balance result (mm): -0.00196771
#> Change in soil water content (mm): -15.8347
#> Soil water balance result (mm): -15.8347
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 370
#>   Infiltration (mm) 402 Infiltration excess (mm) 20 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 24  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 245
#>   Plant extraction from soil (mm) 245  Plant water balance (mm) -0 Hydraulic redistribution (mm) 3
#>   Runoff (mm) 20 Deep drainage (mm) 136

#Evaluate fire hazard
F1 <- fireHazard(S1, SpParamsMED, exampleforest)
# }
```
