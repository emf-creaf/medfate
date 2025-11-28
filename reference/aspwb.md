# Simulation in agricultural areas

Function `aspwb_day` performs water balance for a single day in an
agriculture location. Function `aspwb` performs water balance for
multiple days in an agriculture location.

## Usage

``` r
aspwbInput(crop_factor, control, soil)

aspwb_day(
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

aspwb(
  x,
  meteo,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  waterTableDepth = NA_real_
)
```

## Arguments

- crop_factor:

  Agriculture crop factor.

- control:

  A list with default control parameters (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- soil:

  An object of class
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md).

- x:

  An object of class `aspwbInput`.

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

- meteo:

  A data frame with daily meteorological data series (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

## See also

[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
control <- defaultControl()
examplesoil <- defaultSoilParams(4)

x <- aspwbInput(0.75, control, examplesoil)

# Day to be simulated
d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

#Call simulation function for a single days
sd <- aspwb_day(x, date, meteovec,  
               latitude = 41.82592, elevation = 100) 
#> Package 'meteoland' [ver. 2.2.4]

#Call simulation function for multiple days
S <- aspwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial soil water content (mm): 287.448
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:....................................
#> 
#> Final soil water content (mm): 244.902
#> Final snowpack content (mm): 0
#> Change in soil water content (mm): -42.5459
#> Soil water balance result (mm): -42.5459
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 0
#> Water balance components:
#>   Precipitation (mm) 513
#>   Rain (mm) 462 Snow (mm) 51
#>   Infiltration (mm) 507 Infiltration excess (mm) 6 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 4 Transpiration (mm) 446
#>   Runoff (mm) 6 Deep drainage (mm) 100
```
