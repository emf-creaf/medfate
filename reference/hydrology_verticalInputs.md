# Water vertical inputs

High-level functions to define water inputs into the soil of a stand:

- Function `hydrology_waterInputs` performs canopy water interception
  and snow accumulation/melt.

- Function `hydrology_snowMelt` estimates snow melt using a simple
  energy balance, according to Kergoat (1998).

## Usage

``` r
hydrology_snowMelt(tday, rad, LgroundSWR, elevation)

hydrology_waterInputs(
  x,
  prec,
  rainfallIntensity,
  pet,
  tday,
  rad,
  elevation,
  Cm,
  LgroundPAR,
  LgroundSWR,
  modifyInput = TRUE
)
```

## Arguments

- tday:

  Average day temperature (ºC).

- rad:

  Solar radiation (in MJ/m2/day).

- LgroundSWR:

  Percentage of short-wave radiation (SWR) reaching the ground.

- elevation:

  Altitude above sea level (m).

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  or
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- prec:

  Precipitation for the given day (mm)

- rainfallIntensity:

  Rainfall intensity rate (mm/h).

- pet:

  Potential evapotranspiration for the given day (mm)

- Cm:

  Canopy water storage capacity.

- LgroundPAR:

  Percentage of photosynthetically-active radiation (PAR) reaching the
  ground.

- modifyInput:

  Boolean flag to indicate that the input `x` object should be modified
  during the simulation.

## Value

Function `hydrology_waterInputs` returns a named vector with the
following elements, all in mm:

- Rain:

  Precipitation as rainfall.

- Snow:

  Precipitation as snow.

- Interception:

  Rainfall water intercepted by the canopy and evaporated.

- Snowmelt:

  Snow melted during the day, and added to the water infiltrated.

- NetRain:

  Rainfall reaching the ground.

## Details

The function simulates different vertical hydrological processes, which
are described separately in other functions. If `modifyInput = TRUE` the
function will modify the `x` object (including both soil moisture and
the snowpack on its surface) as a result of simulating hydrological
processes.

## References

Kergoat L. (1998). A model for hydrological equilibrium of leaf area
index on a global scale. Journal of Hydrology 212–213: 268–286.

## See also

[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
[`hydrology_rainInterception`](https://emf-creaf.github.io/medfate/reference/hydrology_interception.md),
[`hydrology_soilEvaporation`](https://emf-creaf.github.io/medfate/reference/hydrology_soilEvaporation.md)

## Author

Miquel De Cáceres Ainsa, CREAF
