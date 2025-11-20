# Bare soil evaporation and herbaceous transpiration

Functions:

- Function `hydrology_soilEvaporationAmount` calculates the amount of
  evaporation from bare soil, following Ritchie (1972).

- Function `hydrology_soilEvaporation` calculates the amount of
  evaporation from bare soil and distributes it among soil layers.

- Function `hydrology_herbaceousTranspiration` calculates the amount of
  transpiration due to herbaceous plants.

## Usage

``` r
hydrology_soilEvaporationAmount(DEF, PETs, Gsoil)

hydrology_soilEvaporation(
  soil,
  snowpack,
  soilFunctions,
  pet,
  LgroundSWR,
  modifySoil = TRUE
)

hydrology_herbaceousTranspiration(
  pet,
  LherbSWR,
  herbLAI,
  soil,
  soilFunctions,
  modifySoil = TRUE
)
```

## Arguments

- DEF:

  Water deficit in the (topsoil) layer.

- PETs:

  Potential evapotranspiration at the soil surface.

- Gsoil:

  Gamma parameter (maximum daily evaporation).

- soil:

  An object of class
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md).

- snowpack:

  The amount of snow (in water equivalents, mm) in the snow pack.

- soilFunctions:

  Soil water retention curve and conductivity functions, either 'SX'
  (for Saxton) or 'VG' (for Van Genuchten).

- pet:

  Potential evapotranspiration for a given day (mm)

- LgroundSWR:

  Percentage of short-wave radiation (SWR) reaching the ground.

- modifySoil:

  Boolean flag to indicate that the input `soil` object should be
  modified during the simulation.

- LherbSWR:

  Percentage of short-wave radiation (SWR) reaching the herbaceous
  layer.

- herbLAI:

  Leaf area index of the herbaceous layer.

## Value

Function `hydrology_soilEvaporationAmount` returns the amount of water
evaporated from the soil.

Function `hydrology_soilEvaporation` returns a vector of water
evaporated from each soil layer.

## References

Ritchie (1972). Model for predicting evaporation from a row crop with
incomplete cover. - Water resources research.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`hydrology_waterInputs`](https://emf-creaf.github.io/medfate/reference/hydrology_verticalInputs.md),
[`hydrology_infiltration`](https://emf-creaf.github.io/medfate/reference/hydrology_infiltration.md)

## Author

Miquel De Cáceres Ainsa, CREAF
