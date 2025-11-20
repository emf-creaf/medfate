# Herbaceous description functions

Functions to calculate attributes of the herbaceous component of a
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object

## Usage

``` r
herb_foliarBiomass(x, SpParams)

herb_fuelLoading(x, SpParams)

herb_LAI(x, SpParams)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

## Value

A single scalar:

- `herb_foliarBiomass`: Herbaceous biomass of leaves (in kg/m2).

- `herb_fuelLoading`: Herbaceous fine fuel loading (in kg/m2).

- `herb_LAI`: Herbaceous leaf area index (m2/m2).

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`plant_basalArea`](https://emf-creaf.github.io/medfate/reference/plant_values.md),
[`summary.forest`](https://emf-creaf.github.io/medfate/reference/summary.forest.md)

## Author

Miquel De Cáceres Ainsa, CREAF
