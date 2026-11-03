# Stand values

Functions to calculate stand attributes of a
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object.

## Usage

``` r
stand_basalArea(x, minDBH = 7.5)

stand_foliarBiomass(x, SpParams, gdd = NA_real_)

stand_fuelLoading(x, SpParams, gdd = NA_real_, includeDead = TRUE)

stand_shrubVolume(x, SpParams)

stand_LAI(x, SpParams, gdd = NA_real_, bounded = TRUE)

stand_dominantTreeDiameter(x, minDBH = 7.5)

stand_treeDensity(x, minDBH = 7.5)

stand_meanTreeHeight(x, minDBH = 7.5)

stand_dominantTreeHeight(x, minDBH = 7.5)

stand_hartBeckingIndex(x, minDBH = 7.5)

stand_quadraticMeanTreeDiameter(x, minDBH = 7.5)

stand_dominantTreeSpecies(x, SpParams)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- minDBH:

  Minimum diameter at breast height (in cm) to include in estimation.

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- gdd:

  Growth degree days (to account for leaf phenology effects).

- includeDead:

  A flag to indicate that standing dead fuels (dead branches) are
  included.

- bounded:

  A boolean flag to indicate that extreme values should be prevented
  (maximum tree LAI = 7 and maximum shrub LAI = 3)

## Value

- `stand_basalArea`: Stand basal area (m2/ha).

- `stand_treeDensity`: Stand tree density (in ind/ha).

- `stand_dominantTreeDiameter`: Dominant tree diameter, i.e the average
  diameter of the 100 widest trees (in cm).

- `stand_quadraticMeanTreeDiameter`: Quadratic mean tree diameter, i.e.
  the diameter value corresponding to the current basal area and
  density.

- `stand_meanTreeHeight`: Mean tree height (in cm).

- `stand_dominantTreeHeight`: Dominant tree height, i.e the average
  height of the 100 tallest trees (in cm).

- `stand_dominantTreeSpecies`: Dominant tree species name, determined in
  terms of basal area (and considering all tree sizes).

- `stand_hartBeckingIndex`: Hart-Becking index.

- `stand_foliarBiomass`: Standing biomass of leaves (in kg/m2).

- `stand_fuel`: Stand fine fuel load (in kg/m2).

- `stand_LAI`: Stand leaf area index (m2/m2).

- `stand_shrubVolume`: Stand shrub phytovolume (m3/m2).

## See also

[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`plant_basalArea`](https://emf-creaf.github.io/medfate/reference/plant_values.md),
[`summary.forest`](https://emf-creaf.github.io/medfate/reference/summary.forest.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Default species parameterization
data(SpParamsMED)
  
#Load example plot
data(exampleforest)
    
#A short way to obtain total basal area
stand_basalArea(exampleforest)
#> [1] 25.0333
    
```
