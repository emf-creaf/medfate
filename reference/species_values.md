# Species description functions

Functions to calculate attributes of a
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object by species or to extract species parameters from a species
parameter table
([`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

## Usage

``` r
species_basalArea(x, SpParams)

species_cover(x, SpParams)

species_density(x, SpParams)

species_foliarBiomass(x, SpParams, gdd = NA_real_)

species_fuelLoading(x, SpParams, gdd = NA_real_, includeDead = TRUE)

species_LAI(x, SpParams, gdd = NA_real_, bounded = TRUE)

species_characterParameter(species, SpParams, parName)

species_parameter(
  species,
  SpParams,
  parName,
  fillMissing = TRUE,
  fillWithGenus = TRUE
)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

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

- species:

  A character vector of species names.

- parName:

  A string with a parameter name.

- fillMissing:

  A boolean flag to try imputation on missing values.

- fillWithGenus:

  A boolean flag to try imputation of missing values using genus values.

## Value

A vector with values for each species in `SpParams`:

- `species_basalArea`: Species basal area (m2/ha).

- `species_cover`: Shrub cover (in percent).

- `species_density`: Plant density (ind/ha). Tree density is directly
  taken from the forest object, while the shrub density is estimated
  from cover and height by calculating the area of a single individual.

- `species_foliarBiomass`: Standing biomass of leaves (in kg/m2).

- `species_fuel`: Fine fuel load (in kg/m2).

- `species_LAI`: Leaf area index (m2/m2).

- `species_phytovolume`: Shrub phytovolume (m3/m2).

- `species_parameter`: A numeric vector with the parameter values of
  each input species.

- `species_characterParameter`: A character vector with the parameter
  values of each input species.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`plant_basalArea`](https://emf-creaf.github.io/medfate/reference/plant_values.md),
[`summary.forest`](https://emf-creaf.github.io/medfate/reference/summary.forest.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
# Default species parameterization
data(SpParamsMED)

# Load example plot
data(exampleforest)

# Species basal area in the forest plot
species_basalArea(exampleforest, SpParamsMED)
#>  Pinus halepensis      Quercus ilex Quercus coccifera 
#>         18.604547          6.428755          0.000000 
  
# Value of parameter "Psi_Extract" for two species
species_parameter(c("Pinus halepensis", "Quercus ilex"), SpParamsMED, "Psi_Extract")
#> [1] -0.9218219 -1.9726871
    
```
