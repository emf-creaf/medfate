# Woody plant cohort description functions

Functions to calculate attributes of woody plants in a
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object.

## Usage

``` r
plant_ID(x, SpParams, treeOffset = 0L, shrubOffset = 0L)

plant_basalArea(x, SpParams)

plant_largerTreeBasalArea(x, SpParams, self_proportion = 0.5)

plant_cover(x, SpParams)

plant_species(x, SpParams)

plant_speciesName(x, SpParams)

plant_density(x, SpParams)

plant_height(x, SpParams)

plant_individualArea(x, SpParams)

plant_crownRatio(x, SpParams)

plant_crownBaseHeight(x, SpParams)

plant_crownLength(x, SpParams)

plant_foliarBiomass(x, SpParams, gdd = NA_real_, competitionEffect = TRUE)

plant_fuelLoading(x, SpParams, gdd = NA_real_, includeDead = TRUE)

plant_equilibriumLeafLitter(x, SpParams, AET = 800)

plant_equilibriumSmallBranchLitter(
  x,
  SpParams,
  smallBranchDecompositionRate = 0.81
)

plant_phytovolume(x, SpParams)

plant_LAI(
  x,
  SpParams,
  gdd = NA_real_,
  bounded = TRUE,
  competitionEffect = TRUE
)

plant_characterParameter(x, SpParams, parName)

plant_parameter(x, SpParams, parName, fillMissing = TRUE, fillWithGenus = TRUE)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- treeOffset, shrubOffset:

  Integers to offset cohort IDs.

- self_proportion:

  Proportion of the target cohort included in the assessment

- gdd:

  Growth degree days (to account for leaf phenology effects).

- competitionEffect:

  Logical flag to indicate the inclusion of competition effect on LAI
  estimates.

- includeDead:

  A flag to indicate that standing dead fuels (dead branches) are
  included.

- AET:

  Actual annual evapotranspiration (in mm).

- smallBranchDecompositionRate:

  Decomposition rate of small branches.

- bounded:

  A boolean flag to indicate that extreme values should be prevented
  (maximum tree LAI = 7 and maximum shrub LAI = 3)

- parName:

  A string with a parameter name.

- fillMissing:

  A boolean flag to try imputation on missing values.

- fillWithGenus:

  A boolean flag to try imputation of missing values using genus values.

## Value

A vector with values for each woody plant cohort of the input
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object:

- `plant_basalArea`: Tree basal area (m2/ha).

- `plant_largerTreeBasalArea`: Basal area (m2/ha) of trees larger (in
  diameter) than the tree. Half of the trees of the same record are
  included.

- `plant_characterParameter`: The parameter values of each plant, as
  strings.

- `plant_cover`: Shrub cover (in percent).

- `plant_crownBaseHeight`: The height corresponding to the start of the
  crown (in cm).

- `plant_crownLength`: The difference between crown base height and
  total height (in cm).

- `plant_crownRatio`: The ratio between crown length and total height
  (between 0 and 1).

- `plant_density`: Plant density (ind/ha). Tree density is directly
  taken from the forest object, while the shrub density is estimated
  from cover and height by calculating the area of a single individual.

- `plant_equilibriumLeafLitter`: Litter biomass of leaves at equilibrium
  (in kg/m2).

- `plant_equilibriumSmallBranchLitter`: Litter biomass of small branches
  (\< 6.35 mm diameter) at equilibrium (in kg/m2).

- `plant_foliarBiomass`: Standing biomass of leaves (in kg/m2).

- `plant_fuelLoading`: Fine fuel load (in kg/m2).

- `plant_height`: Total height (in cm).

- `plant_ID`: Cohort coding for simulation functions (concatenation of
  'T' (Trees) or 'S' (Shrub), cohort index and species index).

- `plant_LAI`: Leaf area index (m2/m2).

- `plant_individualArea`: Area (m2) occupied by a shrub individual.

- `plant_parameter`: The parameter values of each plant, as numeric.

- `plant_phytovolume`: Shrub phytovolume (m3/m2).

- `plant_species`: Species identity integer (indices start with 0).

- `plant_speciesName`: String with species taxonomic name (or a
  functional group).

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`summary.forest`](https://emf-creaf.github.io/medfate/reference/summary.forest.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Default species parameterization
data(SpParamsMED)

#Load example plot
data(exampleforest)

#A plant-level way to obtain stand basal area
sum(plant_basalArea(exampleforest, SpParamsMED), na.rm=TRUE)
#> [1] 25.0333

#The analogous plant-level function for LAI
sum(plant_LAI(exampleforest, SpParamsMED))
#> [1] 1.584948
  
#The analogous plant-level function for fuel loading
sum(plant_fuelLoading(exampleforest, SpParamsMED))
#> [1] 0.5395798
      
#Summary function for 'forest' objects can be also used
summary(exampleforest, SpParamsMED)
#> Tree BA (m2/ha): 25.0333016  adult trees: 25.0333016  saplings: 0 
#> Density (ind/ha) adult trees: 552  saplings: 0  shrubs (estimated): 749.4923076 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 3.75  herbs: 10 
#> LAI (m2/m2) total: 1.7585845  adult trees: 1.5543216  saplings: 0  shrubs: 0.030626  herbs: 0.1736369 
#> Fuel loading (kg/m2) total: 0.5588728  adult trees: 0.5255004  saplings: 0  shrubs: 0.0140795  herbs: 0.019293 
#> PAR ground (%): 40.0075402  SWR ground (%): 50.7329667 

#Cohort IDs in the models
plant_ID(exampleforest, SpParamsMED)
#> [1] "T1_148" "T2_168" "S1_165"
      
```
