# Plant regeneration

Functions to simulate annual plant regeneration from seed recruitment or
from resprouting

## Usage

``` r
regeneration_seedproduction(forest, SpParams, control)

regeneration_seedrefill(seedBank, refillSpecies, refillPercent = NULL)

regeneration_seedmortality(seedBank, SpParams, minPercent = 1)

regeneration_germination(forest, SpParams, control)

regeneration_seedlings_daily(
  forest,
  SpParams,
  control,
  growthResult,
  verbose = FALSE
)

regeneration_seedlings(
  forest,
  SpParams,
  control,
  minMonthTemp,
  moistureIndex,
  verbose = FALSE
)

regeneration_recruitment(
  forest,
  SpParams,
  control,
  minMonthTemp,
  moistureIndex,
  verbose = FALSE
)

regeneration_recruitment_daily(
  forest,
  SpParams,
  control,
  growthResult,
  verbose = FALSE
)

regeneration_resprouting(
  forest,
  internalMortality,
  SpParams,
  control,
  management_results = NULL
)
```

## Arguments

- forest:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  and
  [`SpParamsDefinition`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- control:

  A list with default control parameters (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- seedBank:

  A data frame with columns 'Species' and 'Percent', describing a seed
  bank.

- refillSpecies:

  A string vector of species names corresponding to seed rain to refill
  seed bank.

- refillPercent:

  A numeric vector of indicating the percentage of seed bank refilling
  (if missing then seed bank is set to 100%).

- minPercent:

  A minimum percent of seed bank to retain entry in `seedBank` element
  of `forest`.

- growthResult:

  An object of class 'growth'.

- verbose:

  Boolean flag to indicate console output during calculations.

- minMonthTemp:

  Minimum month temperature.

- moistureIndex:

  Moisture index (annual precipitation over annual potential
  evapotranspiration).

- internalMortality:

  A data frame with mortality occurred in the last year of simulation.

- management_results:

  The result of calling a management function (see
  [`defaultManagementFunction`](https://emf-creaf.github.io/medfate/reference/defaultManagementFunction.md)).

## Value

- `regeneration_seedproduction` returns a list of species names

- `regeneration_seedrefill` and `regeneration_seedmortality` return a
  copy of the input `data.frame` object with an update seed bank.

- `regeneration_resprouting` and `regeneration_recruitment` return a new
  object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
  with the new plant cohorts.

## Details

- `regeneration_seedproduction` evaluates if reproductive individuals
  (i.e. sufficiently tall individuals) are present.

- `regeneration_seedrefill` fills seed bank of input `forest` object
  with seed rain.

- `regeneration_seedmortality` updates seed bank of input `forest`
  object according to annual seed mortality.

- `regeneration_recruitment` evaluates recruitment from the seed bank
  (or local seed production if seed bank is missing). Minimum month
  temperature and moisture index values are used to determine if
  recruitment was successful. Species also require a minimum amount of
  light at the ground level.

- `regeneration_resprouting` evaluates resprouting occurs after
  “mortality” from die-back (including drought- or pathogen-induced
  dessication), cutting or burning of the aerial part in a species with
  resprouting ability, but not after carbon starvation or baseline
  mortality (unspecific mortality causes).

## See also

[`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control <- defaultControl("Granier")
control$recruitmentMode = "annual/deterministic" 

#Recruitment limits
plant_parameter(exampleforest, SpParamsMED, "MinTempRecr")
#>    T1_148    T2_168    S1_165 
#>  1.083300 -3.744526  1.669536 
plant_parameter(exampleforest, SpParamsMED, "MinMoistureRecr")
#>     T1_148     T2_168     S1_165 
#> 0.10154153 0.09657161 0.22301894 

#Compare seed recruitment outcomes
regeneration_recruitment(exampleforest, SpParamsMED, control, 0, 0.25)
#> $treeData
#> [1] Species DBH     Height  N       Z50     Z95     Z100    Age     ObsID  
#> <0 rows> (or 0-length row.names)
#> 
#> $shrubData
#> [1] Species Height  Cover   Z50     Z95     Z100    Age     ObsID  
#> <0 rows> (or 0-length row.names)
#> 
#> $seedlingBank
#> [1] Species Percent Age     Z50     Z95     Z100   
#> <0 rows> (or 0-length row.names)
#> 
#> attr(,"class")
#> [1] "forest" "list"  
regeneration_recruitment(exampleforest, SpParamsMED, control, 3, 0.25)
#> $treeData
#> [1] Species DBH     Height  N       Z50     Z95     Z100    Age     ObsID  
#> <0 rows> (or 0-length row.names)
#> 
#> $shrubData
#> [1] Species Height  Cover   Z50     Z95     Z100    Age     ObsID  
#> <0 rows> (or 0-length row.names)
#> 
#> $seedlingBank
#> [1] Species Percent Age     Z50     Z95     Z100   
#> <0 rows> (or 0-length row.names)
#> 
#> attr(,"class")
#> [1] "forest" "list"  
```
