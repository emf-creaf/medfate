# Forest dynamics

Function `fordyn` implements a forest dynamics model that simulates
growth, mortality, recruitment and (optionally) management actions in a
given forest stand during a period specified in the input climatic data.

## Usage

``` r
fordyn(
  forest,
  soil,
  SpParams,
  meteo,
  control,
  latitude,
  elevation = NA,
  slope = NA,
  aspect = NA,
  CO2ByYear = numeric(0),
  management_function = NULL,
  management_args = NULL
)
```

## Arguments

- forest:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).
  Alternatively, the output of a previous run (an object of class
  `fordyn`), if continuing a previous simulation.

- soil:

  An object of class
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md).

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  and
  [`SpParamsDefinition`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- meteo:

  A data frame with daily weather data series (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- control:

  A list with default control parameters (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- latitude:

  Latitude (in degrees).

- elevation, slope, aspect:

  Elevation above sea level (in m), slope (in degrees) and aspect (in
  degrees from North).

- CO2ByYear:

  A named numeric vector with years as names and atmospheric CO2
  concentration (in ppm) as values. Used to specify annual changes in
  CO2 concentration along the simulation (as an alternative to
  specifying daily values in `meteo`).

- management_function:

  A function that implements forest management actions (see details).

- management_args:

  A list of additional arguments to be passed to the
  `management_function`.

## Value

A list of class `fordyn` with the following elements:

- `"StandSummary"`: A data frame with stand-level summaries (tree basal
  area, tree density, shrub cover, etc.) at the beginning of the
  simulation and after each simulated year.

- `"SpeciesSummary"`: A data frame with species-level summaries (tree
  basal area, tree density, shrub cover, etc.) at the beginning of the
  simulation and after each simulated year.

- `"CohortSummary"`: A data frame with cohort-level summaries (tree
  basal area, tree density, shrub cover, etc.) at the beginning of the
  simulation and after each simulated year.

- `"TreeTable"`: A data frame with tree-cohort data (species, density,
  diameter, height, etc.) at the beginning of the simulation (if any)
  and after each simulated year.

- `"DeadTreeTable"`: A data frame with dead tree-cohort data (species,
  density, diameter, height, etc.) at the beginning of the simulation
  and after each simulated year.

- `"CutTreeTable"`: A data frame with cut tree data (species, density,
  diameter, height, etc.) after each simulated year.

- `"ShrubTable"`: A data frame with shrub-cohort data (species, density,
  cover, height, etc.) at the beginning of the simulation and after each
  simulated year.

- `"DeadShrubTable"`: A data frame with dead shrub-cohort data (species,
  density, cover, height, etc.) at the beginning of the simulation (if
  any) and after each simulated year.

- `"CutShrubTable"`: A data frame with cut shrub data (species, density,
  cover, height, etc.) after each simulated year.

- `"ForestStructures"`: A list with the
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
  object of the stand at the beginning of the simulation and after each
  simulated year.

- `"GrowthResults"`: A list with the results of calling function
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)
  for each simulated year.

- `"ManagementArgs"`: A list of management arguments to be used in
  another call to `fordyn`.

- `"NextInputObject"`: An object of class `growthInput` to be used in a
  subsequent simulation.

- `"NextForestObject"`: An object of class `forest` to be used in a
  subsequent simulation.

## Details

Function `fordyn` simulates forest dynamics for annual time steps,
building on other simulation functions. For each simulated year, the
function performs the following steps:

1.  Calls function
    [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)
    to simulate daily water/carbon balance, growth and mortality
    processes and update the forest object.

2.  If required, calls function `management_function`, using as
    parameters the forest object and `management_args`, which may result
    in a density reduction for existing plant cohorts and/or a set of
    new planted cohorts.

3.  Simulate natural recruitment (for species present in the stand or
    given in a seed rain input).

4.  Prepares the input of function
    [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)
    for the next annual time step.

5.  Store forest status, management arguments, and summaries.

To enable forest management, the user needs to provide a function that
implements it, which is passed to `fordyn` via its argument
`management_function`. Such function should have the following
arguments:

- `"x"`: the
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
  object representing the stand to be managed.

- `"args"`: a list of parameters regulating the behavior of the
  management function.

- `"verbose"`: a logical flag to enable console output during the
  execution of the management function.

and return a list with the following elements:

- `"action"`: A string identifying the action performed (e.g.
  "thinning").

- `"N_tree_cut"`: A vector with the density of trees removed.

- `"Cover_shrub_cut"`: A vector with the cover of shrubs removed.

- `"planted_forest"`: An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
  with the new plant cohorts resulting from tree/shrub planting.

- `"management_args"`: A list of management arguments to be used in the
  next call to the management function.

An example of management function is provided in
[`defaultManagementFunction`](https://emf-creaf.github.io/medfate/reference/defaultManagementFunction.md).

## References

De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini
M, García-Valdés R, Nadal-Sala D, Sabaté S, Martin-StPaul N, Morin X,
D'Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A trait-enabled
model to simulate Mediterranean forest function and dynamics at regional
scales. Geoscientific Model Development 16: 3165-3201
(https://doi.org/10.5194/gmd-16-3165-2023).

## See also

[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
[`regeneration`](https://emf-creaf.github.io/medfate/reference/regeneration.md),
[`plot.growth`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md),
[`defaultManagementFunction`](https://emf-creaf.github.io/medfate/reference/defaultManagementFunction.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
# \donttest{
#Load example daily meteorological data
data(examplemeteo)
#Prepare a two-year meteorological data with half precipitation during 
#the second year
meteo2001 <- examplemeteo
meteo2002 <- examplemeteo
meteo2002$Precipitation <- meteo2002$Precipitation/2
meteo2002$dates <- seq(as.Date("2002-01-01"), 
                           as.Date("2002-12-31"), by="day")
meteo_01_02 <- rbind(meteo2001, meteo2002)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control <- defaultControl("Granier")

#Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

#Call simulation function
fd<-fordyn(exampleforest, examplesoil, 
           SpParamsMED, meteo_01_02, control,
           latitude = 41.82592, elevation = 100)
#> Simulating year 2001 (1/2):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
#> Simulating year 2002 (2/2):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1

#Stand-level summaries
fd$StandSummary
#>   Step NumTreeSpecies NumTreeCohorts NumShrubSpecies NumShrubCohorts
#> 1    0              2              2               1               1
#> 2    1              2              2               1               1
#> 3    2              2              2               1               1
#>   TreeDensityLive TreeBasalAreaLive DominantTreeHeight DominantTreeDiameter
#> 1        552.0000          25.03330           800.0000             37.55000
#> 2        551.3681          25.09418           804.2767             37.62962
#> 3        550.7349          25.15483           808.5257             37.70902
#>   QuadraticMeanTreeDiameter HartBeckingIndex ShrubCoverLive BasalAreaDead
#> 1                  24.02949         53.20353       3.750000    0.00000000
#> 2                  24.07248         52.95094       3.052800    0.03892163
#> 3                  24.11540         52.70294       3.067691    0.03916609
#>   ShrubCoverDead BasalAreaCut ShrubCoverCut
#> 1    0.000000000            0             0
#> 2    0.005284283            0             0
#> 3    0.004705828            0             0

#Tree table by annual steps
fd$TreeTable
#>   Step Year Cohort          Species      DBH   Height        N Z50  Z95 Z100
#> 1    0   NA T1_148 Pinus halepensis 37.55000 800.0000 168.0000 100  300   NA
#> 2    0   NA T2_168     Quercus ilex 14.60000 660.0000 384.0000 300 1000   NA
#> 3    1 2001 T1_148 Pinus halepensis 37.62962 804.2767 167.7002 100  300   NA
#> 4    1 2001 T2_168     Quercus ilex 14.62362 660.7351 383.6680 300 1000   NA
#> 5    2 2002 T1_148 Pinus halepensis 37.70902 808.5257 167.3997 100  300   NA
#> 6    2 2002 T2_168     Quercus ilex 14.64746 661.4781 383.3352 300 1000   NA
#>   Age ObsID
#> 1  NA  <NA>
#> 2  NA  <NA>
#> 3  NA    NA
#> 4  NA    NA
#> 5  NA    NA
#> 6  NA    NA

#Dead tree table by annual steps
fd$DeadTreeTable
#>   Step Year Cohort          Species      DBH   Height         N N_starvation
#> 1    1 2001 T1_148 Pinus halepensis 37.62962 804.2767 0.2998357            0
#> 2    1 2001 T2_168     Quercus ilex 14.62362 660.7351 0.3320154            0
#> 3    2 2002 T1_148 Pinus halepensis 37.70902 808.5257 0.3004794            0
#> 4    2 2002 T2_168     Quercus ilex 14.64746 661.4781 0.3328175            0
#>   N_dessication N_burnt N_resprouting_stumps Z50  Z95 Z100 Age ObsID
#> 1             0       0                    0 100  300   NA  NA    NA
#> 2             0       0                    0 300 1000   NA  NA    NA
#> 3             0       0                    0 100  300   NA  NA    NA
#> 4             0       0                    0 300 1000   NA  NA    NA
# }
```
