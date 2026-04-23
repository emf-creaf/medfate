# Summarize simulation results

Function `summary` summarizes the model's output in different temporal
steps (i.e. weekly, annual, ...).

## Usage

``` r
# S3 method for class 'spwb'
summary(
  object,
  freq = "years",
  output = "WaterBalance",
  FUN = sum,
  bySpecies = FALSE,
  months = NULL,
  ...
)

# S3 method for class 'aspwb'
summary(
  object,
  freq = "years",
  output = "WaterBalance",
  FUN = sum,
  bySpecies = FALSE,
  months = NULL,
  ...
)

# S3 method for class 'pwb'
summary(
  object,
  freq = "years",
  output = "WaterBalance",
  FUN = sum,
  bySpecies = FALSE,
  months = NULL,
  ...
)

# S3 method for class 'growth'
summary(
  object,
  freq = "years",
  output = "WaterBalance",
  FUN = sum,
  bySpecies = FALSE,
  months = NULL,
  ...
)

# S3 method for class 'fordyn'
summary(
  object,
  freq = "years",
  output = "WaterBalance",
  FUN = sum,
  bySpecies = FALSE,
  months = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `spwb`, `aspwb`, `pwb`, `growth` or `fordyn`.

- freq:

  Frequency of summary statistics (see
  [`cut.Date`](https://rdrr.io/r/base/cut.POSIXt.html)).

- output:

  The data table to be summarized. Accepted values are the path to data
  tables in `object`, such as 'WaterBalance', 'Soil', 'Stand' or
  'Plants\$LAI'. It is also possible to use strings like 'Transpiration'
  and the function will interpret it as 'Plants\$Transpiration'.

- FUN:

  The function to summarize results (e.g., `sum`, `mean`, ...)

- bySpecies:

  Allows aggregating output by species before calculating summaries
  (only has an effect with some values of `output`). Aggregation can
  involve a sum (as for plant lai or transpiration) or a LAI-weighted
  mean (as for plant stress or plant water potential).

- months:

  A vector of month numbers (1 to 12) to subset the season where
  summaries apply.

- ...:

  Additional parameters for function `summary`.

## Value

A matrix with dates as row names and the desired summaries in columns

## Note

When applied to
[`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
objects, the summary function can be used to gather the results of
different yearly steps into a single table while keeping a daily
resolution (i.e. using `freq = "days"`.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`pwb`](https://emf-creaf.github.io/medfate/reference/pwb.md),
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
[`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md),
[`plot.spwb`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md),
[`extract`](https://emf-creaf.github.io/medfate/reference/extract.md)

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
x <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
S1<-spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant water content (mm): 4.69853
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 4.69673
#> Final soil water content (mm): 275.597
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00180112
#> Plant water balance result (mm): -0.00180112
#> Change in soil water content (mm): -15.2779
#> Soil water balance result (mm): -15.2779
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 83 Net rainfall (mm) 379
#>   Infiltration (mm) 409 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 25  Herbaceous transpiration (mm) 0 Woody plant transpiration (mm) 248
#>   Plant extraction from soil (mm) 248  Plant water balance (mm) -0 Hydraulic redistribution (mm) 2
#>   Runoff (mm) 21 Deep drainage (mm) 152

#Queries the tables in 'Soil'
names(S1$Soil)
#> [1] "SWC"            "RWC"            "REW"            "ML"            
#> [5] "Psi"            "PlantExt"       "HydraulicInput"

#Monthly summary (averages) of soil relative water content
summary(S1, freq="months",FUN=mean, output="RWC")
#>                    1         2         3         4   Overall
#> 2001-01-01 0.9939254 0.9995584 0.9998243 1.0264495 1.0011086
#> 2001-02-01 0.9403087 0.9896970 0.9954058 0.9924257 0.9798925
#> 2001-03-01 0.9571119 0.9940778 1.0012073 1.0513909 0.9932407
#> 2001-04-01 0.8662313 0.9680515 0.9843775 0.9750693 0.9491700
#> 2001-05-01 0.8960120 0.9777122 0.9919273 0.9871979 0.9632537
#> 2001-06-01 0.6581725 0.8825682 0.9575542 0.9406417 0.8555524
#> 2001-07-01 0.8954626 0.9388493 0.9594759 0.9200050 0.9320891
#> 2001-08-01 0.9172623 0.9812995 0.9902977 0.9826257 0.9687722
#> 2001-09-01 0.9229054 0.9838555 0.9920716 0.9860795 0.9719416
#> 2001-10-01 0.9539247 0.9925547 0.9940312 1.0103994 0.9857441
#> 2001-11-01 0.9480897 0.9914491 1.0062875 1.0920833 0.9956336
#> 2001-12-01 0.8921194 0.9803405 0.9937419 0.9898401 0.9641432

#Queries the tables in 'Plants'
names(S1$Plants)
#>  [1] "LAI"                 "LAIlive"             "FPAR"               
#>  [4] "AbsorbedSWRFraction" "Transpiration"       "GrossPhotosynthesis"
#>  [7] "PlantPsi"            "LeafPLC"             "StemPLC"            
#> [10] "PlantWaterBalance"   "LeafRWC"             "StemRWC"            
#> [13] "LFMC"                "PlantStress"        

#Monthly summary (averages) of plant stress
summary(S1, freq="months",FUN=mean, output="PlantStress", 
        bySpecies = TRUE)
#>            Pinus halepensis Quercus coccifera Quercus ilex
#> 2001-01-01      0.004784746       0.003156729  0.006375382
#> 2001-02-01      0.006989065       0.004082373  0.007644662
#> 2001-03-01      0.006473653       0.003851819  0.007315582
#> 2001-04-01      0.012180853       0.006145362  0.010386440
#> 2001-05-01      0.010558457       0.005455470  0.009418937
#> 2001-06-01      0.071051366       0.025234052  0.031630680
#> 2001-07-01      0.022169351       0.009365977  0.014110460
#> 2001-08-01      0.008900049       0.004851299  0.008679593
#> 2001-09-01      0.008389474       0.004649397  0.008411343
#> 2001-10-01      0.006314736       0.003811202  0.007283814
#> 2001-11-01      0.006928037       0.004034940  0.007565379
#> 2001-12-01      0.009788701       0.005198114  0.009118803
# }
```
