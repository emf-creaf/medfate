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
#> Final plant water content (mm): 4.69674
#> Final soil water content (mm): 275.757
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00178912
#> Plant water balance result (mm): -0.00180604
#> Change in soil water content (mm): -15.1184
#> Soil water balance result (mm): -15.1184
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 83 Net rainfall (mm) 380
#>   Infiltration (mm) 410 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 25  Herbaceous transpiration (mm) 0 Woody plant transpiration (mm) 246
#>   Plant extraction from soil (mm) 246  Plant water balance (mm) -0 Hydraulic redistribution (mm) 2
#>   Runoff (mm) 21 Deep drainage (mm) 154

#Queries the tables in 'Soil'
names(S1$Soil)
#> [1] "SWC"            "RWC"            "REW"            "ML"            
#> [5] "Psi"            "PlantExt"       "HydraulicInput"

#Monthly summary (averages) of soil relative water content
summary(S1, freq="months",FUN=mean, output="RWC")
#>                    1         2         3         4   Overall
#> 2001-01-01 0.9939327 0.9995597 0.9998246 1.0264658 1.0011126
#> 2001-02-01 0.9403668 0.9897262 0.9954148 0.9924381 0.9799215
#> 2001-03-01 0.9571712 0.9941056 1.0012242 1.0515212 0.9932837
#> 2001-04-01 0.8665224 0.9682320 0.9844315 0.9751415 0.9493322
#> 2001-05-01 0.8961655 0.9778100 0.9919579 0.9873223 0.9633498
#> 2001-06-01 0.6586943 0.8832971 0.9578311 0.9409453 0.8560707
#> 2001-07-01 0.8965367 0.9398255 0.9621061 0.9253273 0.9339717
#> 2001-08-01 0.9193542 0.9823236 0.9911082 0.9889506 0.9705431
#> 2001-09-01 0.9237217 0.9843918 0.9924843 0.9863422 0.9724828
#> 2001-10-01 0.9544451 0.9930806 0.9957817 1.0112287 0.9866193
#> 2001-11-01 0.9483943 0.9917490 1.0070386 1.0931390 0.9961309
#> 2001-12-01 0.8924066 0.9807768 0.9938735 0.9900085 0.9644373

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
#> 2001-01-01      0.004784526       0.003156632  0.006375245
#> 2001-02-01      0.006985942       0.004081034  0.007642733
#> 2001-03-01      0.006470491       0.003850444  0.007313595
#> 2001-04-01      0.012155169       0.006135056  0.010372283
#> 2001-05-01      0.010542298       0.005449125  0.009410390
#> 2001-06-01      0.070664646       0.025113737  0.031500404
#> 2001-07-01      0.021951608       0.009283859  0.014002860
#> 2001-08-01      0.008735273       0.004782823  0.008585625
#> 2001-09-01      0.008335547       0.004625986  0.008376971
#> 2001-10-01      0.006289725       0.003797446  0.007260573
#> 2001-11-01      0.006906538       0.004024817  0.007549637
#> 2001-12-01      0.009766188       0.005187575  0.009101858
# }
```
