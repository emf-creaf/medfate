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
#> Initial plant water content (mm): 4.67289
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 4.67119
#> Final soil water content (mm): 276.182
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00169593
#> Plant water balance result (mm): -0.00169593
#> Change in soil water content (mm): -14.693
#> Soil water balance result (mm): -14.693
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 79 Net rainfall (mm) 383
#>   Infiltration (mm) 413 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 26  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 235  Mistletoe transpiration (mm) 0
#>   Plant extraction from soil (mm) 235  Plant water balance (mm) -0 Hydraulic redistribution (mm) 1
#>   Runoff (mm) 21 Deep drainage (mm) 166

#Queries the tables in 'Soil'
names(S1$Soil)
#> [1] "SWC"            "RWC"            "REW"            "ML"            
#> [5] "Psi"            "PlantExt"       "HydraulicInput"

#Monthly summary (averages) of soil relative water content
summary(S1, freq="months",FUN=mean, output="RWC")
#>                    1         2         3         4   Overall
#> 2001-01-01 0.9942564 0.9996099 0.9998356 1.0272478 1.0012929
#> 2001-02-01 0.9426428 0.9906869 0.9957033 0.9928530 0.9809585
#> 2001-03-01 0.9589266 0.9947910 1.0017022 1.0558848 0.9945490
#> 2001-04-01 0.8741375 0.9725389 0.9856468 0.9770097 0.9533553
#> 2001-05-01 0.8992756 0.9794703 0.9924614 0.9881092 0.9649577
#> 2001-06-01 0.6729591 0.8903824 0.9611528 0.9448683 0.8635292
#> 2001-07-01 0.9033783 0.9487505 0.9712450 0.9525160 0.9443666
#> 2001-08-01 0.9225763 0.9845985 0.9931668 0.9956903 0.9734488
#> 2001-09-01 0.9265513 0.9856813 0.9929336 0.9877134 0.9739230
#> 2001-10-01 0.9564556 0.9938994 0.9965750 1.0145129 0.9879675
#> 2001-11-01 0.9495456 0.9924047 1.0089759 1.0958715 0.9974528
#> 2001-12-01 0.8935899 0.9817154 0.9941632 0.9904044 0.9652077

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
#> 2001-01-01      0.004774718       0.004636100  0.006369228
#> 2001-02-01      0.006904397       0.005927707  0.007545434
#> 2001-03-01      0.006409395       0.005602524  0.007232488
#> 2001-04-01      0.012499747       0.008542719  0.009432663
#> 2001-05-01      0.011048040       0.007780924  0.008844053
#> 2001-06-01      0.121924036       0.027267909  0.018126328
#> 2001-07-01      0.034436494       0.010339992  0.009829380
#> 2001-08-01      0.008815849       0.006853804  0.008221264
#> 2001-09-01      0.008362553       0.006680760  0.008134317
#> 2001-10-01      0.006201583       0.005522331  0.007196549
#> 2001-11-01      0.006887895       0.005870680  0.007464402
#> 2001-12-01      0.009964233       0.007569394  0.008874725
# }
```
