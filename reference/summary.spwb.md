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
#> Initial plant water content (mm): 6.27649
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 6.27316
#> Final soil water content (mm): 273.017
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00332414
#> Plant water balance result (mm): -0.00332414
#> Change in soil water content (mm): -17.8575
#> Soil water balance result (mm): -17.8575
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 98 Net rainfall (mm) 364
#>   Infiltration (mm) 398 Infiltration excess (mm) 17 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 19  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 285  Mistletoe transpiration (mm) 0
#>   Plant extraction from soil (mm) 285  Plant water balance (mm) -0 Hydraulic redistribution (mm) 4
#>   Runoff (mm) 17 Deep drainage (mm) 111

#Queries the tables in 'Soil'
names(S1$Soil)
#> [1] "SWC"            "RWC"            "REW"            "ML"            
#> [5] "Psi"            "PlantExt"       "HydraulicInput"

#Monthly summary (averages) of soil relative water content
summary(S1, freq="months",FUN=mean, output="RWC")
#>                    1         2         3         4   Overall
#> 2001-01-01 0.9922841 0.9995369 0.9998012 1.0224347 1.0002910
#> 2001-02-01 0.9237283 0.9891183 0.9949362 0.9913448 0.9755395
#> 2001-03-01 0.9460333 0.9936564 0.9998809 1.0376451 0.9886978
#> 2001-04-01 0.8127859 0.9623238 0.9823793 0.9711539 0.9334184
#> 2001-05-01 0.8430044 0.9648168 0.9692545 0.9552007 0.9364210
#> 2001-06-01 0.5924615 0.8492615 0.9325044 0.9114991 0.8172349
#> 2001-07-01 0.8608872 0.8904710 0.9043581 0.8831814 0.8863944
#> 2001-08-01 0.8700562 0.9425732 0.9085965 0.9059218 0.9129111
#> 2001-09-01 0.8814436 0.9715974 0.9235005 0.9277876 0.9333788
#> 2001-10-01 0.9100092 0.9754789 0.9540769 0.9526840 0.9521865
#> 2001-11-01 0.9335624 0.9901299 0.9969592 1.0725127 0.9872319
#> 2001-12-01 0.8717185 0.9795592 0.9929063 0.9880597 0.9586453

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
#> 2001-01-01      0.004834264       0.007521984  0.006402859
#> 2001-02-01      0.008044764       0.010495855  0.007994937
#> 2001-03-01      0.007242272       0.009717724  0.007564517
#> 2001-04-01      0.021103273       0.018316846  0.010984798
#> 2001-05-01      0.019491910       0.016906520  0.010478273
#> 2001-06-01      0.249291506       0.069747516  0.024170174
#> 2001-07-01      0.060761231       0.026505968  0.013319477
#> 2001-08-01      0.014524272       0.016143190  0.010528102
#> 2001-09-01      0.012194086       0.013968830  0.009515885
#> 2001-10-01      0.009350481       0.011775533  0.008598226
#> 2001-11-01      0.007971032       0.010438269  0.007895301
#> 2001-12-01      0.011877493       0.013596654  0.009459432
# }
```
