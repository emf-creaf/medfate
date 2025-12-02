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
#> Final plant water content (mm): 4.69659
#> Final soil water content (mm): 275.04
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00193896
#> Plant water balance result (mm): -0.00196771
#> Change in soil water content (mm): -15.8347
#> Soil water balance result (mm): -15.8347
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 370
#>   Infiltration (mm) 402 Infiltration excess (mm) 20 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 24  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 245
#>   Plant extraction from soil (mm) 245  Plant water balance (mm) -0 Hydraulic redistribution (mm) 3
#>   Runoff (mm) 20 Deep drainage (mm) 136

#Queries the tables in 'Soil'
names(S1$Soil)
#> [1] "SWC"            "RWC"            "REW"            "ML"            
#> [5] "Psi"            "PlantExt"       "HydraulicInput"

#Monthly summary (averages) of soil relative water content
summary(S1, freq="months",FUN=mean, output="RWC")
#>                    1         2         3         4   Overall
#> 2001-01-01 0.9933229 0.9995287 0.9998170 1.0245305 1.0007536
#> 2001-02-01 0.9356686 0.9893121 0.9953523 0.9923279 0.9786265
#> 2001-03-01 0.9536151 0.9937073 1.0006835 1.0457591 0.9915495
#> 2001-04-01 0.8486464 0.9652030 0.9840110 0.9744356 0.9437447
#> 2001-05-01 0.8866067 0.9741910 0.9776506 0.9647627 0.9535823
#> 2001-06-01 0.6340446 0.8731299 0.9442963 0.9259298 0.8411272
#> 2001-07-01 0.8801459 0.9300465 0.9214577 0.8980834 0.9127849
#> 2001-08-01 0.9060927 0.9765833 0.9556000 0.9452654 0.9513354
#> 2001-09-01 0.9159341 0.9834925 0.9915814 0.9859966 0.9700258
#> 2001-10-01 0.9495742 0.9905925 0.9895699 1.0014769 0.9818463
#> 2001-11-01 0.9445896 0.9913419 1.0034212 1.0885972 0.9936607
#> 2001-12-01 0.8872943 0.9803752 0.9937536 0.9899409 0.9630398

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
#> 2001-01-01      0.004802950       0.003164506  0.006385834
#> 2001-02-01      0.007239412       0.004179872  0.007769254
#> 2001-03-01      0.006670538       0.003929707  0.007416149
#> 2001-04-01      0.013859818       0.006761170  0.011142126
#> 2001-05-01      0.011572852       0.005847889  0.009924219
#> 2001-06-01      0.090837952       0.031064334  0.037361358
#> 2001-07-01      0.028534409       0.011392976  0.016305436
#> 2001-08-01      0.009790152       0.005243674  0.009227532
#> 2001-09-01      0.008852733       0.004824900  0.008629474
#> 2001-10-01      0.006527921       0.003908823  0.007428933
#> 2001-11-01      0.007135652       0.004116428  0.007668880
#> 2001-12-01      0.010149324       0.005330473  0.009277930
# }
```
