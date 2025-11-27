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
#> Initial plant water content (mm): 4.73001
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 4.72839
#> Final soil water content (mm): 274.93
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00162134
#> Plant water balance result (mm): -0.00163359
#> Change in soil water content (mm): -15.9454
#> Soil water balance result (mm): -15.9454
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 0
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 370
#>   Infiltration (mm) 401 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 25  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 247
#>   Plant extraction from soil (mm) 247  Plant water balance (mm) -0 Hydraulic redistribution (mm) 3
#>   Runoff (mm) 21 Deep drainage (mm) 131

#Queries the tables in 'Soil'
names(S1$Soil)
#> [1] "SWC"            "RWC"            "REW"            "ML"            
#> [5] "Psi"            "PlantExt"       "HydraulicInput"

#Monthly summary (averages) of soil relative water content
summary(S1, freq="months",FUN=mean, output="RWC")
#>                    1         2         3         4   Overall
#> 2001-01-01 0.9937121 0.9994108 0.9997816 1.0244415 1.0007791
#> 2001-02-01 0.9406420 0.9870166 0.9946293 0.9916723 0.9786150
#> 2001-03-01 0.9575952 0.9921181 1.0001649 1.0447712 0.9916070
#> 2001-04-01 0.8691697 0.9568931 0.9798945 0.9704885 0.9437391
#> 2001-05-01 0.8954893 0.9699323 0.9764303 0.9598696 0.9531287
#> 2001-06-01 0.6735878 0.8567177 0.9321367 0.9132172 0.8393223
#> 2001-07-01 0.9006140 0.9195416 0.8981414 0.8760523 0.9050023
#> 2001-08-01 0.9128670 0.9708158 0.9368690 0.9315899 0.9443077
#> 2001-09-01 0.9251964 0.9790296 0.9894022 0.9839707 0.9696283
#> 2001-10-01 0.9535367 0.9876684 0.9867276 0.9973065 0.9804278
#> 2001-11-01 0.9497432 0.9887678 1.0016620 1.0866380 0.9931741
#> 2001-12-01 0.8954726 0.9764038 0.9922932 0.9885181 0.9628353

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
#> 2001-01-01      0.004773662       0.003160351  0.006381252
#> 2001-02-01      0.006785896       0.004089571  0.007678385
#> 2001-03-01      0.006289085       0.003851156  0.007334555
#> 2001-04-01      0.011384560       0.006125331  0.010487139
#> 2001-05-01      0.010089242       0.005515647  0.009590412
#> 2001-06-01      0.055799936       0.022439959  0.029519026
#> 2001-07-01      0.018312570       0.008788778  0.014056669
#> 2001-08-01      0.008964251       0.005109164  0.009146389
#> 2001-09-01      0.007960909       0.004620249  0.008426610
#> 2001-10-01      0.006222902       0.003854981  0.007389951
#> 2001-11-01      0.006657749       0.004018635  0.007573598
#> 2001-12-01      0.009156016       0.005132602  0.009078819
# }
```
