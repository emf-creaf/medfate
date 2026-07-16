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
#> Initial plant water content (mm): 4.20666
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:............
#> 
#> Final plant water content (mm): 4.20517
#> Final soil water content (mm): 276.623
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.0014939
#> Plant water balance result (mm): -0.0014939
#> Change in soil water content (mm): -14.2514
#> Soil water balance result (mm): -14.2514
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 0
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 76 Net rainfall (mm) 386
#>   Infiltration (mm) 415 Infiltration excess (mm) 22 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 27  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 226  Mistletoe transpiration (mm) 0
#>   Plant extraction from soil (mm) 226  Plant water balance (mm) -0 Hydraulic redistribution (mm) 1
#>   Runoff (mm) 22 Deep drainage (mm) 176

#Queries the tables in 'Soil'
names(S1$Soil)
#> [1] "SWC"            "RWC"            "REW"            "ML"            
#> [5] "Psi"            "PlantExt"       "HydraulicInput"

#Monthly summary (averages) of soil relative water content
summary(S1, freq="months",FUN=mean, output="RWC")
#>                    1         2         3         4   Overall
#> 2001-01-01 0.9945457 0.9996149 0.9998392 1.0281651 1.0014592
#> 2001-02-01 0.9453621 0.9907991 0.9957983 0.9930637 0.9816875
#> 2001-03-01 0.9609838 0.9948783 1.0019575 1.0590769 0.9954652
#> 2001-04-01 0.8844007 0.9736230 0.9860373 0.9779399 0.9563941
#> 2001-05-01 0.9037000 0.9797443 0.9926215 0.9887058 0.9662084
#> 2001-06-01 0.6857383 0.8935220 0.9621365 0.9463854 0.8681976
#> 2001-07-01 0.9098364 0.9521706 0.9734077 0.9584585 0.9484340
#> 2001-08-01 0.9273797 0.9849498 0.9933994 0.9975872 0.9749737
#> 2001-09-01 0.9314025 0.9860404 0.9932042 0.9893668 0.9754467
#> 2001-10-01 0.9594396 0.9943663 0.9974132 1.0186764 0.9895061
#> 2001-11-01 0.9519208 0.9925199 1.0106735 1.0982826 0.9987489
#> 2001-12-01 0.8967461 0.9818268 0.9943001 0.9906832 0.9660568

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
#> 2001-01-01      0.004766131       0.007453863  0.006364386
#> 2001-02-01      0.006758337       0.009409180  0.007483040
#> 2001-03-01      0.006292729       0.008911940  0.007182275
#> 2001-04-01      0.011513809       0.013080402  0.009192987
#> 2001-05-01      0.010537544       0.012169904  0.008721668
#> 2001-06-01      0.110666204       0.040593515  0.017511721
#> 2001-07-01      0.031461048       0.015551968  0.009544880
#> 2001-08-01      0.008449392       0.010740038  0.008107148
#> 2001-09-01      0.008036716       0.010480940  0.008017727
#> 2001-10-01      0.006066233       0.008755170  0.007122774
#> 2001-11-01      0.006744639       0.009319030  0.007403308
#> 2001-12-01      0.009727859       0.011981015  0.008795514
# }
```
