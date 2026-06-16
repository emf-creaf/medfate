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
#> Final plant water content (mm): 4.69657
#> Final soil water content (mm): 275.597
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.00195716
#> Plant water balance result (mm): -0.00195716
#> Change in soil water content (mm): -15.278
#> Soil water balance result (mm): -15.278
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 83 Net rainfall (mm) 379
#>   Infiltration (mm) 410 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 26  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 246  Mistletoe transpiration (mm) 0
#>   Plant extraction from soil (mm) 246  Plant water balance (mm) -0 Hydraulic redistribution (mm) 1
#>   Runoff (mm) 21 Deep drainage (mm) 153

#Queries the tables in 'Soil'
names(S1$Soil)
#> [1] "SWC"            "RWC"            "REW"            "ML"            
#> [5] "Psi"            "PlantExt"       "HydraulicInput"

#Monthly summary (averages) of soil relative water content
summary(S1, freq="months",FUN=mean, output="RWC")
#>                    1         2         3         4   Overall
#> 2001-01-01 0.9939254 0.9995584 0.9998243 1.0264495 1.0011086
#> 2001-02-01 0.9402986 0.9897019 0.9954072 0.9924270 0.9798926
#> 2001-03-01 0.9571058 0.9940808 1.0012082 1.0513934 0.9932410
#> 2001-04-01 0.8654939 0.9684146 0.9844284 0.9752993 0.9491800
#> 2001-05-01 0.8944921 0.9767795 0.9884316 0.9762043 0.9604679
#> 2001-06-01 0.6637610 0.8802067 0.9547155 0.9357566 0.8546670
#> 2001-07-01 0.8978875 0.9416216 0.9602171 0.9274586 0.9347387
#> 2001-08-01 0.9165364 0.9816983 0.9892815 0.9796820 0.9681907
#> 2001-09-01 0.9226619 0.9838782 0.9919070 0.9867313 0.9719187
#> 2001-10-01 0.9541977 0.9923357 0.9931837 1.0107775 0.9855389
#> 2001-11-01 0.9480521 0.9914508 1.0062545 1.0920598 0.9956144
#> 2001-12-01 0.8920457 0.9803679 0.9937509 0.9898497 0.9641402

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
#> 2001-02-01      0.007042305       0.004092019  0.007613252
#> 2001-03-01      0.006521342       0.003859944  0.007287650
#> 2001-04-01      0.013636609       0.006245863  0.009694499
#> 2001-05-01      0.011813046       0.005583891  0.009022613
#> 2001-06-01      0.144419334       0.022875553  0.019025585
#> 2001-07-01      0.039841455       0.008460900  0.010286129
#> 2001-08-01      0.009450930       0.004939353  0.008407302
#> 2001-09-01      0.008676338       0.004683011  0.008253989
#> 2001-10-01      0.006310566       0.003809713  0.007276993
#> 2001-11-01      0.007003698       0.004047575  0.007521592
#> 2001-12-01      0.010135072       0.005245079  0.008934919
# }
```
