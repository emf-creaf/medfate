# Basic water balance

## About this vignette

The present document describes how to run the soil plant water balance
model described in De Cáceres et al. (2015) using package `medfate`. The
document illustrates how to prepare the inputs, use the simulation
functions and inspect the outputs. All the details of the model design
and formulation can be found at the
[medfatebook](https://emf-creaf.github.io/medfatebook/index.html).
Because it introduces many basic features of simulations with package
`medfate`, this document should be read before addressing advanced
topics of water balance simulations or growth simulations.

## Preparing model inputs

Model inputs are explained in greater detail in vignettes
[*Understanding model
inputs*](https://emf-creaf.github.io/medfate/articles/intro/UnderstandingInputs.html)
and [*Preparing model
inputs*](https://emf-creaf.github.io/medfate/articles/intro/PreparingInputs.html).
Here we only review the different steps required to run function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md).

### Soil, vegetation, meteorology and species data

Soil information needs to be entered as a `data frame` with soil layers
in rows and physical attributes in columns. Soil physical attributes can
be initialized to default values, for a given number of layers, using
function
[`defaultSoilParams()`](https://emf-creaf.github.io/medfate/reference/defaultSoilParams.md):

``` r

examplesoil <- defaultSoilParams(4)
examplesoil
```

    ##   widths clay sand om nitrogen ph  bd rfc
    ## 1    300   25   25 NA       NA NA 1.5  25
    ## 2    700   25   25 NA       NA NA 1.5  45
    ## 3   1000   25   25 NA       NA NA 1.5  75
    ## 4   2000   25   25 NA       NA NA 1.5  95

As explained in the package overview, models included in `medfate` were
primarily designed to be ran on **forest inventory plots**. Here we use
the example object provided with the package:

``` r

data(exampleforest)
exampleforest
```

    ## $treeData
    ##            Species   DBH Height   N Z50  Z95
    ## 1 Pinus halepensis 37.55    800 168 100  300
    ## 2     Quercus ilex 14.60    660 384 300 1000
    ## 
    ## $shrubData
    ##             Species Height Cover Z50  Z95
    ## 1 Quercus coccifera     80  3.75 200 1000
    ## 
    ## attr(,"class")
    ## [1] "forest" "list"

Importantly, a data frame with daily weather for the period to be
simulated is required. Here we use the default data frame included with
the package:

``` r

data(examplemeteo)
head(examplemeteo)
```

    ##        dates MinTemperature MaxTemperature Precipitation MinRelativeHumidity
    ## 1 2001-01-01     -0.5934215       6.287950      4.869109            65.15411
    ## 2 2001-01-02     -2.3662458       4.569737      2.498292            57.43761
    ## 3 2001-01-03     -3.8541036       2.661951      0.000000            58.77432
    ## 4 2001-01-04     -1.8744860       3.097705      5.796973            66.84256
    ## 5 2001-01-05      0.3288287       7.551532      1.884401            62.97656
    ## 6 2001-01-06      0.5461322       7.186784     13.359801            74.25754
    ##   MaxRelativeHumidity Radiation WindSpeed
    ## 1           100.00000  12.89251  2.000000
    ## 2            94.71780  13.03079  7.662544
    ## 3            94.66823  16.90722  2.000000
    ## 4            95.80950  11.07275  2.000000
    ## 5           100.00000  13.45205  7.581347
    ## 6           100.00000  12.84841  6.570501

Finally, simulations in `medfate` require a data frame with species
parameter values, which we load using defaults for Catalonia (NE Spain):

``` r

data("SpParamsMED")
```

### Simulation control

Apart from data inputs, the behaviour of simulation models can be
controlled using a set of global parameters. The default
parameterization is obtained using function
[`defaultControl()`](https://emf-creaf.github.io/medfate/reference/defaultControl.md):

``` r

control <- defaultControl("Granier")
```

Some parameters deserve explanation here:

1.  Console output can be turned off by setting `verbose = FALSE`.
2.  The soil water retention curves can be switched between Saxton’s and
    Van Genuchten’s using parameter `soilFunctions`.
3.  The complexity of the soil water balance calculations will be very
    different if we set `transpirationMode = "Sperry"` or
    `transpirationMode = "Sureau"`, instead of
    `transpirationMode = "Granier"`.

### Water balance input object

A last object is needed before calling simulation functions, called
`spwbInput`. It consists in the compilation of aboveground and
belowground parameters and the specification of additional parameter
values for each plant cohort. The object can be generated using function
[`spwbInput()`](https://emf-creaf.github.io/medfate/reference/modelInput.md):

``` r

x <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
```

Different parameter variables will be drawn depending on the value of
`transpirationMode`. For the basic water balance model
(`transpirationMode = "Granier"`), relatively few parameters are needed.
All the input information for forest data and species parameter values
can be inspected by accessing the different elements of this object,
whose names are.

``` r

names(x)
```

    ##  [1] "control"                 "soil"                   
    ##  [3] "snowpack"                "canopy"                 
    ##  [5] "herbLAI"                 "herbLAImax"             
    ##  [7] "cohorts"                 "above"                  
    ##  [9] "below"                   "belowLayers"            
    ## [11] "paramsPhenology"         "paramsAnatomy"          
    ## [13] "paramsInterception"      "paramsTranspiration"    
    ## [15] "paramsWaterStorage"      "internalPhenology"      
    ## [17] "internalWater"           "internalLAIDistribution"
    ## [19] "internalFCCS"            "version"

Finally, note that users can set cohort-specific parameters for soil
water balance (instead of using species-level values) by modifying
manually the parameter values in this object. Since some parameters may
be coordinated by design, however, it is better to use specific package
functions for this purpose.

## Executing the soil water balance model

### Water balance for a single day

Soil water balance simulations will normally span periods of several
months or years, but since the model operates at a daily temporal scale,
it is possible to perform soil water balance for one day only. This is
done using function
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md).
In the following code we select day 100 from the meteorological input
data and perform soil water balance for that day only:

``` r

d <- 100
date <- examplemeteo$dates[d]
meteovec <- unlist(examplemeteo[d,])
sd1<-spwb_day(x, date, meteovec,  
             latitude = 41.82592, elevation = 100, slope= 0, aspect = 0)
```

Function
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
is most useful when working with the complex transpiration model. This
is why so many meteorological variables are required. The output of
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
is a list with five elements:

``` r

names(sd1)
```

    ## [1] "cohorts"      "topography"   "weather"      "WaterBalance" "Soil"        
    ## [6] "Stand"        "Plants"

- **cohorts**: Table with the species code and species name of each
  cohort.
- **WaterBalance**: Contains the soil water balance flows
  (precipitation, infiltration, transpiration, …)
- **Soil**: Contains output values by soil layer (i.e. water evaporated
  from each soil layer, water transpired from each soil layer and the
  final soil water potential).
- **Stand**: A list with stand LAI (expanded and dead leaves), canopy
  water retention capacity and the proportion of light (SWR or PAR)
  reaching the ground.
- **Plants**: Contains output values by plant cohort (i.e. LAI values,
  transpiration, water potential, drought stress index, …).

``` r

sd1
```

    ## $cohorts
    ##         SP              Name
    ## T1_148 148  Pinus halepensis
    ## T2_168 168      Quercus ilex
    ## S1_165 165 Quercus coccifera
    ## 
    ## $topography
    ## elevation     slope    aspect 
    ##       100         0         0 
    ## 
    ## $weather
    ##        tday        prec        tmin        tmax       rhmin       rhmax 
    ##   6.2323731   0.0000000   0.3881289  10.0320962  42.0207334  82.3036989 
    ##         rad        wind        Catm        Patm         pet        rint 
    ##  28.7201692   3.3228840 386.0000000         NaN   3.9023342   1.5000000 
    ## 
    ## $WaterBalance
    ##                     PET                    Rain                    Snow 
    ##                3.902334                0.000000                0.000000 
    ##                 NetRain                Snowmelt                   Runon 
    ##                0.000000                0.000000                0.000000 
    ##            Infiltration      InfiltrationExcess        SaturationExcess 
    ##                0.000000                0.000000                0.000000 
    ##                  Runoff            DeepDrainage         CapillarityRise 
    ##                0.000000                0.000000                0.000000 
    ##         SoilEvaporation       HerbTranspiration         PlantExtraction 
    ##                0.500000                0.000000                1.064905 
    ##           Transpiration  MistletoeTranspiration HydraulicRedistribution 
    ##                1.064905                0.000000                0.000000 
    ## 
    ## $Soil
    ##           Psi HerbTranspiration HydraulicInput HydraulicOutput PlantExtraction
    ## 1 -0.03631332                 0              0     0.876840481     0.876840481
    ## 2 -0.03323042                 0              0     0.171326921     0.171326921
    ## 3 -0.03302880                 0              0     0.013949411     0.013949411
    ## 4 -0.03301439                 0              0     0.002788117     0.002788117
    ## 
    ## $Stand
    ##          LAI      LAIherb      LAIlive  LAIexpanded      LAIdead LAImistletoe 
    ##     1.826041     0.000000     1.826041     1.826041     0.000000     0.000000 
    ##           Cm   LgroundPAR   LgroundSWR 
    ##     1.517698    38.912485    49.700665 
    ## 
    ## $Plants
    ##               LAI    LAIlive     FPAR AbsorbedSWRFraction Extraction
    ## T1_148 1.20935408 1.20935408 88.57701           44.432640  0.7497994
    ## T2_168 0.56790878 0.56790878 69.05317           24.578915  0.2770367
    ## S1_165 0.04877848 0.04877848 39.89406            3.156504  0.0380688
    ##        Transpiration MistletoeTranspiration GrossPhotosynthesis PlantPsi
    ## T1_148     0.7497994                      0           4.9357104   -0.033
    ## T2_168     0.2770367                      0           2.5171799   -0.033
    ## S1_165     0.0380688                      0           0.2240176   -0.033
    ##                DDS   StemRWC   LeafRWC      LFMC StemPLC      LeafPLC
    ## T1_148 0.004613739 0.9998346 0.9981628 122.02311       0 0.0000000000
    ## T2_168 0.006282261 0.9997405 0.9991208 109.14738       0 0.0009945146
    ## S1_165 0.007310922 0.9995007 0.9987189  98.29725       0 0.0000000000
    ##         WaterBalance
    ## T1_148 -8.881784e-16
    ## T2_168  0.000000e+00
    ## S1_165  0.000000e+00
    ## 
    ## attr(,"class")
    ## [1] "spwb_day" "list"

### Water balance for multiple days

Most often, users will use function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md) to run
the soil water balance model. This function requires the `spwbInput`
object and the meteorological data frame. However, function
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
by default modifies the state variables of the input objects. In
particular, the values of soil moisture are now:

``` r

x$soil$W
```

    ## [1] 0.9798040 0.9985313 0.9998158 0.9999080

We simply reset state variables to their default values so that new
simulations are not affected by the end state of the previous
simulation:

``` r

resetInputs(x)
x$soil$W
```

    ## [1] 1 1 1 1

Now we are ready to call function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md):

``` r

S <- spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
```

    ## Initial plant water content (mm): 6.27649
    ## Initial soil water content (mm): 290.875
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2001]:............
    ## 
    ## Final plant water content (mm): 6.27316
    ## Final soil water content (mm): 273.017
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.00332414
    ## Plant water balance result (mm): -0.00332414
    ## Change in soil water content (mm): -17.8575
    ## Soil water balance result (mm): -17.8575
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): -7.10543e-15
    ## Water balance components:
    ##   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
    ##   Interception (mm) 98 Net rainfall (mm) 364
    ##   Infiltration (mm) 398 Infiltration excess (mm) 17 Saturation excess (mm) 0 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 19  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 285  Mistletoe transpiration (mm) 0
    ##   Plant extraction from soil (mm) 285  Plant water balance (mm) -0 Hydraulic redistribution (mm) 4
    ##   Runoff (mm) 17 Deep drainage (mm) 111

Function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md)
returns an object of class with the same name, actually a list:

``` r

class(S)
```

    ## [1] "spwb" "list"

If we inspect its elements, we realize that the output is arranged
differently than in
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md):

``` r

names(S)
```

    ##  [1] "latitude"     "topography"   "weather"      "spwbInput"    "spwbOutput"  
    ##  [6] "WaterBalance" "Soil"         "Snow"         "Stand"        "Plants"

In particular, element `spwbInput` contains a copy of the input
parameters that were used to run the model:

``` r

names(S$spwbInput)
```

    ##  [1] "control"                 "soil"                   
    ##  [3] "snowpack"                "canopy"                 
    ##  [5] "herbLAI"                 "herbLAImax"             
    ##  [7] "cohorts"                 "above"                  
    ##  [9] "below"                   "belowLayers"            
    ## [11] "paramsPhenology"         "paramsAnatomy"          
    ## [13] "paramsInterception"      "paramsTranspiration"    
    ## [15] "paramsWaterStorage"      "internalPhenology"      
    ## [17] "internalWater"           "internalLAIDistribution"
    ## [19] "internalFCCS"            "version"

As before, `WaterBalance` contains water balance components, but in this
case in form of a data frame with days in rows:

``` r

head(S$WaterBalance)
```

    ##                  PET Precipitation      Rain Snow    NetRain Snowmelt
    ## 2001-01-01 0.8828475      4.869109  4.869109    0  3.2972315        0
    ## 2001-01-02 1.6375337      2.498292  2.498292    0  0.9721473        0
    ## 2001-01-03 1.3017026      0.000000  0.000000    0  0.0000000        0
    ## 2001-01-04 0.5690790      5.796973  5.796973    0  4.2353558        0
    ## 2001-01-05 1.6760567      1.884401  1.884401    0  0.7332674        0
    ## 2001-01-06 1.2077028     13.359801 13.359801    0 11.5940646        0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2001-01-01    3.2972315                  0                0      0    2.6175309
    ## 2001-01-02    0.9721473                  0                0      0    0.1881079
    ## 2001-01-03    0.0000000                  0                0      0    0.0000000
    ## 2001-01-04    4.2353558                  0                0      0    3.0529310
    ## 2001-01-05    0.7332674                  0                0      0    0.1206369
    ## 2001-01-06   11.5940646                  0                0      0    4.1214138
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2001-01-01               0          2.2515780     1.571877       0.4387811
    ## 2001-01-02               0          2.4730095     1.526144       0.5000000
    ## 2001-01-03               0          0.7182154     0.000000       0.3630184
    ## 2001-01-04               0          1.8630006     1.561617       0.1461396
    ## 2001-01-05               0          2.1085118     1.151134       0.5000000
    ## 2001-01-06               0          2.3578418     1.765736       0.2625835
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01                 0       0.2409196     0.2409196
    ## 2001-01-02                 0       0.4468653     0.4468653
    ## 2001-01-03                 0       0.3551970     0.3551970
    ## 2001-01-04                 0       0.1552438     0.1552438
    ## 2001-01-05                 0       0.4573778     0.4573778
    ## 2001-01-06                 0       0.3295223     0.3295223
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2001-01-01                      0                       0
    ## 2001-01-02                      0                       0
    ## 2001-01-03                      0                       0
    ## 2001-01-04                      0                       0
    ## 2001-01-05                      0                       0
    ## 2001-01-06                      0                       0

Element `Plants` is in turn a list with several dataframes with plant
output variables, for example plant water potentials are in:

``` r

head(S$Plants$PlantPsi)
```

    ##                 T1_148      T2_168      S1_165
    ## 2001-01-01 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-02 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-03 -0.03335674 -0.03318817 -0.03325548
    ## 2001-01-04 -0.03475656 -0.03396887 -0.03428379
    ## 2001-01-05 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-06 -0.03376019 -0.03340110 -0.03354466

## Inspecting model outputs

### Plots

Package `medfate` provides a simple `plot` function for objects of class
`spwb`. It can be used to show meteorological inputs, snow dynamics, and
different components of the water balance:

``` r

plot(S, type = "PET_Precipitation")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-19-1.png)

``` r

plot(S, type = "Snow")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-19-2.png)

``` r

plot(S, type = "Export")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-19-3.png)

``` r

plot(S, type = "Evapotranspiration")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-20-1.png)

Function `plot` is also allows displaying soil moisture dynamics by
layer, which can be done in four different ways (the first two only
imply a change in axis units):

``` r

plot(S, type="SoilTheta")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-21-1.png)

``` r

plot(S, type="SoilRWC")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-21-2.png)

``` r

plot(S, type="SoilPsi")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-21-3.png)

``` r

plot(S, type="SoilVol")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-21-4.png)

Finally, the same function can also be used to draw the dynamics of
plant variables by cohorts, such as transpiration, gross photosynthesis
or water potential:

``` r

plot(S, type="Transpiration")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-22-1.png)

``` r

plot(S, type="GrossPhotosynthesis")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-22-2.png)

``` r

plot(S, type="PlantPsi")
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-22-3.png)

Finally, one can interactively create plots using function `shinyplot`,
e.g.:

``` r

shinyplot(S)
```

### Extracting output

Simulation outputs in form of lists have a nested structure that is not
easy to handle. Functions are provided to extract model outputs as
`data.frame` objects. The following code extracts daily series of
stand-level variables, including their units:

``` r

extract(S, "forest", addunits = TRUE) |>
  tibble::as_tibble()
```

    ## # A tibble: 365 × 31
    ##    date           PET Precipitation    Rain   Snow NetRain Snowmelt Infiltration
    ##    <date>     [L/m^2]       [L/m^2] [L/m^2] [L/m^… [L/m^2]  [L/m^2]      [L/m^2]
    ##  1 2001-01-01   0.883          4.87    4.87   0      3.30      0           3.30 
    ##  2 2001-01-02   1.64           2.50    2.50   0      0.972     0           0.972
    ##  3 2001-01-03   1.30           0       0      0      0         0           0    
    ##  4 2001-01-04   0.569          5.80    5.80   0      4.24      0           4.24 
    ##  5 2001-01-05   1.68           1.88    1.88   0      0.733     0           0.733
    ##  6 2001-01-06   1.21          13.4    13.4    0     11.6       0          11.6  
    ##  7 2001-01-07   0.637          5.38    0      5.38   0         0           0    
    ##  8 2001-01-08   0.832          0       0      0      0         0           0    
    ##  9 2001-01-09   1.98           0       0      0      0         0           0    
    ## 10 2001-01-10   0.829          5.12    5.12   0      3.55      5.38        8.92 
    ## # ℹ 355 more rows
    ## # ℹ 23 more variables: InfiltrationExcess [L/m^2], SaturationExcess [L/m^2],
    ## #   Runoff [L/m^2], DeepDrainage [L/m^2], CapillarityRise [L/m^2],
    ## #   Evapotranspiration [L/m^2], Interception [L/m^2], SoilEvaporation [L/m^2],
    ## #   HerbTranspiration [L/m^2], PlantExtraction [L/m^2], Transpiration [L/m^2],
    ## #   MistletoeTranspiration [L/m^2], HydraulicRedistribution [L/m^2],
    ## #   LAI [m^2/m^2], LAIherb [m^2/m^2], LAIlive [m^2/m^2], …

And a similar code can be used to daily series of cohort-level
variables:

``` r

extract(S, "cohort", addunits = TRUE) |>
  tibble::as_tibble()
```

    ## # A tibble: 1,095 × 17
    ##    date       cohort species               LAI LAIlive FPAR AbsorbedSWRFraction
    ##    <date>     <chr>  <chr>            [m^2/m^… [m^2/m…  [%]               <dbl>
    ##  1 2001-01-01 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ##  2 2001-01-02 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ##  3 2001-01-03 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ##  4 2001-01-04 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ##  5 2001-01-05 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ##  6 2001-01-06 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ##  7 2001-01-07 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ##  8 2001-01-08 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ##  9 2001-01-09 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ## 10 2001-01-10 T1_148 Pinus halepensis     1.21    1.21 88.6                44.4
    ## # ℹ 1,085 more rows
    ## # ℹ 10 more variables: Transpiration [L/m^2], GrossPhotosynthesis [L/m^2],
    ## #   PlantPsi [MPa], LeafPLC <dbl>, StemPLC <dbl>, PlantWaterBalance [L/m^2],
    ## #   LeafRWC [%], StemRWC [%], LFMC [%], PlantStress <dbl>

### Temporal summaries

While the simulation model uses daily steps, users will normally be
interested in outputs at larger time scales. The package provides a
`summary` for objects of class `spwb`. This function can be used to
summarize the model’s output at different temporal steps (i.e. weekly,
annual, …). For example, to obtain the water balance by months one can
use:

``` r

summary(S, freq="months",FUN=mean, output="WaterBalance")
```

    ##                 PET Precipitation       Rain      Snow    NetRain   Snowmelt
    ## 2001-01-01 1.011397    2.41127383 1.87415609 0.5371177 1.31145572 0.42235503
    ## 2001-02-01 2.278646    0.17855109 0.08778069 0.0907704 0.03415765 0.19831578
    ## 2001-03-01 2.368035    2.41917349 2.41917349 0.0000000 1.91458541 0.01762496
    ## 2001-04-01 3.086567    0.63056064 0.29195973 0.3386009 0.12854028 0.33860091
    ## 2001-05-01 3.662604    0.76337345 0.76337345 0.0000000 0.57061353 0.00000000
    ## 2001-06-01 5.265359    0.21959509 0.21959509 0.0000000 0.15317183 0.00000000
    ## 2001-07-01 4.443053    3.27810591 3.27810591 0.0000000 2.78038699 0.00000000
    ## 2001-08-01 4.463242    1.92222891 1.92222891 0.0000000 1.52567092 0.00000000
    ## 2001-09-01 3.453891    1.30651303 1.30651303 0.0000000 1.03132466 0.00000000
    ## 2001-10-01 2.405506    1.33598175 1.33598175 0.0000000 1.03574616 0.00000000
    ## 2001-11-01 1.716591    2.20566281 1.47764599 0.7280168 1.32071471 0.72801682
    ## 2001-12-01 1.608082    0.05046181 0.05046181 0.0000000 0.01963594 0.00000000
    ##            Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2001-01-01   1.73381075         0.00000000                0 0.00000000
    ## 2001-02-01   0.23247343         0.00000000                0 0.00000000
    ## 2001-03-01   1.93221036         0.00000000                0 0.00000000
    ## 2001-04-01   0.46714119         0.00000000                0 0.00000000
    ## 2001-05-01   0.57061353         0.00000000                0 0.00000000
    ## 2001-06-01   0.15317183         0.00000000                0 0.00000000
    ## 2001-07-01   2.55343876         0.22694822                0 0.22694822
    ## 2001-08-01   1.49948275         0.02618817                0 0.02618817
    ## 2001-09-01   1.03132466         0.00000000                0 0.00000000
    ## 2001-10-01   0.95026994         0.08547622                0 0.08547622
    ## 2001-11-01   1.81733582         0.23139571                0 0.23139571
    ## 2001-12-01   0.01963594         0.00000000                0 0.00000000
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2001-01-01  1.382178687               0          0.9989094   0.56270036
    ## 2001-02-01  0.008777884               0          0.7121687   0.05362304
    ## 2001-03-01  1.041106893               0          1.2590795   0.50458808
    ## 2001-04-01  0.000000000               0          1.0045101   0.16341945
    ## 2001-05-01  0.000000000               0          1.2035039   0.19275992
    ## 2001-06-01  0.000000000               0          1.2353774   0.06642325
    ## 2001-07-01  0.000000000               0          1.7309780   0.49771893
    ## 2001-08-01  0.000000000               0          1.6433062   0.39655800
    ## 2001-09-01  0.000000000               0          1.2394083   0.27518837
    ## 2001-10-01  0.032325770               0          0.9986147   0.30023559
    ## 2001-11-01  1.162945478               0          0.6816138   0.15693127
    ## 2001-12-01  0.000000000               0          0.4836883   0.03082587
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01     0.160266514                 0       0.2759425     0.2759425
    ## 2001-02-01     0.038658527                 0       0.6198871     0.6198871
    ## 2001-03-01     0.109850392                 0       0.6446411     0.6446411
    ## 2001-04-01     0.010088380                 0       0.8310022     0.8310022
    ## 2001-05-01     0.026066728                 0       0.9846772     0.9846772
    ## 2001-06-01     0.004851102                 0       1.1641031     1.1641031
    ## 2001-07-01     0.078073983                 0       1.1551851     1.1551851
    ## 2001-08-01     0.039191568                 0       1.2075566     1.2075566
    ## 2001-09-01     0.027618004                 0       0.9366019     0.9366019
    ## 2001-10-01     0.044785147                 0       0.6535939     0.6535939
    ## 2001-11-01     0.057721755                 0       0.4669608     0.4669608
    ## 2001-12-01     0.016697111                 0       0.4361653     0.4361653
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2001-01-01                      0             0.001297018
    ## 2001-02-01                      0             0.000000000
    ## 2001-03-01                      0             0.001345449
    ## 2001-04-01                      0             0.000000000
    ## 2001-05-01                      0             0.005356283
    ## 2001-06-01                      0             0.000000000
    ## 2001-07-01                      0             0.056468643
    ## 2001-08-01                      0             0.040628214
    ## 2001-09-01                      0             0.027179247
    ## 2001-10-01                      0             0.008441485
    ## 2001-11-01                      0             0.002010522
    ## 2001-12-01                      0             0.000000000

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of plant cohorts by months:

``` r

summary(S, freq="months",FUN=mean, output="PlantStress")
```

    ##                 T1_148      T2_168      S1_165
    ## 2001-01-01 0.004834264 0.006402859 0.007521984
    ## 2001-02-01 0.008044764 0.007994937 0.010495855
    ## 2001-03-01 0.007242272 0.007564517 0.009717724
    ## 2001-04-01 0.021103273 0.010984798 0.018316846
    ## 2001-05-01 0.019491910 0.010478273 0.016906520
    ## 2001-06-01 0.249291506 0.024170174 0.069747516
    ## 2001-07-01 0.060761231 0.013319477 0.026505968
    ## 2001-08-01 0.014524272 0.010528102 0.016143190
    ## 2001-09-01 0.012194086 0.009515885 0.013968830
    ## 2001-10-01 0.009350481 0.008598226 0.011775533
    ## 2001-11-01 0.007971032 0.007895301 0.010438269
    ## 2001-12-01 0.011877493 0.009459432 0.013596654

The `summary` function can be also used to aggregate the output by
species. In this case, the values of plant cohorts belonging to the same
species will be averaged using LAI values as weights. For example, we
may average the daily drought stress across cohorts of the same species
(here there is only one cohort by species, so this does not modify the
output):

``` r

head(summary(S, freq="day", output="PlantStress", bySpecies = TRUE))
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01      0.004613739       0.007310922  0.006282261
    ## 2001-01-02      0.004613739       0.007310922  0.006282261
    ## 2001-01-03      0.004688806       0.007384316  0.006323310
    ## 2001-01-04      0.004987215       0.007681368  0.006493968
    ## 2001-01-05      0.004613739       0.007310922  0.006282261
    ## 2001-01-06      0.004774186       0.007467586  0.006369800

Or we can combine the aggregation by species with a temporal aggregation
(here monthly averages):

``` r

summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01      0.004834264       0.007521984  0.006402859
    ## 2001-02-01      0.008044764       0.010495855  0.007994937
    ## 2001-03-01      0.007242272       0.009717724  0.007564517
    ## 2001-04-01      0.021103273       0.018316846  0.010984798
    ## 2001-05-01      0.019491910       0.016906520  0.010478273
    ## 2001-06-01      0.249291506       0.069747516  0.024170174
    ## 2001-07-01      0.060761231       0.026505968  0.013319477
    ## 2001-08-01      0.014524272       0.016143190  0.010528102
    ## 2001-09-01      0.012194086       0.013968830  0.009515885
    ## 2001-10-01      0.009350481       0.011775533  0.008598226
    ## 2001-11-01      0.007971032       0.010438269  0.007895301
    ## 2001-12-01      0.011877493       0.013596654  0.009459432

### Specific output functions

The package provides some functions to extract or transform specific
outputs from soil plant water balance simulations. In particular,
function
[`droughtStress()`](https://emf-creaf.github.io/medfate/reference/droughtStress.md)
allows calculating several plant stress indices, such as the number of
days with drought stress \> 0.5 or the maximum drought stress:

``` r

droughtStress(S, index = "NDD", freq = "years", draw=FALSE)
```

    ##            T1_148 T2_168 S1_165
    ## 2001-01-01      0      0      0

``` r

droughtStress(S, index = "MDS", freq = "years", draw=FALSE)
```

    ##             T1_148     T2_168    S1_165
    ## 2001-01-01 0.49145 0.03411611 0.1178586

As the general summary function,
[`droughtStress()`](https://emf-creaf.github.io/medfate/reference/droughtStress.md)
allows calculating stress indices at several temporal scales. For
example the water stress index (integral of water potential values) can
be calculated and drawn for every month:

``` r

droughtStress(S, index = "WSI", freq = "months", draw=TRUE)
```

![](BasicWaterBalance_files/figure-html/unnamed-chunk-31-1.png)

Another specific summary function is
[`waterUseEfficiency()`](https://emf-creaf.github.io/medfate/reference/waterUseEfficiency.md).
This is most useful with advanced water and energy balance modeling, but
for simple water balance it calculates the ratio between photosynthesis
and transpiration at the desired scale. In this case it is equal to the
value of the input species parameter `WUE`:

``` r

waterUseEfficiency(S, type = "Stand Ag/E", freq = "months", draw=FALSE)
```

    ##            Stand Ag/E
    ## 2001-01-01   9.752461
    ## 2001-02-01   7.882513
    ## 2001-03-01   8.637627
    ## 2001-04-01   8.420369
    ## 2001-05-01   7.936004
    ## 2001-06-01   6.274635
    ## 2001-07-01   6.767208
    ## 2001-08-01   6.086707
    ## 2001-09-01   7.159872
    ## 2001-10-01   7.340686
    ## 2001-11-01   8.379917
    ## 2001-12-01   8.059213

## References

- De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
  R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
  forest inventory data to predict drought stress: the role of forest
  structural changes vs. climate changes. Agricultural and Forest
  Meteorology 213: 77-90
  (<https://doi.org/10.1016/j.agrformet.2015.06.012>).
