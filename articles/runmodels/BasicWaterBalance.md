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

    ##   widths clay sand om nitrogen  bd rfc
    ## 1    300   25   25 NA       NA 1.5  25
    ## 2    700   25   25 NA       NA 1.5  45
    ## 3   1000   25   25 NA       NA 1.5  75
    ## 4   2000   25   25 NA       NA 1.5  95

As explained in the package overview, models included in `medfate` were
primarily designed to be ran on **forest inventory plots**. Here we use
the example object provided with the package:

``` r
data(exampleforest)
exampleforest
```

    ## $treeData
    ##            Species   N   DBH Height Z50  Z95
    ## 1 Pinus halepensis 168 37.55    800 100  600
    ## 2     Quercus ilex 384 14.60    660 300 1000
    ## 
    ## $shrubData
    ##             Species Cover Height Z50  Z95
    ## 1 Quercus coccifera  3.75     80 200 1000
    ## 
    ## $herbCover
    ## [1] 10
    ## 
    ## $herbHeight
    ## [1] 20
    ## 
    ## $seedBank
    ## [1] Species Percent
    ## <0 rows> (or 0-length row.names)
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

    ## Package 'meteoland' [ver. 2.2.4]

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
    ##  28.7201692   3.3228840 386.0000000          NA   3.9023342   1.5000000 
    ## 
    ## $WaterBalance
    ##                     PET                    Rain                    Snow 
    ##              3.90233421              0.00000000              0.00000000 
    ##                 NetRain                Snowmelt                   Runon 
    ##              0.00000000              0.00000000              0.00000000 
    ##            Infiltration      InfiltrationExcess        SaturationExcess 
    ##              0.00000000              0.00000000              0.00000000 
    ##                  Runoff            DeepDrainage         CapillarityRise 
    ##              0.00000000              0.00000000              0.00000000 
    ##         SoilEvaporation       HerbTranspiration         PlantExtraction 
    ##              0.50000000              0.04872542              0.89968355 
    ##           Transpiration HydraulicRedistribution 
    ##              0.89968355              0.00000000 
    ## 
    ## $Soil
    ##           Psi HerbTranspiration HydraulicInput HydraulicOutput PlantExtraction
    ## 1 -0.03586605      0.0444001775              0     0.654144225     0.654144225
    ## 2 -0.03329772      0.0034620610              0     0.217680714     0.217680714
    ## 3 -0.03304737      0.0006078123              0     0.022325608     0.022325608
    ## 4 -0.03302988      0.0002553696              0     0.005533003     0.005533003
    ## 
    ## $Stand
    ##         LAI     LAIherb     LAIlive LAIexpanded     LAIdead          Cm 
    ##   1.7585845   0.1736369   1.5849476   1.5849476   0.0000000   1.3904846 
    ##  LgroundPAR  LgroundSWR 
    ##  40.0075402  50.7329667 
    ## 
    ## $Plants
    ##               LAI    LAIlive     FPAR AbsorbedSWRFraction Extraction
    ## T1_148 0.84874773 0.84874773 92.18285           35.076344 0.55258221
    ## T2_168 0.70557382 0.70557382 72.36365           30.444383 0.32034759
    ## S1_165 0.03062604 0.03062604 44.32407            2.366131 0.02675376
    ##        Transpiration GrossPhotosynthesis PlantPsi         DDS   StemRWC
    ## T1_148    0.55258221           3.7143297   -0.033 0.004613739 0.9998350
    ## T2_168    0.32034759           2.9300190   -0.033 0.006282261 0.9997500
    ## S1_165    0.02675376           0.1635903   -0.033 0.003088161 0.9983684
    ##          LeafRWC      LFMC      StemPLC      LeafPLC WaterBalance
    ## T1_148 0.9979679 125.86544 0.000000e+00 0.000000e+00            0
    ## T2_168 0.9986797  93.05186 3.681226e-09 2.905929e-03            0
    ## S1_165 0.9987207  96.38420 4.394058e-08 3.473351e-05            0
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

    ## [1] 0.9824193 0.9981043 0.9996972 0.9998090

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

    ## Initial plant water content (mm): 4.73001
    ## Initial soil water content (mm): 290.875
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2001]:............
    ## 
    ## Final plant water content (mm): 4.72839
    ## Final soil water content (mm): 274.93
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.00162134
    ## Plant water balance result (mm): -0.00163359
    ## Change in soil water content (mm): -15.9454
    ## Soil water balance result (mm): -15.9454
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): 0
    ## Water balance components:
    ##   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
    ##   Interception (mm) 92 Net rainfall (mm) 370
    ##   Infiltration (mm) 401 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 25  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 247
    ##   Plant extraction from soil (mm) 247  Plant water balance (mm) -0 Hydraulic redistribution (mm) 3
    ##   Runoff (mm) 21 Deep drainage (mm) 131

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
    ## 2001-01-01 0.8828475      4.869109  4.869109    0  3.4241795        0
    ## 2001-01-02 1.6375337      2.498292  2.498292    0  1.0728002        0
    ## 2001-01-03 1.3017026      0.000000  0.000000    0  0.0000000        0
    ## 2001-01-04 0.5690790      5.796973  5.796973    0  4.3636213        0
    ## 2001-01-05 1.6760567      1.884401  1.884401    0  0.7547843        0
    ## 2001-01-06 1.2077028     13.359801 13.359801    0 11.7252817        0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2001-01-01    3.4241795                  0                0      0    2.7617207
    ## 2001-01-02    1.0728002                  0                0      0    0.1899867
    ## 2001-01-03    0.0000000                  0                0      0    0.0000000
    ## 2001-01-04    4.3636213                  0                0      0    3.2495151
    ## 2001-01-05    0.7547843                  0                0      0    0.1001772
    ## 2001-01-06   11.7252817                  0                0      0    4.1214138
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2001-01-01               0          2.1073881     1.444929       0.4478948
    ## 2001-01-02               0          2.3231089     1.425491       0.5000000
    ## 2001-01-03               0          0.8014863     0.000000       0.4854151
    ## 2001-01-04               0          1.7311676     1.433352       0.1596697
    ## 2001-01-05               0          2.0365726     1.129617       0.5000000
    ## 2001-01-06               0          2.2354533     1.634519       0.3077228
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01       0.011023432       0.2035406     0.2035406
    ## 2001-01-02       0.020463709       0.3771538     0.3771538
    ## 2001-01-03       0.016266929       0.2998042     0.2998042
    ## 2001-01-04       0.007111348       0.1310349     0.1310349
    ## 2001-01-05       0.020945746       0.3860097     0.3860097
    ## 2001-01-06       0.015092529       0.2781190     0.2781190
    ##            HydraulicRedistribution
    ## 2001-01-01                       0
    ## 2001-01-02                       0
    ## 2001-01-03                       0
    ## 2001-01-04                       0
    ## 2001-01-05                       0
    ## 2001-01-06                       0

Element `Plants` is in turn a list with several dataframes with plant
output variables, for example plant water potentials are in:

``` r
head(S$Plants$PlantPsi)
```

    ##                 T1_148      T2_168      S1_165
    ## 2001-01-01 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-02 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-03 -0.03302923 -0.03301702 -0.03302310
    ## 2001-01-04 -0.03443797 -0.03388137 -0.03415838
    ## 2001-01-05 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-06 -0.03350265 -0.03329260 -0.03339730

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
stand-level variables:

``` r
df <- extract(S, "forest")
head(df)
```

    ##         date       PET Precipitation      Rain Snow    NetRain Snowmelt
    ## 1 2001-01-01 0.8828475      4.869109  4.869109    0  3.4241795        0
    ## 2 2001-01-02 1.6375337      2.498292  2.498292    0  1.0728002        0
    ## 3 2001-01-03 1.3017026      0.000000  0.000000    0  0.0000000        0
    ## 4 2001-01-04 0.5690790      5.796973  5.796973    0  4.3636213        0
    ## 5 2001-01-05 1.6760567      1.884401  1.884401    0  0.7547843        0
    ## 6 2001-01-06 1.2077028     13.359801 13.359801    0 11.7252817        0
    ##   Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 1    3.4241795                  0                0      0    2.7617207
    ## 2    1.0728002                  0                0      0    0.1899867
    ## 3    0.0000000                  0                0      0    0.0000000
    ## 4    4.3636213                  0                0      0    3.2495151
    ## 5    0.7547843                  0                0      0    0.1001772
    ## 6   11.7252817                  0                0      0    4.1214138
    ##   CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 1               0          2.1073881     1.444929       0.4478948
    ## 2               0          2.3231089     1.425491       0.5000000
    ## 3               0          0.8014863     0.000000       0.4854151
    ## 4               0          1.7311676     1.433352       0.1596697
    ## 5               0          2.0365726     1.129617       0.5000000
    ## 6               0          2.2354533     1.634519       0.3077228
    ##   HerbTranspiration PlantExtraction Transpiration HydraulicRedistribution
    ## 1       0.011023432       0.2035406     0.2035406                       0
    ## 2       0.020463709       0.3771538     0.3771538                       0
    ## 3       0.016266929       0.2998042     0.2998042                       0
    ## 4       0.007111348       0.1310349     0.1310349                       0
    ## 5       0.020945746       0.3860097     0.3860097                       0
    ## 6       0.015092529       0.2781190     0.2781190                       0
    ##        LAI   LAIherb  LAIlive LAIexpanded LAIdead       Cm LgroundPAR
    ## 1 1.758585 0.1736369 1.584948    1.584948       0 1.390485   40.00754
    ## 2 1.756533 0.1736369 1.584948    1.582896       0 1.389459   40.05271
    ## 3 1.756533 0.1736369 1.584948    1.582896       0 1.389459   40.05271
    ## 4 1.756532 0.1736369 1.584948    1.582895       0 1.389458   40.05274
    ## 5 1.756459 0.1736369 1.584948    1.582823       0 1.389422   40.05433
    ## 6 1.756459 0.1736369 1.584948    1.582823       0 1.389422   40.05433
    ##   LgroundSWR SWE
    ## 1   50.73297   0
    ## 2   50.77539   0
    ## 3   50.77539   0
    ## 4   50.77541   0
    ## 5   50.77691   0
    ## 6   50.77691   0

And a similar code can be used to daily series of cohort-level
variables:

``` r
df <- extract(S, "cohort")
head(df)
```

    ##         date cohort          species       LAI   LAIlive     FPAR
    ## 1 2001-01-01 T1_148 Pinus halepensis 0.8487477 0.8487477 92.18285
    ## 2 2001-01-02 T1_148 Pinus halepensis 0.8487477 0.8487477 92.18285
    ## 3 2001-01-03 T1_148 Pinus halepensis 0.8487477 0.8487477 92.18285
    ## 4 2001-01-04 T1_148 Pinus halepensis 0.8487477 0.8487477 92.18285
    ## 5 2001-01-05 T1_148 Pinus halepensis 0.8487477 0.8487477 92.18285
    ## 6 2001-01-06 T1_148 Pinus halepensis 0.8487477 0.8487477 92.18285
    ##   AbsorbedSWRFraction Transpiration GrossPhotosynthesis    PlantPsi LeafPLC
    ## 1            35.07634    0.12501386            1.213416 -0.03300000       0
    ## 2            35.08269    0.23181486            2.103487 -0.03300000       0
    ## 3            35.08269    0.18427233            1.773967 -0.03302923       0
    ## 4            35.08269    0.08053601            0.813088 -0.03443797       0
    ## 5            35.08269    0.23725811            2.183562 -0.03300000       0
    ## 6            35.08269    0.17094100            1.790601 -0.03350265       0
    ##   StemPLC PlantWaterBalance   LeafRWC   StemRWC     LFMC PlantStress
    ## 1       0      0.000000e+00 0.9979679 0.9998350 125.8654 0.004613739
    ## 2       0      0.000000e+00 0.9979679 0.9998350 125.8654 0.004613739
    ## 3       0     -1.319457e-06 0.9979661 0.9998348 125.8653 0.004619875
    ## 4       0     -6.358830e-05 0.9978794 0.9998278 125.8582 0.004918765
    ## 5       0      6.490775e-05 0.9979679 0.9998350 125.8654 0.004613739
    ## 6       0     -2.268898e-05 0.9979370 0.9998325 125.8629 0.004719625

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
    ## 2001-01-01 1.011397    2.41127383 1.87415609 0.5371177 1.34650405 0.42235503
    ## 2001-02-01 2.278646    0.17855109 0.08778069 0.0907704 0.03516937 0.19831578
    ## 2001-03-01 2.368035    2.41917349 2.41917349 0.0000000 1.93993344 0.01762496
    ## 2001-04-01 3.086567    0.63056064 0.29195973 0.3386009 0.13498040 0.33860091
    ## 2001-05-01 3.662604    0.76337345 0.76337345 0.0000000 0.57914176 0.00000000
    ## 2001-06-01 5.265359    0.21959509 0.21959509 0.0000000 0.15793410 0.00000000
    ## 2001-07-01 4.443053    3.27810591 3.27810591 0.0000000 2.82143775 0.00000000
    ## 2001-08-01 4.463242    1.92222891 1.92222891 0.0000000 1.55485742 0.00000000
    ## 2001-09-01 3.453891    1.30651303 1.30651303 0.0000000 1.05153159 0.00000000
    ## 2001-10-01 2.405506    1.33598175 1.33598175 0.0000000 1.05623072 0.00000000
    ## 2001-11-01 1.716591    2.20566281 1.47764599 0.7280168 1.33465670 0.72801682
    ## 2001-12-01 1.608082    0.05046181 0.05046181 0.0000000 0.02043637 0.00000000
    ##            Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2001-01-01   1.76885908         0.00000000                0 0.00000000
    ## 2001-02-01   0.23348515         0.00000000                0 0.00000000
    ## 2001-03-01   1.95755840         0.00000000                0 0.00000000
    ## 2001-04-01   0.47358131         0.00000000                0 0.00000000
    ## 2001-05-01   0.57914176         0.00000000                0 0.00000000
    ## 2001-06-01   0.15793410         0.00000000                0 0.00000000
    ## 2001-07-01   2.54076456         0.28067318                0 0.28067318
    ## 2001-08-01   1.51515475         0.03970267                0 0.03970267
    ## 2001-09-01   1.05153159         0.00000000                0 0.00000000
    ## 2001-10-01   0.94759983         0.10863089                0 0.10863089
    ## 2001-11-01   1.81629174         0.24638178                0 0.24638178
    ## 2001-12-01   0.02043637         0.00000000                0 0.00000000
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2001-01-01  1.425282878               0          0.9464146   0.52765204
    ## 2001-02-01  0.004853845               0          0.6633175   0.05261132
    ## 2001-03-01  1.132339712               0          1.1732922   0.47924005
    ## 2001-04-01  0.000000000               0          0.9159723   0.15697933
    ## 2001-05-01  0.000000000               0          1.1345364   0.18423170
    ## 2001-06-01  0.000000000               0          1.2880297   0.06166099
    ## 2001-07-01  0.000000000               0          1.6092562   0.45666816
    ## 2001-08-01  0.002505740               0          1.4982164   0.36737150
    ## 2001-09-01  0.073121104               0          1.1456698   0.25498144
    ## 2001-10-01  0.402239462               0          0.9356021   0.27975103
    ## 2001-11-01  1.239729839               0          0.6193004   0.14298929
    ## 2001-12-01  0.000000000               0          0.4357714   0.03002544
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01     0.173216782        0.01263891       0.2329068     0.2329068
    ## 2001-02-01     0.058593311        0.02847436       0.5236385     0.5236385
    ## 2001-03-01     0.120264376        0.02960281       0.5441850     0.5441850
    ## 2001-04-01     0.014325253        0.03857039       0.7060974     0.7060974
    ## 2001-05-01     0.066876539        0.04578116       0.8376470     0.8376470
    ## 2001-06-01     0.005860066        0.06530098       1.1552077     1.1552077
    ## 2001-07-01     0.096013097        0.05580089       1.0007740     1.0007740
    ## 2001-08-01     0.060662936        0.05620329       1.0139786     1.0139786
    ## 2001-09-01     0.061783297        0.04349952       0.7854056     0.7854056
    ## 2001-10-01     0.077791871        0.03030255       0.5477566     0.5477566
    ## 2001-11-01     0.063988112        0.02162232       0.3907007     0.3907007
    ## 2001-12-01     0.020160923        0.02024972       0.3653354     0.3653354
    ##            HydraulicRedistribution
    ## 2001-01-01            0.0005510393
    ## 2001-02-01            0.0002024280
    ## 2001-03-01            0.0006028065
    ## 2001-04-01            0.0000000000
    ## 2001-05-01            0.0105120588
    ## 2001-06-01            0.0000000000
    ## 2001-07-01            0.0498345905
    ## 2001-08-01            0.0285237013
    ## 2001-09-01            0.0012627675
    ## 2001-10-01            0.0042802377
    ## 2001-11-01            0.0005463483
    ## 2001-12-01            0.0000000000

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of plant cohorts by months:

``` r
summary(S, freq="months",FUN=mean, output="PlantStress")
```

    ##                 T1_148      T2_168      S1_165
    ## 2001-01-01 0.004773662 0.006381252 0.003160351
    ## 2001-02-01 0.006785896 0.007678385 0.004089571
    ## 2001-03-01 0.006289085 0.007334555 0.003851156
    ## 2001-04-01 0.011384560 0.010487139 0.006125331
    ## 2001-05-01 0.010089242 0.009590412 0.005515647
    ## 2001-06-01 0.055799936 0.029519026 0.022439959
    ## 2001-07-01 0.018312570 0.014056669 0.008788778
    ## 2001-08-01 0.008964251 0.009146389 0.005109164
    ## 2001-09-01 0.007960909 0.008426610 0.004620249
    ## 2001-10-01 0.006222902 0.007389951 0.003854981
    ## 2001-11-01 0.006657749 0.007573598 0.004018635
    ## 2001-12-01 0.009156016 0.009078819 0.005132602

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
    ## 2001-01-01      0.004613739       0.003088161  0.006282261
    ## 2001-01-02      0.004613739       0.003088161  0.006282261
    ## 2001-01-03      0.004619875       0.003090967  0.006285972
    ## 2001-01-04      0.004918765       0.003229591  0.006474811
    ## 2001-01-05      0.004613739       0.003088161  0.006282261
    ## 2001-01-06      0.004719625       0.003136505  0.006346104

Or we can combine the aggregation by species with a temporal aggregation
(here monthly averages):

``` r
summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01      0.004773662       0.003160351  0.006381252
    ## 2001-02-01      0.006785896       0.004089571  0.007678385
    ## 2001-03-01      0.006289085       0.003851156  0.007334555
    ## 2001-04-01      0.011384560       0.006125331  0.010487139
    ## 2001-05-01      0.010089242       0.005515647  0.009590412
    ## 2001-06-01      0.055799936       0.022439959  0.029519026
    ## 2001-07-01      0.018312570       0.008788778  0.014056669
    ## 2001-08-01      0.008964251       0.005109164  0.009146389
    ## 2001-09-01      0.007960909       0.004620249  0.008426610
    ## 2001-10-01      0.006222902       0.003854981  0.007389951
    ## 2001-11-01      0.006657749       0.004018635  0.007573598
    ## 2001-12-01      0.009156016       0.005132602  0.009078819

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

    ##               T1_148     T2_168     S1_165
    ## 2001-01-01 0.1074927 0.04785827 0.03979383

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
    ## 2001-01-01  10.152140
    ## 2001-02-01   8.237171
    ## 2001-03-01   9.035142
    ## 2001-04-01   8.836011
    ## 2001-05-01   8.274313
    ## 2001-06-01   6.390615
    ## 2001-07-01   7.024143
    ## 2001-08-01   6.304413
    ## 2001-09-01   7.476309
    ## 2001-10-01   7.675817
    ## 2001-11-01   8.756408
    ## 2001-12-01   8.435627

## References

- De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
  R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
  forest inventory data to predict drought stress: the role of forest
  structural changes vs. climate changes. Agricultural and Forest
  Meteorology 213: 77-90
  (<https://doi.org/10.1016/j.agrformet.2015.06.012>).
