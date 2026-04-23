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
    ##               3.9023342               0.0000000               0.0000000 
    ##                 NetRain                Snowmelt                   Runon 
    ##               0.0000000               0.0000000               0.0000000 
    ##            Infiltration      InfiltrationExcess        SaturationExcess 
    ##               0.0000000               0.0000000               0.0000000 
    ##                  Runoff            DeepDrainage         CapillarityRise 
    ##               0.0000000               0.0000000               0.0000000 
    ##         SoilEvaporation       HerbTranspiration         PlantExtraction 
    ##               0.5000000               0.0000000               0.8996836 
    ##           Transpiration HydraulicRedistribution 
    ##               0.8996836               0.0000000 
    ## 
    ## $Soil
    ##           Psi HerbTranspiration HydraulicInput HydraulicOutput PlantExtraction
    ## 1 -0.03587845                 0              0     0.703516265     0.703516265
    ## 2 -0.03323977                 0              0     0.178256444     0.178256444
    ## 3 -0.03303089                 0              0     0.014956372     0.014956372
    ## 4 -0.03301525                 0              0     0.002954471     0.002954471
    ## 
    ## $Stand
    ##         LAI     LAIherb     LAIlive LAIexpanded     LAIdead          Cm 
    ##    1.584948    0.000000    1.584948    1.584948    0.000000    1.216848 
    ##  LgroundPAR  LgroundSWR 
    ##   43.636170   54.102798 
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
    ##          LeafRWC      LFMC StemPLC      LeafPLC WaterBalance
    ## T1_148 0.9979679 125.86544       0 0.0000000000 0.000000e+00
    ## T2_168 0.9986797  93.05186       0 0.0009945146 0.000000e+00
    ## S1_165 0.9987207  96.38420       0 0.0000000000 3.469447e-18
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

    ## [1] 0.9823464 0.9984719 0.9998026 0.9999025

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

    ## Initial plant water content (mm): 4.69853
    ## Initial soil water content (mm): 290.875
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2001]:............
    ## 
    ## Final plant water content (mm): 4.69673
    ## Final soil water content (mm): 275.597
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.00180112
    ## Plant water balance result (mm): -0.00180112
    ## Change in soil water content (mm): -15.2779
    ## Soil water balance result (mm): -15.2779
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): -7.10543e-15
    ## Water balance components:
    ##   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
    ##   Interception (mm) 83 Net rainfall (mm) 379
    ##   Infiltration (mm) 409 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 25  Herbaceous transpiration (mm) 0 Woody plant transpiration (mm) 248
    ##   Plant extraction from soil (mm) 248  Plant water balance (mm) -0 Hydraulic redistribution (mm) 2
    ##   Runoff (mm) 21 Deep drainage (mm) 152

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
    ## 2001-01-01 0.8828475      4.869109  4.869109    0  3.6000026        0
    ## 2001-01-02 1.6375337      2.498292  2.498292    0  1.2454970        0
    ## 2001-01-03 1.3017026      0.000000  0.000000    0  0.0000000        0
    ## 2001-01-04 0.5690790      5.796973  5.796973    0  4.5381440        0
    ## 2001-01-05 1.6760567      1.884401  1.884401    0  0.8222806        0
    ## 2001-01-06 1.2077028     13.359801 13.359801    0 11.9109817        0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2001-01-01    3.6000026                  0                0      0    2.9188167
    ## 2001-01-02    1.2454970                  0                0      0    0.3679634
    ## 2001-01-03    0.0000000                  0                0      0    0.0000000
    ## 2001-01-04    4.5381440                  0                0      0    3.4528705
    ## 2001-01-05    0.8222806                  0                0      0    0.1416065
    ## 2001-01-06   11.9109817                  0                0      0    4.1214138
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2001-01-01               0          1.9502921     1.269106       0.4776452
    ## 2001-01-02               0          2.1303282     1.252795       0.5000000
    ## 2001-01-03               0          0.8001077     0.000000       0.5000000
    ## 2001-01-04               0          1.5439948     1.258829       0.1540028
    ## 2001-01-05               0          1.9485359     1.062121       0.5000000
    ## 2001-01-06               0          2.0621670     1.448819       0.3349340
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01                 0       0.2035406     0.2035406
    ## 2001-01-02                 0       0.3775336     0.3775336
    ## 2001-01-03                 0       0.3001077     0.3001077
    ## 2001-01-04                 0       0.1311630     0.1311630
    ## 2001-01-05                 0       0.3864150     0.3864150
    ## 2001-01-06                 0       0.2784142     0.2784142
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
    ## 2001-01-03 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-04 -0.03464644 -0.03390356 -0.03420012
    ## 2001-01-05 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-06 -0.03345144 -0.03323815 -0.03332335

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

    ## # A tibble: 365 × 29
    ##    date           PET Precipitation    Rain   Snow NetRain Snowmelt Infiltration
    ##    <date>     [L/m^2]       [L/m^2] [L/m^2] [L/m^… [L/m^2]  [L/m^2]      [L/m^2]
    ##  1 2001-01-01   0.883          4.87    4.87   0      3.60      0           3.60 
    ##  2 2001-01-02   1.64           2.50    2.50   0      1.25      0           1.25 
    ##  3 2001-01-03   1.30           0       0      0      0         0           0    
    ##  4 2001-01-04   0.569          5.80    5.80   0      4.54      0           4.54 
    ##  5 2001-01-05   1.68           1.88    1.88   0      0.822     0           0.822
    ##  6 2001-01-06   1.21          13.4    13.4    0     11.9       0          11.9  
    ##  7 2001-01-07   0.637          5.38    0      5.38   0         0           0    
    ##  8 2001-01-08   0.832          0       0      0      0         0           0    
    ##  9 2001-01-09   1.98           0       0      0      0         0           0    
    ## 10 2001-01-10   0.829          5.12    5.12   0      3.85      5.38        9.23 
    ## # ℹ 355 more rows
    ## # ℹ 21 more variables: InfiltrationExcess [L/m^2], SaturationExcess [L/m^2],
    ## #   Runoff [L/m^2], DeepDrainage [L/m^2], CapillarityRise [L/m^2],
    ## #   Evapotranspiration [L/m^2], Interception [L/m^2], SoilEvaporation [L/m^2],
    ## #   HerbTranspiration [L/m^2], PlantExtraction [L/m^2], Transpiration [L/m^2],
    ## #   HydraulicRedistribution [L/m^2], LAI [m^2/m^2], LAIherb [m^2/m^2],
    ## #   LAIlive [m^2/m^2], LAIexpanded [m^2/m^2], LAIdead [m^2/m^2], Cm [L/m^2], …

And a similar code can be used to daily series of cohort-level
variables:

``` r
extract(S, "cohort", addunits = TRUE) |>
  tibble::as_tibble()
```

    ## # A tibble: 1,095 × 17
    ##    date       cohort species               LAI LAIlive FPAR AbsorbedSWRFraction
    ##    <date>     <chr>  <chr>            [m^2/m^… [m^2/m…  [%]               <dbl>
    ##  1 2001-01-01 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ##  2 2001-01-02 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ##  3 2001-01-03 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ##  4 2001-01-04 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ##  5 2001-01-05 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ##  6 2001-01-06 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ##  7 2001-01-07 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ##  8 2001-01-08 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ##  9 2001-01-09 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
    ## 10 2001-01-10 T1_148 Pinus halepensis    0.849   0.849 92.2                35.1
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
    ## 2001-01-01 1.011397    2.41127383 1.87415609 0.5371177 1.40082592 0.42235503
    ## 2001-02-01 2.278646    0.17855109 0.08778069 0.0907704 0.03830413 0.19831578
    ## 2001-03-01 2.368035    2.41917349 2.41917349 0.0000000 1.98241650 0.01762496
    ## 2001-04-01 3.086567    0.63056064 0.29195973 0.3386009 0.14714336 0.33860091
    ## 2001-05-01 3.662604    0.76337345 0.76337345 0.0000000 0.59589515 0.00000000
    ## 2001-06-01 5.265359    0.21959509 0.21959509 0.0000000 0.16400406 0.00000000
    ## 2001-07-01 4.443053    3.27810591 3.27810591 0.0000000 2.87125513 0.00000000
    ## 2001-08-01 4.463242    1.92222891 1.92222891 0.0000000 1.59138259 0.00000000
    ## 2001-09-01 3.453891    1.30651303 1.30651303 0.0000000 1.07683791 0.00000000
    ## 2001-10-01 2.405506    1.33598175 1.33598175 0.0000000 1.08372293 0.00000000
    ## 2001-11-01 1.716591    2.20566281 1.47764599 0.7280168 1.35115693 0.72801682
    ## 2001-12-01 1.608082    0.05046181 0.05046181 0.0000000 0.02201960 0.00000000
    ##            Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2001-01-01    1.8231810         0.00000000                0 0.00000000
    ## 2001-02-01    0.2366199         0.00000000                0 0.00000000
    ## 2001-03-01    2.0000415         0.00000000                0 0.00000000
    ## 2001-04-01    0.4857443         0.00000000                0 0.00000000
    ## 2001-05-01    0.5958952         0.00000000                0 0.00000000
    ## 2001-06-01    0.1640041         0.00000000                0 0.00000000
    ## 2001-07-01    2.5907436         0.28051151                0 0.28051151
    ## 2001-08-01    1.5477774         0.04360518                0 0.04360518
    ## 2001-09-01    1.0768379         0.00000000                0 0.00000000
    ## 2001-10-01    0.9725580         0.11116492                0 0.11116492
    ## 2001-11-01    1.8306931         0.24848066                0 0.24848066
    ## 2001-12-01    0.0220196         0.00000000                0 0.00000000
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2001-01-01   1.48164455               0          0.8855736   0.47333016
    ## 2001-02-01   0.01041079               0          0.6344160   0.04947656
    ## 2001-03-01   1.20180007               0          1.1032069   0.43675699
    ## 2001-04-01   0.00000000               0          0.8660697   0.14481637
    ## 2001-05-01   0.01086453               0          1.0792371   0.16747830
    ## 2001-06-01   0.00000000               0          1.2107556   0.05559103
    ## 2001-07-01   0.00000000               0          1.5101116   0.40685078
    ## 2001-08-01   0.33990559               0          1.4202530   0.33084632
    ## 2001-09-01   0.11327541               0          1.0846412   0.22967512
    ## 2001-10-01   0.49130531               0          0.8840700   0.25225882
    ## 2001-11-01   1.28544450               0          0.5785263   0.12648906
    ## 2001-12-01   0.00000000               0          0.4172253   0.02844221
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01     0.179100203                 0       0.2331432     0.2331432
    ## 2001-02-01     0.060665383                 0       0.5242741     0.5242741
    ## 2001-03-01     0.121389865                 0       0.5450600     0.5450600
    ## 2001-04-01     0.013992087                 0       0.7072612     0.7072612
    ## 2001-05-01     0.072295145                 0       0.8394637     0.8394637
    ## 2001-06-01     0.005619203                 0       1.1495454     1.1495454
    ## 2001-07-01     0.094276069                 0       1.0089848     1.0089848
    ## 2001-08-01     0.064058755                 0       1.0253479     1.0253479
    ## 2001-09-01     0.061111233                 0       0.7938548     0.7938548
    ## 2001-10-01     0.078036293                 0       0.5537749     0.5537749
    ## 2001-11-01     0.057106071                 0       0.3949312     0.3949312
    ## 2001-12-01     0.019569384                 0       0.3692137     0.3692137
    ##            HydraulicRedistribution
    ## 2001-01-01            9.590242e-04
    ## 2001-02-01            0.000000e+00
    ## 2001-03-01            1.312838e-03
    ## 2001-04-01            0.000000e+00
    ## 2001-05-01            1.313720e-03
    ## 2001-06-01            1.401205e-05
    ## 2001-07-01            4.718319e-02
    ## 2001-08-01            5.858901e-03
    ## 2001-09-01            1.613232e-03
    ## 2001-10-01            5.009414e-03
    ## 2001-11-01            9.687727e-04
    ## 2001-12-01            0.000000e+00

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of plant cohorts by months:

``` r
summary(S, freq="months",FUN=mean, output="PlantStress")
```

    ##                 T1_148      T2_168      S1_165
    ## 2001-01-01 0.004784746 0.006375382 0.003156729
    ## 2001-02-01 0.006989065 0.007644662 0.004082373
    ## 2001-03-01 0.006473653 0.007315582 0.003851819
    ## 2001-04-01 0.012180853 0.010386440 0.006145362
    ## 2001-05-01 0.010558457 0.009418937 0.005455470
    ## 2001-06-01 0.071051366 0.031630680 0.025234052
    ## 2001-07-01 0.022169351 0.014110460 0.009365977
    ## 2001-08-01 0.008900049 0.008679593 0.004851299
    ## 2001-09-01 0.008389474 0.008411343 0.004649397
    ## 2001-10-01 0.006314736 0.007283814 0.003811202
    ## 2001-11-01 0.006928037 0.007565379 0.004034940
    ## 2001-12-01 0.009788701 0.009118803 0.005198114

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
    ## 2001-01-03      0.004613739       0.003088161  0.006282261
    ## 2001-01-04      0.004963522       0.003234714  0.006479670
    ## 2001-01-05      0.004613739       0.003088161  0.006282261
    ## 2001-01-06      0.004708802       0.003127494  0.006334217

Or we can combine the aggregation by species with a temporal aggregation
(here monthly averages):

``` r
summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01      0.004784746       0.003156729  0.006375382
    ## 2001-02-01      0.006989065       0.004082373  0.007644662
    ## 2001-03-01      0.006473653       0.003851819  0.007315582
    ## 2001-04-01      0.012180853       0.006145362  0.010386440
    ## 2001-05-01      0.010558457       0.005455470  0.009418937
    ## 2001-06-01      0.071051366       0.025234052  0.031630680
    ## 2001-07-01      0.022169351       0.009365977  0.014110460
    ## 2001-08-01      0.008900049       0.004851299  0.008679593
    ## 2001-09-01      0.008389474       0.004649397  0.008411343
    ## 2001-10-01      0.006314736       0.003811202  0.007283814
    ## 2001-11-01      0.006928037       0.004034940  0.007565379
    ## 2001-12-01      0.009788701       0.005198114  0.009118803

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
    ## 2001-01-01 0.1468479 0.05409005 0.04769724

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
    ## 2001-01-01  10.153011
    ## 2001-02-01   8.238213
    ## 2001-03-01   9.036727
    ## 2001-04-01   8.838613
    ## 2001-05-01   8.276795
    ## 2001-06-01   6.399879
    ## 2001-07-01   7.034622
    ## 2001-08-01   6.309183
    ## 2001-09-01   7.485318
    ## 2001-10-01   7.685658
    ## 2001-11-01   8.767267
    ## 2001-12-01   8.447074

## References

- De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
  R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
  forest inventory data to predict drought stress: the role of forest
  structural changes vs. climate changes. Agricultural and Forest
  Meteorology 213: 77-90
  (<https://doi.org/10.1016/j.agrformet.2015.06.012>).
