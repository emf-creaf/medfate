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
    ##               0.5000000               0.0000000               0.8249229 
    ##           Transpiration  MistletoeTranspiration HydraulicRedistribution 
    ##               0.8249229               0.0000000               0.0000000 
    ## 
    ## $Soil
    ##           Psi HerbTranspiration HydraulicInput HydraulicOutput PlantExtraction
    ## 1 -0.03572510                 0              0     0.641887425     0.641887425
    ## 2 -0.03322327                 0              0     0.166030731     0.166030731
    ## 3 -0.03302923                 0              0     0.014152858     0.014152858
    ## 4 -0.03301472                 0              0     0.002851839     0.002851839
    ## 
    ## $Stand
    ##          LAI      LAIherb      LAIlive  LAIexpanded      LAIdead LAImistletoe 
    ##     1.451663     0.000000     1.451663     1.451663     0.000000     0.000000 
    ##           Cm   LgroundPAR   LgroundSWR 
    ##     1.102945    46.733767    56.922251 
    ## 
    ## $Plants
    ##               LAI    LAIlive     FPAR AbsorbedSWRFraction Extraction
    ## T1_148 0.75422783 0.75422783 92.71417           32.805782  0.4933618
    ## T2_168 0.64411380 0.64411380 76.59611           29.007524  0.2913785
    ## S1_165 0.05332129 0.05332129 48.02393            3.738549  0.0401826
    ##        Transpiration MistletoeTranspiration GrossPhotosynthesis PlantPsi
    ## T1_148     0.4933618                      0           3.3262642   -0.033
    ## T2_168     0.2913785                      0           2.6865369   -0.033
    ## S1_165     0.0401826                      0           0.2529852   -0.033
    ##                DDS   StemRWC   LeafRWC      LFMC StemPLC      LeafPLC
    ## T1_148 0.004613739 0.9998346 0.9981628 122.02311       0 0.0000000000
    ## T2_168 0.006282261 0.9997405 0.9991208 109.14738       0 0.0009945146
    ## S1_165 0.007310922 0.9995007 0.9987189  98.29725       0 0.0000000000
    ##         WaterBalance
    ## T1_148 -4.440892e-16
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

    ## [1] 0.9832504 0.9985767 0.9998132 0.9999059

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

    ## Initial plant water content (mm): 4.20666
    ## Initial soil water content (mm): 290.875
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2001]:............
    ## 
    ## Final plant water content (mm): 4.20517
    ## Final soil water content (mm): 276.623
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.0014939
    ## Plant water balance result (mm): -0.0014939
    ## Change in soil water content (mm): -14.2514
    ## Soil water balance result (mm): -14.2514
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): 0
    ## Water balance components:
    ##   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
    ##   Interception (mm) 76 Net rainfall (mm) 386
    ##   Infiltration (mm) 415 Infiltration excess (mm) 22 Saturation excess (mm) 0 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 27  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 226  Mistletoe transpiration (mm) 0
    ##   Plant extraction from soil (mm) 226  Plant water balance (mm) -0 Hydraulic redistribution (mm) 1
    ##   Runoff (mm) 22 Deep drainage (mm) 176

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
    ## 2001-01-01 0.8828475      4.869109  4.869109    0  3.7161954        0
    ## 2001-01-02 1.6375337      2.498292  2.498292    0  1.3602886        0
    ## 2001-01-03 1.3017026      0.000000  0.000000    0  0.0000000        0
    ## 2001-01-04 0.5690790      5.796973  5.796973    0  4.6539798        0
    ## 2001-01-05 1.6760567      1.884401  1.884401    0  0.8806518        0
    ## 2001-01-06 1.2077028     13.359801 13.359801    0 12.0368347        0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2001-01-01    3.7161954                  0                0      0    3.0295684
    ## 2001-01-02    1.3602886                  0                0      0    0.5141268
    ## 2001-01-03    0.0000000                  0                0      0    0.0000000
    ## 2001-01-04    4.6539798                  0                0      0    3.6009039
    ## 2001-01-05    0.8806518                  0                0      0    0.1632798
    ## 2001-01-06   12.0368347                  0                0      0    4.1214138
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2001-01-01               0          1.8395405     1.152913       0.5000000
    ## 2001-01-02               0          1.9841648     1.138003       0.5000000
    ## 2001-01-03               0          0.7751697     0.000000       0.5000000
    ## 2001-01-04               0          1.4208993     1.142993       0.1576421
    ## 2001-01-05               0          1.8580549     1.003750       0.5000000
    ## 2001-01-06               0          1.9597298     1.322966       0.3814785
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01                 0       0.1866270     0.1866270
    ## 2001-01-02                 0       0.3461618     0.3461618
    ## 2001-01-03                 0       0.2751697     0.2751697
    ## 2001-01-04                 0       0.1202641     0.1202641
    ## 2001-01-05                 0       0.3543052     0.3543052
    ## 2001-01-06                 0       0.2552853     0.2552853
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
    ## 2001-01-03 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-04 -0.03459905 -0.03387617 -0.03416471
    ## 2001-01-05 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-06 -0.03329973 -0.03315810 -0.03321465

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
    ##  1 2001-01-01   0.883          4.87    4.87   0      3.72      0           3.72 
    ##  2 2001-01-02   1.64           2.50    2.50   0      1.36      0           1.36 
    ##  3 2001-01-03   1.30           0       0      0      0         0           0    
    ##  4 2001-01-04   0.569          5.80    5.80   0      4.65      0           4.65 
    ##  5 2001-01-05   1.68           1.88    1.88   0      0.881     0           0.881
    ##  6 2001-01-06   1.21          13.4    13.4    0     12.0       0          12.0  
    ##  7 2001-01-07   0.637          5.38    0      5.38   0         0           0    
    ##  8 2001-01-08   0.832          0       0      0      0         0           0    
    ##  9 2001-01-09   1.98           0       0      0      0         0           0    
    ## 10 2001-01-10   0.829          5.12    5.12   0      3.96      5.38        9.34 
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
    ##  1 2001-01-01 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ##  2 2001-01-02 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ##  3 2001-01-03 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ##  4 2001-01-04 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ##  5 2001-01-05 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ##  6 2001-01-06 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ##  7 2001-01-07 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ##  8 2001-01-08 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ##  9 2001-01-09 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
    ## 10 2001-01-10 T1_148 Pinus halepensis    0.754   0.754 92.7                32.8
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
    ## 2001-01-01 1.011397    2.41127383 1.87415609 0.5371177 1.43878978 0.42235503
    ## 2001-02-01 2.278646    0.17855109 0.08778069 0.0907704 0.04102322 0.19831578
    ## 2001-03-01 2.368035    2.41917349 2.41917349 0.0000000 2.01380403 0.01762496
    ## 2001-04-01 3.086567    0.63056064 0.29195973 0.3386009 0.15663121 0.33860091
    ## 2001-05-01 3.662604    0.76337345 0.76337345 0.0000000 0.60951681 0.00000000
    ## 2001-06-01 5.265359    0.21959509 0.21959509 0.0000000 0.16851902 0.00000000
    ## 2001-07-01 4.443053    3.27810591 3.27810591 0.0000000 2.90709159 0.00000000
    ## 2001-08-01 4.463242    1.92222891 1.92222891 0.0000000 1.61854736 0.00000000
    ## 2001-09-01 3.453891    1.30651303 1.30651303 0.0000000 1.09567120 0.00000000
    ## 2001-10-01 2.405506    1.33598175 1.33598175 0.0000000 1.10512424 0.00000000
    ## 2001-11-01 1.716591    2.20566281 1.47764599 0.7280168 1.36280045 0.72801682
    ## 2001-12-01 1.608082    0.05046181 0.05046181 0.0000000 0.02358271 0.00000000
    ##            Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2001-01-01   1.86114482         0.00000000                0 0.00000000
    ## 2001-02-01   0.23933900         0.00000000                0 0.00000000
    ## 2001-03-01   2.03142899         0.00000000                0 0.00000000
    ## 2001-04-01   0.49523212         0.00000000                0 0.00000000
    ## 2001-05-01   0.60951681         0.00000000                0 0.00000000
    ## 2001-06-01   0.16851902         0.00000000                0 0.00000000
    ## 2001-07-01   2.60732111         0.29977048                0 0.29977048
    ## 2001-08-01   1.57225299         0.04629436                0 0.04629436
    ## 2001-09-01   1.09567120         0.00000000                0 0.00000000
    ## 2001-10-01   0.99215768         0.11296655                0 0.11296655
    ## 2001-11-01   1.83659147         0.25422580                0 0.25422580
    ## 2001-12-01   0.02358271         0.00000000                0 0.00000000
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2001-01-01   1.52395257               0          0.8371713   0.43536630
    ## 2001-02-01   0.01786488               0          0.5933041   0.04675746
    ## 2001-03-01   1.27401087               0          1.0324278   0.40536945
    ## 2001-04-01   0.00000000               0          0.8025344   0.13532852
    ## 2001-05-01   0.11978072               0          1.0015633   0.15385664
    ## 2001-06-01   0.00000000               0          1.0872480   0.05107606
    ## 2001-07-01   0.21624027               0          1.3955023   0.37101432
    ## 2001-08-01   0.47541500               0          1.3192700   0.30368156
    ## 2001-09-01   0.18109630               0          1.0060257   0.21084183
    ## 2001-10-01   0.57324383               0          0.8219592   0.23085752
    ## 2001-11-01   1.34218565               0          0.5304826   0.11484554
    ## 2001-12-01   0.00000000               0          0.3856804   0.02687911
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01      0.18803212                 0       0.2137728     0.2137728
    ## 2001-02-01      0.06576432                 0       0.4807823     0.4807823
    ## 2001-03-01      0.12722866                 0       0.4998297     0.4998297
    ## 2001-04-01      0.01842982                 0       0.6487760     0.6487760
    ## 2001-05-01      0.07798656                 0       0.7697201     0.7697201
    ## 2001-06-01      0.00642879                 0       1.0297432     1.0297432
    ## 2001-07-01      0.10479231                 0       0.9196956     0.9196956
    ## 2001-08-01      0.07516169                 0       0.9404267     0.9404267
    ## 2001-09-01      0.06715962                 0       0.7280242     0.7280242
    ## 2001-10-01      0.08325073                 0       0.5078510     0.5078510
    ## 2001-11-01      0.05348927                 0       0.3621478     0.3621478
    ## 2001-12-01      0.02027398                 0       0.3385274     0.3385274
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2001-01-01                      0            8.632154e-04
    ## 2001-02-01                      0            0.000000e+00
    ## 2001-03-01                      0            1.159206e-03
    ## 2001-04-01                      0            1.331529e-03
    ## 2001-05-01                      0            9.685021e-05
    ## 2001-06-01                      0            0.000000e+00
    ## 2001-07-01                      0            1.375693e-02
    ## 2001-08-01                      0            2.406171e-03
    ## 2001-09-01                      0            1.132185e-03
    ## 2001-10-01                      0            7.540357e-04
    ## 2001-11-01                      0            8.536563e-04
    ## 2001-12-01                      0            0.000000e+00

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of plant cohorts by months:

``` r

summary(S, freq="months",FUN=mean, output="PlantStress")
```

    ##                 T1_148      T2_168      S1_165
    ## 2001-01-01 0.004766131 0.006364386 0.007453863
    ## 2001-02-01 0.006758337 0.007483040 0.009409180
    ## 2001-03-01 0.006292729 0.007182275 0.008911940
    ## 2001-04-01 0.011513809 0.009192987 0.013080402
    ## 2001-05-01 0.010537544 0.008721668 0.012169904
    ## 2001-06-01 0.110666204 0.017511721 0.040593515
    ## 2001-07-01 0.031461048 0.009544880 0.015551968
    ## 2001-08-01 0.008449392 0.008107148 0.010740038
    ## 2001-09-01 0.008036716 0.008017727 0.010480940
    ## 2001-10-01 0.006066233 0.007122774 0.008755170
    ## 2001-11-01 0.006744639 0.007403308 0.009319030
    ## 2001-12-01 0.009727859 0.008795514 0.011981015

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
    ## 2001-01-03      0.004613739       0.007310922  0.006282261
    ## 2001-01-04      0.004953334       0.007646834  0.006473673
    ## 2001-01-05      0.004613739       0.007310922  0.006282261
    ## 2001-01-06      0.004676784       0.007372574  0.006316747

Or we can combine the aggregation by species with a temporal aggregation
(here monthly averages):

``` r

summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01      0.004766131       0.007453863  0.006364386
    ## 2001-02-01      0.006758337       0.009409180  0.007483040
    ## 2001-03-01      0.006292729       0.008911940  0.007182275
    ## 2001-04-01      0.011513809       0.013080402  0.009192987
    ## 2001-05-01      0.010537544       0.012169904  0.008721668
    ## 2001-06-01      0.110666204       0.040593515  0.017511721
    ## 2001-07-01      0.031461048       0.015551968  0.009544880
    ## 2001-08-01      0.008449392       0.010740038  0.008107148
    ## 2001-09-01      0.008036716       0.010480940  0.008017727
    ## 2001-10-01      0.006066233       0.008755170  0.007122774
    ## 2001-11-01      0.006744639       0.009319030  0.007403308
    ## 2001-12-01      0.009727859       0.011981015  0.008795514

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
    ## 2001-01-01 0.2549964 0.02343213 0.06518573

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
    ## 2001-01-01  10.182667
    ## 2001-02-01   8.265741
    ## 2001-03-01   9.067584
    ## 2001-04-01   8.872904
    ## 2001-05-01   8.304651
    ## 2001-06-01   6.443255
    ## 2001-07-01   7.065821
    ## 2001-08-01   6.325558
    ## 2001-09-01   7.510940
    ## 2001-10-01   7.712734
    ## 2001-11-01   8.797458
    ## 2001-12-01   8.478022

## References

- De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
  R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
  forest inventory data to predict drought stress: the role of forest
  structural changes vs. climate changes. Agricultural and Forest
  Meteorology 213: 77-90
  (<https://doi.org/10.1016/j.agrformet.2015.06.012>).
