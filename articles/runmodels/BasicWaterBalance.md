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
    ## T1_158 158  Pinus halepensis
    ## T2_179 179      Quercus ilex
    ## S1_176 176 Quercus coccifera
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
    ##               0.5000000               0.0000000               0.8571652 
    ##           Transpiration  MistletoeTranspiration HydraulicRedistribution 
    ##               0.8571652               0.0000000               0.0000000 
    ## 
    ## $Soil
    ##           Psi HerbTranspiration HydraulicInput HydraulicOutput PlantExtraction
    ## 1 -0.03580806                 0              0     0.675263291     0.675263291
    ## 2 -0.03322204                 0              0     0.165117101     0.165117101
    ## 3 -0.03302885                 0              0     0.013970988     0.013970988
    ## 4 -0.03301452                 0              0     0.002813863     0.002813863
    ## 
    ## $Stand
    ##          LAI      LAIherb      LAIlive  LAIexpanded      LAIdead LAImistletoe 
    ##     1.497716     0.000000     1.497716     1.497716     0.000000     0.000000 
    ##           Cm   LgroundPAR   LgroundSWR 
    ##     1.160807    45.723900    56.008543 
    ## 
    ## $Plants
    ##               LAI    LAIlive     FPAR AbsorbedSWRFraction Extraction
    ## T1_158 0.82389823 0.82389823 92.06855            34.81224 0.53176459
    ## T2_179 0.62107792 0.62107792 75.49574            27.96302 0.28530096
    ## S1_176 0.05274013 0.05274013 46.97224             3.66289 0.04009969
    ##        Transpiration MistletoeTranspiration GrossPhotosynthesis PlantPsi
    ## T1_158    0.53176459                      0           3.5720759   -0.033
    ## T2_179    0.28530096                      0           2.6251315   -0.033
    ## S1_176    0.04009969                      0           0.2504349   -0.033
    ##                DDS   StemRWC   LeafRWC      LFMC StemPLC LeafPLC WaterBalance
    ## T1_158 0.004613739 0.9998346 0.9981628 122.02311       0       0            0
    ## T2_179 0.006282261 0.9997405 0.9991208 109.14738       0       0            0
    ## S1_176 0.004541797 0.9995054 0.9987742  98.30066       0       0            0
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

    ## [1] 0.9827608 0.9985845 0.9998156 0.9999071

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

    ## Initial plant water content (mm): 4.67289
    ## Initial soil water content (mm): 290.875
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2001]:............
    ## 
    ## Final plant water content (mm): 4.67119
    ## Final soil water content (mm): 276.182
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.00169593
    ## Plant water balance result (mm): -0.00169593
    ## Change in soil water content (mm): -14.693
    ## Soil water balance result (mm): -14.693
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): -7.10543e-15
    ## Water balance components:
    ##   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
    ##   Interception (mm) 79 Net rainfall (mm) 383
    ##   Infiltration (mm) 413 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 26  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 235  Mistletoe transpiration (mm) 0
    ##   Plant extraction from soil (mm) 235  Plant water balance (mm) -0 Hydraulic redistribution (mm) 1
    ##   Runoff (mm) 21 Deep drainage (mm) 166

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
    ## 2001-01-01 0.8828475      4.869109  4.869109    0  3.6578430        0
    ## 2001-01-02 1.6375337      2.498292  2.498292    0  1.3026155        0
    ## 2001-01-03 1.3017026      0.000000  0.000000    0  0.0000000        0
    ## 2001-01-04 0.5690790      5.796973  5.796973    0  4.5956523        0
    ## 2001-01-05 1.6760567      1.884401  1.884401    0  0.8616218        0
    ## 2001-01-06 1.2077028     13.359801 13.359801    0 11.9754283        0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2001-01-01    3.6578430                  0                0      0    2.9694516
    ## 2001-01-02    1.3026155                  0                0      0    0.4429239
    ## 2001-01-03    0.0000000                  0                0      0    0.0000000
    ## 2001-01-04    4.5956523                  0                0      0    3.5291110
    ## 2001-01-05    0.8616218                  0                0      0    0.1585395
    ## 2001-01-06   11.9754283                  0                0      0    4.1214138
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2001-01-01               0          1.8996573     1.211266       0.4944700
    ## 2001-01-02               0          2.0553677     1.195676       0.5000000
    ## 2001-01-03               0          0.7859248     0.000000       0.5000000
    ## 2001-01-04               0          1.4819371     1.201321       0.1556520
    ## 2001-01-05               0          1.8909330     1.022780       0.5000000
    ## 2001-01-06               0          2.0111056     1.384372       0.3614728
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01                 0       0.1939214     0.1939214
    ## 2001-01-02                 0       0.3596916     0.3596916
    ## 2001-01-03                 0       0.2859248     0.2859248
    ## 2001-01-04                 0       0.1249645     0.1249645
    ## 2001-01-05                 0       0.3681534     0.3681534
    ## 2001-01-06                 0       0.2652605     0.2652605
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

    ##                 T1_158      T2_179      S1_176
    ## 2001-01-01 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-02 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-03 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-04 -0.03462455 -0.03388946 -0.03418295
    ## 2001-01-05 -0.03300000 -0.03300000 -0.03300000
    ## 2001-01-06 -0.03336168 -0.03319078 -0.03325903

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
    ##  1 2001-01-01   0.883          4.87    4.87   0      3.66      0           3.66 
    ##  2 2001-01-02   1.64           2.50    2.50   0      1.30      0           1.30 
    ##  3 2001-01-03   1.30           0       0      0      0         0           0    
    ##  4 2001-01-04   0.569          5.80    5.80   0      4.60      0           4.60 
    ##  5 2001-01-05   1.68           1.88    1.88   0      0.862     0           0.862
    ##  6 2001-01-06   1.21          13.4    13.4    0     12.0       0          12.0  
    ##  7 2001-01-07   0.637          5.38    0      5.38   0         0           0    
    ##  8 2001-01-08   0.832          0       0      0      0         0           0    
    ##  9 2001-01-09   1.98           0       0      0      0         0           0    
    ## 10 2001-01-10   0.829          5.12    5.12   0      3.91      5.38        9.28 
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
    ##  1 2001-01-01 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ##  2 2001-01-02 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ##  3 2001-01-03 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ##  4 2001-01-04 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ##  5 2001-01-05 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ##  6 2001-01-06 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ##  7 2001-01-07 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ##  8 2001-01-08 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ##  9 2001-01-09 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
    ## 10 2001-01-10 T1_158 Pinus halepensis    0.824   0.824 92.1                34.8
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
    ## 2001-01-01 1.011397    2.41127383 1.87415609 0.5371177 1.42114546 0.42235503
    ## 2001-02-01 2.278646    0.17855109 0.08778069 0.0907704 0.04013675 0.19831578
    ## 2001-03-01 2.368035    2.41917349 2.41917349 0.0000000 2.00031636 0.01762496
    ## 2001-04-01 3.086567    0.63056064 0.29195973 0.3386009 0.15286403 0.33860091
    ## 2001-05-01 3.662604    0.76337345 0.76337345 0.0000000 0.60312887 0.00000000
    ## 2001-06-01 5.265359    0.21959509 0.21959509 0.0000000 0.16640364 0.00000000
    ## 2001-07-01 4.443053    3.27810591 3.27810591 0.0000000 2.88946818 0.00000000
    ## 2001-08-01 4.463242    1.92222891 1.92222891 0.0000000 1.60573214 0.00000000
    ## 2001-09-01 3.453891    1.30651303 1.30651303 0.0000000 1.08679369 0.00000000
    ## 2001-10-01 2.405506    1.33598175 1.33598175 0.0000000 1.09482469 0.00000000
    ## 2001-11-01 1.716591    2.20566281 1.47764599 0.7280168 1.35693642 0.72801682
    ## 2001-12-01 1.608082    0.05046181 0.05046181 0.0000000 0.02307311 0.00000000
    ##            Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2001-01-01   1.84350049         0.00000000                0 0.00000000
    ## 2001-02-01   0.23845253         0.00000000                0 0.00000000
    ## 2001-03-01   2.01794132         0.00000000                0 0.00000000
    ## 2001-04-01   0.49146494         0.00000000                0 0.00000000
    ## 2001-05-01   0.60312887         0.00000000                0 0.00000000
    ## 2001-06-01   0.16640364         0.00000000                0 0.00000000
    ## 2001-07-01   2.60039972         0.28906846                0 0.28906846
    ## 2001-08-01   1.56079752         0.04493462                0 0.04493462
    ## 2001-09-01   1.08679369         0.00000000                0 0.00000000
    ## 2001-10-01   0.98276002         0.11206467                0 0.11206467
    ## 2001-11-01   1.83385259         0.25110065                0 0.25110065
    ## 2001-12-01   0.02307311         0.00000000                0 0.00000000
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2001-01-01   1.50422719               0          0.8593736   0.45301063
    ## 2001-02-01   0.01261102               0          0.6125900   0.04764393
    ## 2001-03-01   1.24367981               0          1.0624895   0.41885712
    ## 2001-04-01   0.00000000               0          0.8289131   0.13909570
    ## 2001-05-01   0.07363909               0          1.0347856   0.16024458
    ## 2001-06-01   0.00000000               0          1.1186614   0.05319145
    ## 2001-07-01   0.12858424               0          1.4441387   0.38863773
    ## 2001-08-01   0.43505198               0          1.3644548   0.31649677
    ## 2001-09-01   0.15065663               0          1.0399323   0.21971934
    ## 2001-10-01   0.53858933               0          0.8490664   0.24115707
    ## 2001-11-01   1.31985118               0          0.5502709   0.12070957
    ## 2001-12-01   0.00000000               0          0.3988745   0.02738870
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01     0.184236164                 0       0.2221268     0.2221268
    ## 2001-02-01     0.065423904                 0       0.4995222     0.4995222
    ## 2001-03-01     0.124310900                 0       0.5193215     0.5193215
    ## 2001-04-01     0.016155926                 0       0.6736615     0.6736615
    ## 2001-05-01     0.075090922                 0       0.7994501     0.7994501
    ## 2001-06-01     0.006163322                 0       1.0593067     1.0593067
    ## 2001-07-01     0.102316292                 0       0.9531847     0.9531847
    ## 2001-08-01     0.071023070                 0       0.9769349     0.9769349
    ## 2001-09-01     0.063895690                 0       0.7563173     0.7563173
    ## 2001-10-01     0.080263992                 0       0.5276454     0.5276454
    ## 2001-11-01     0.053296432                 0       0.3762649     0.3762649
    ## 2001-12-01     0.019776793                 0       0.3517090     0.3517090
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2001-01-01                      0            0.0009234183
    ## 2001-02-01                      0            0.0000000000
    ## 2001-03-01                      0            0.0012544086
    ## 2001-04-01                      0            0.0005414322
    ## 2001-05-01                      0            0.0000000000
    ## 2001-06-01                      0            0.0000000000
    ## 2001-07-01                      0            0.0165760485
    ## 2001-08-01                      0            0.0037513475
    ## 2001-09-01                      0            0.0013013133
    ## 2001-10-01                      0            0.0014121583
    ## 2001-11-01                      0            0.0008617183
    ## 2001-12-01                      0            0.0000000000

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of plant cohorts by months:

``` r

summary(S, freq="months",FUN=mean, output="PlantStress")
```

    ##                 T1_158      T2_179      S1_176
    ## 2001-01-01 0.004774718 0.006369228 0.004636100
    ## 2001-02-01 0.006904397 0.007545434 0.005927707
    ## 2001-03-01 0.006409395 0.007232488 0.005602524
    ## 2001-04-01 0.012499747 0.009432663 0.008542719
    ## 2001-05-01 0.011048040 0.008844053 0.007780924
    ## 2001-06-01 0.121924036 0.018126328 0.027267909
    ## 2001-07-01 0.034436494 0.009829380 0.010339992
    ## 2001-08-01 0.008815849 0.008221264 0.006853804
    ## 2001-09-01 0.008362553 0.008134317 0.006680760
    ## 2001-10-01 0.006201583 0.007196549 0.005522331
    ## 2001-11-01 0.006887895 0.007464402 0.005870680
    ## 2001-12-01 0.009964233 0.008874725 0.007569394

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
    ## 2001-01-01      0.004613739       0.004541797  0.006282261
    ## 2001-01-02      0.004613739       0.004541797  0.006282261
    ## 2001-01-03      0.004613739       0.004541797  0.006282261
    ## 2001-01-04      0.004958814       0.004754072  0.006476584
    ## 2001-01-05      0.004613739       0.004541797  0.006282261
    ## 2001-01-06      0.004689850       0.004588090  0.006323880

Or we can combine the aggregation by species with a temporal aggregation
(here monthly averages):

``` r

summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01      0.004774718       0.004636100  0.006369228
    ## 2001-02-01      0.006904397       0.005927707  0.007545434
    ## 2001-03-01      0.006409395       0.005602524  0.007232488
    ## 2001-04-01      0.012499747       0.008542719  0.009432663
    ## 2001-05-01      0.011048040       0.007780924  0.008844053
    ## 2001-06-01      0.121924036       0.027267909  0.018126328
    ## 2001-07-01      0.034436494       0.010339992  0.009829380
    ## 2001-08-01      0.008815849       0.006853804  0.008221264
    ## 2001-09-01      0.008362553       0.006680760  0.008134317
    ## 2001-10-01      0.006201583       0.005522331  0.007196549
    ## 2001-11-01      0.006887895       0.005870680  0.007464402
    ## 2001-12-01      0.009964233       0.007569394  0.008874725

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

    ##            T1_158 T2_179 S1_176
    ## 2001-01-01      0      0      0

``` r

droughtStress(S, index = "MDS", freq = "years", draw=FALSE)
```

    ##               T1_158     T2_179     S1_176
    ## 2001-01-01 0.2782897 0.02441786 0.04453265

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
    ## 2001-01-01  10.103225
    ## 2001-02-01   8.193833
    ## 2001-03-01   8.986614
    ## 2001-04-01   8.786425
    ## 2001-05-01   8.234128
    ## 2001-06-01   6.403861
    ## 2001-07-01   7.008689
    ## 2001-08-01   6.282414
    ## 2001-09-01   7.444944
    ## 2001-10-01   7.642402
    ## 2001-11-01   8.718791
    ## 2001-12-01   8.398522

## References

- De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
  R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
  forest inventory data to predict drought stress: the role of forest
  structural changes vs. climate changes. Agricultural and Forest
  Meteorology 213: 77-90
  (<https://doi.org/10.1016/j.agrformet.2015.06.012>).
