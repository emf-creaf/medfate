# Advanced water and energy balance

## About this vignette

This document describes how to run a water and energy balance model that
uses a more detailed approach for hydraulics and stomatal regulation
described in De Cáceres et al. (2021) and Ruffault et al. (2022). We
recommend reading vignette [*Basic water
balance*](https://emf-creaf.github.io/medfate/articles/runmodels/BasicWaterBalance.html)
before this one for a more accessible introduction to soil water balance
modelling. This vignette is meant to teach users to run the simulation
model within R. All the details of the model design and formulation can
be found at the
[medfatebook](https://emf-creaf.github.io/medfatebook/index.html).

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

Apart from data inputs, the behaviour of simulation models is controlled
using a set of global parameters. The default parameterization is
obtained using function
[`defaultControl()`](https://emf-creaf.github.io/medfate/reference/defaultControl.md):

``` r

control <- defaultControl("Sperry")
```

To use the advanced water balance model we must change the values of
`transpirationMode` to switch from `"Granier"` to either `"Sperry"` or
`"Sureau"`.

Since we will be inspecting subdaily results, we need to set the flag to
obtain subdaily output:

``` r

control$subdailyResults <- TRUE
```

### Water balance input object

A last object is needed before calling simulation functions, called
`spwbInput`. It consists in the compilation of aboveground, belowground
parameters and the specification of additional parameter values for each
plant cohort. This is done by calling function
[`spwbInput()`](https://emf-creaf.github.io/medfate/reference/modelInput.md):

``` r

x <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
```

The `spwbInput` object for advanced water and energy balance is similar
to that of simple water balance simulations, but contains more elements.
Information about the cohort species is found in element `cohorts`,
i.e. the cohort code, the species index and species name:

``` r

x$cohorts
```

    ##         SP              Name
    ## T1_148 148  Pinus halepensis
    ## T2_168 168      Quercus ilex
    ## S1_165 165 Quercus coccifera

Element `soil` contains soil layer parameters and state variables
(moisture and temperature):

``` r

x$soil
```

    ##   widths sand clay      usda om nitrogen ph  bd rfc  macro     Ksat VG_alpha
    ## 1    300   25   25 Silt loam NA       NA NA 1.5  25 0.0485 5401.471 89.16112
    ## 2    700   25   25 Silt loam NA       NA NA 1.5  45 0.0485 5401.471 89.16112
    ## 3   1000   25   25 Silt loam NA       NA NA 1.5  75 0.0485 5401.471 89.16112
    ## 4   2000   25   25 Silt loam NA       NA NA 1.5  95 0.0485 5401.471 89.16112
    ##       VG_n VG_theta_res VG_theta_sat W Temp
    ## 1 1.303861        0.041     0.423715 1   NA
    ## 2 1.303861        0.041     0.423715 1   NA
    ## 3 1.303861        0.041     0.423715 1   NA
    ## 4 1.303861        0.041     0.423715 1   NA

As an aside, the columns in `x$soil` that were not present in the input
data frame `examplesoil` are created by an internal call to a soil
initialization function called
[`soil()`](https://emf-creaf.github.io/medfate/reference/soil.md).

Element `canopy` contains state variables within the canopy:

``` r

x$canopy
```

    ##    zlow zmid  zup LAIlive LAIexpanded LAIdead LAImistletoe Tair Cair VPair
    ## 1     0   50  100      NA          NA      NA           NA   NA   NA    NA
    ## 2   100  150  200      NA          NA      NA           NA   NA   NA    NA
    ## 3   200  250  300      NA          NA      NA           NA   NA   NA    NA
    ## 4   300  350  400      NA          NA      NA           NA   NA   NA    NA
    ## 5   400  450  500      NA          NA      NA           NA   NA   NA    NA
    ## 6   500  550  600      NA          NA      NA           NA   NA   NA    NA
    ## 7   600  650  700      NA          NA      NA           NA   NA   NA    NA
    ## 8   700  750  800      NA          NA      NA           NA   NA   NA    NA
    ## 9   800  850  900      NA          NA      NA           NA   NA   NA    NA
    ## 10  900  950 1000      NA          NA      NA           NA   NA   NA    NA
    ## 11 1000 1050 1100      NA          NA      NA           NA   NA   NA    NA
    ## 12 1100 1150 1200      NA          NA      NA           NA   NA   NA    NA
    ## 13 1200 1250 1300      NA          NA      NA           NA   NA   NA    NA
    ## 14 1300 1350 1400      NA          NA      NA           NA   NA   NA    NA
    ## 15 1400 1450 1500      NA          NA      NA           NA   NA   NA    NA
    ## 16 1500 1550 1600      NA          NA      NA           NA   NA   NA    NA
    ## 17 1600 1650 1700      NA          NA      NA           NA   NA   NA    NA
    ## 18 1700 1750 1800      NA          NA      NA           NA   NA   NA    NA
    ## 19 1800 1850 1900      NA          NA      NA           NA   NA   NA    NA
    ## 20 1900 1950 2000      NA          NA      NA           NA   NA   NA    NA
    ## 21 2000 2050 2100      NA          NA      NA           NA   NA   NA    NA
    ## 22 2100 2150 2200      NA          NA      NA           NA   NA   NA    NA
    ## 23 2200 2250 2300      NA          NA      NA           NA   NA   NA    NA
    ## 24 2300 2350 2400      NA          NA      NA           NA   NA   NA    NA
    ## 25 2400 2450 2500      NA          NA      NA           NA   NA   NA    NA
    ## 26 2500 2550 2600      NA          NA      NA           NA   NA   NA    NA
    ## 27 2600 2650 2700      NA          NA      NA           NA   NA   NA    NA
    ## 28 2700 2750 2800      NA          NA      NA           NA   NA   NA    NA

Canopy temperature, water vapour pressure and $`CO_2`$ concentration are
state variables needed for canopy energy balance. If the canopy energy
balance assumes a single canopy layer, the same values will be assumed
through the canopy. Variation of within-canopy state variables is
modelled if a multi-canopy energy balance is used (see control parameter
`multiLayerBalance`).

As you may already known, element `above` contains the aboveground
structure data that we already know:

``` r

x$above
```

    ##          H        CR   LAI_live LAI_expanded LAI_dead LAI_mistletoe Age ObsID
    ## T1_148 800 0.6534132 1.20935408   1.20935408        0             0  NA  <NA>
    ## T2_168 660 0.6359169 0.56790878   0.56790878        0             0  NA  <NA>
    ## S1_165  80 0.8032817 0.04877848   0.04877848        0             0  NA  <NA>

Belowground parameters can be seen in `below`:

``` r

x$below
```

    ##        Z50  Z95 Z100 poolProportions
    ## T1_148 100  300   NA      0.66228187
    ## T2_168 300 1000   NA      0.31100544
    ## S1_165 200 1000   NA      0.02671269

and in `belowLayers`:

``` r

x$belowLayers
```

    ## $V
    ##                1          2           3            4
    ## T1_148 0.9498377 0.04811006 0.001774047 0.0002781442
    ## T2_168 0.5008953 0.45059411 0.040648313 0.0078622840
    ## S1_165 0.6799879 0.27379114 0.035676316 0.0105446776
    ## 
    ## $L
    ##               1        2        3        4
    ## T1_148 2086.448 1307.358 2073.244 4045.856
    ## T2_168 1413.593 1807.433 2300.951 4209.276
    ## S1_165 1907.898 1779.971 2349.397 4300.342
    ## 
    ## $VGrhizo_kmax
    ##               1          2          3          4
    ## T1_148 12670486   641770.5   23665.14   3710.341
    ## T2_168 44953618 40439261.1 3648045.36 705612.759
    ## S1_165  9962311  4011237.2  522683.69 154487.113
    ## 
    ## $VCroot_kmax
    ##               1         2           3            4
    ## T1_148 2.889971 0.2336108 0.005432084 0.0004364266
    ## T2_168 2.024110 1.4240808 0.100912866 0.0106697119
    ## S1_165 2.361274 1.0190767 0.100606036 0.0162454141
    ## 
    ## $Wpool
    ##        1 2 3 4
    ## T1_148 1 1 1 1
    ## T2_168 1 1 1 1
    ## S1_165 1 1 1 1
    ## 
    ## $RhizoPsi
    ##             1      2      3      4
    ## T1_148 -0.033 -0.033 -0.033 -0.033
    ## T2_168 -0.033 -0.033 -0.033 -0.033
    ## S1_165 -0.033 -0.033 -0.033 -0.033
    ## 
    ## $RHOP
    ## $RHOP$T1_148
    ##                 1          2          3          4
    ## T1_148 0.66228187 0.66228187 0.66228187 0.66228187
    ## T2_168 0.31100544 0.31100544 0.31100544 0.31100544
    ## S1_165 0.02671269 0.02671269 0.02671269 0.02671269
    ## 
    ## $RHOP$T2_168
    ##                 1          2          3          4
    ## T1_148 0.66228187 0.66228187 0.66228187 0.66228187
    ## T2_168 0.31100544 0.31100544 0.31100544 0.31100544
    ## S1_165 0.02671269 0.02671269 0.02671269 0.02671269
    ## 
    ## $RHOP$S1_165
    ##                 1          2          3          4
    ## T1_148 0.66228187 0.66228187 0.66228187 0.66228187
    ## T2_168 0.31100544 0.31100544 0.31100544 0.31100544
    ## S1_165 0.02671269 0.02671269 0.02671269 0.02671269

The `spwbInput`object also includes cohort parameter values for several
kinds of traits. For example, plant anatomy parameters are described in
`paramsAnatomy`:

``` r

x$paramsAnatomy
```

    ##        Hmed    Al2As      SLA LeafWidth LeafDensity WoodDensity FineRootDensity
    ## T1_148  850 1317.523 5.140523 0.1384772   0.2881200        0.60       0.2881200
    ## T2_168   48 3908.823 6.340000 1.7674359   0.6650000        0.89       0.6650000
    ## S1_165   90 2436.475 7.152702 1.3200000   0.3669009        0.73       0.3669009
    ##        conduit2sapwood       SRL RLD     r635
    ## T1_148       0.9258273 3870.0000  10 1.964226
    ## T2_168       0.6263370  735.7025  10 1.805872
    ## S1_165       0.6263370  996.2760  10 2.289452

Parameters related to plant transpiration and photosynthesis can be seen
in `paramsTranspiration`:

``` r

x$paramsTranspiration
```

    ##             Gswmin    Gswmax  Vmax298  Jmax298 Kmax_stemxylem Kmax_rootxylem
    ## T1_148 0.003086667 0.2850000 72.19617 124.1687           0.15           0.60
    ## T2_168 0.004473333 0.2007222 68.51600 118.7863           0.40           1.60
    ## S1_165 0.010455247 0.1800000 58.79590 104.3488           0.29           1.16
    ##        VCleaf_kmax VCleafapo_kmax VCleaf_slope VCleaf_P50  VCleaf_c  VCleaf_d
    ## T1_148    4.000000          8.000    133.86620  -2.303772 11.137050 -2.380849
    ## T2_168    4.000000          8.000     19.14428  -1.964085  1.339370 -2.582279
    ## S1_165    6.763001         13.526     47.93469  -3.022976  4.412593 -3.284789
    ##        kleaf_symp VCstem_kmax VCstem_slope VCstem_P50  VCstem_c  VCstem_d
    ## T1_148      8.000    1.339563    68.302912  -5.139633 12.709999 -5.290000
    ## T2_168      8.000    4.992013    14.607857  -6.964747  3.560000 -7.720000
    ## S1_165     13.526    7.473393     9.773996  -7.000000  2.132799 -8.312465
    ##        VCroot_kmax VCroot_slope VCroot_P50  VCroot_c  VCroot_d VGrhizo_kmax
    ## T1_148    3.129450    103.96607  -2.966325 11.137050 -3.065569     13339632
    ## T2_168    3.559773     22.32794  -1.684034  1.339370 -2.214081     89746537
    ## S1_165    3.497202     16.46463  -4.389000  2.469094 -5.091343     14650719
    ##        Plant_kmax   FR_leaf   FR_stem   FR_root
    ## T1_148  0.7598454 0.1899614 0.5672339 0.2428047
    ## T2_168  1.3675461 0.3418865 0.2739468 0.3841667
    ## S1_165  1.7617598 0.2604997 0.2357376 0.5037627

Parameters related to pressure-volume curves and water storage capacity
of leaf and stem organs are in `paramsWaterStorage`:

``` r

x$paramsWaterStorage
```

    ##           maxFMC maxMCleaf maxMCstem   LeafPI0  LeafEPS    LeafAF     Vleaf
    ## T1_148 122.16441  142.2755 101.30719 -1.236151 10.52923 0.3447917 0.5488595
    ## T2_168 109.23042  159.3801  47.00007 -2.070000 18.27000 0.4580000 0.1347648
    ## S1_165  98.39177  132.9039  71.62682 -2.370000 17.23000 0.2400000 0.2902653
    ##        StemPI0  StemEPS    StemAF Vsapwood
    ## T1_148 -1.9760 12.82234 0.9258273 6.306451
    ## T2_168 -3.1824 44.33022 0.6263370 1.318710
    ## S1_165 -2.5168 22.41056 0.6263370 1.384287

Finally, remember that one can play with plant-specific parameters for
soil water balance (instead of using species-level values) by modifying
manually the parameter values in this object.

## Static analysis of sub-models

Before using the advanced water and energy balance model, is important
to understand the parameters that influence the different sub-models.
Package `medfate` provides low-level functions corresponding to
sub-models (light extinction, hydraulics, transpiration,
photosynthesis…). In addition, there are several high-level plotting
functions that allow examining several aspects of these processes.

### Vulnerability curves

Given a `spwbInput` object, we can use function
[`hydraulics_vulnerabilityCurvePlot()`](https://emf-creaf.github.io/medfate/reference/hydraulics_conductancefunctions.md)
to inspect **vulnerability curves** (i.e. how hydraulic conductance of a
given segment changes with the water potential) for each plant cohort
and each of the different segments of the soil-plant hydraulic network:
rhizosphere, roots, stems and leaves:

``` r

hydraulics_vulnerabilityCurvePlot(x, type="leaf")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-17-1.png)

``` r

hydraulics_vulnerabilityCurvePlot(x, type="stem")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-17-2.png)

``` r

hydraulics_vulnerabilityCurvePlot(x, type="root")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-17-3.png)

``` r

hydraulics_vulnerabilityCurvePlot(x, examplesoil, type="rhizo")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-17-4.png)

The maximum values and shape of vulnerability curves for leaves and
stems are regulated by parameters in `paramsTranspiration`. Roots have
vulnerability curve parameters in the same data frame, but maximum
conductance values need to be specified for each soil layer and are
given in `belowLayers$VCroot_kmax`. Note that the last call to
[`hydraulics_vulnerabilityCurvePlot()`](https://emf-creaf.github.io/medfate/reference/hydraulics_conductancefunctions.md)
includes a `soil` object. This is because the van Genuchten parameters
that define the shape of the vulnerability curve for the rhizosphere are
stored in this object. Maximum conductance values in the rhizosphere are
given in `belowLayers$VGrhizo_kmax`.

### Supply functions

The vulnerability curves conforming the hydraulic network are used in
the model to build the **supply function**, which relates water flow
(i.e. transpiration) with the drop of water potential along the whole
hydraulic pathway. The supply function contains not only these two
variables, but also the water potential of intermediate nodes in the the
hydraulic network. Function
[`hydraulics_supplyFunctionPlot()`](https://emf-creaf.github.io/medfate/reference/hydraulics_supplyfunctions.md)
can be used to inspect any of this variables:

``` r

hydraulics_supplyFunctionPlot(x, type="E")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-18-1.png)

``` r

hydraulics_supplyFunctionPlot(x, type="ERhizo")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-18-2.png)

``` r

hydraulics_supplyFunctionPlot(x, type="dEdP")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-18-3.png)

``` r

hydraulics_supplyFunctionPlot(x, type="StemPsi")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-18-4.png)

Calls to
[`hydraulics_supplyFunctionPlot()`](https://emf-creaf.github.io/medfate/reference/hydraulics_supplyfunctions.md)
always need both a `spwbInput` object and a `soil` object. The soil
moisture state (i.e. its water potential) is the starting point for the
calculation of the supply function, so different curves will be obtained
for different values of soil moisture.

### Pressure volume curves

``` r

moisture_pressureVolumeCurvePlot(x, segment="leaf", fraction="symplastic")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-20-1.png)

``` r

moisture_pressureVolumeCurvePlot(x, segment="leaf", fraction="apoplastic")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-20-2.png)

``` r

moisture_pressureVolumeCurvePlot(x, segment="stem", fraction="symplastic")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-20-3.png)

``` r

moisture_pressureVolumeCurvePlot(x, segment="stem", fraction="apoplastic")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-20-4.png)

## Water balance for a single day

### Running the model

Soil water balance simulations will normally span periods of several
months or years, but since the model operates at a daily and subdaily
temporal scales, it is possible to perform soil water balance for one
day only. This is done using function
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md).
In the following code we select the same day as before from the
meteorological input data and perform soil water balance for that day
only:

``` r

date <- examplemeteo$dates[d]
meteovec <- unlist(examplemeteo[d,])
sd1<-spwb_day(x, date, meteovec, 
              latitude = 41.82592, elevation = 100, slope= 0, aspect = 0)
```

The output of
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
is a list with several elements:

``` r

names(sd1)
```

    ##  [1] "cohorts"            "topography"         "weather"           
    ##  [4] "WaterBalance"       "EnergyBalance"      "Soil"              
    ##  [7] "Stand"              "Plants"             "SunlitLeaves"      
    ## [10] "ShadeLeaves"        "RhizoPsi"           "ExtractionInst"    
    ## [13] "PlantsInst"         "RadiationInputInst" "SunlitLeavesInst"  
    ## [16] "ShadeLeavesInst"    "LightExtinction"    "LWRExtinction"     
    ## [19] "CanopyTurbulence"

### Water balance output

Element `WaterBalance` contains the soil water balance flows of the day
(precipitation, infiltration, transpiration, …)

``` r

sd1$WaterBalance
```

    ##                     PET                    Rain                    Snow 
    ##                3.902334                0.000000                0.000000 
    ##                 NetRain                Snowmelt                   Runon 
    ##                0.000000                0.000000                0.000000 
    ##            Infiltration      InfiltrationExcess        SaturationExcess 
    ##                0.000000                0.000000                0.000000 
    ##                  Runoff            DeepDrainage         CapillarityRise 
    ##                0.000000                0.000000                0.000000 
    ##         SoilEvaporation       HerbTranspiration         PlantExtraction 
    ##                0.500000                0.000000                0.867061 
    ##           Transpiration  MistletoeTranspiration HydraulicRedistribution 
    ##                0.867061                0.000000                0.000000

And `Soil` contains water evaporated from each soil layer, water
transpired from each soil layer and the final soil water potential:

``` r

sd1$Soil
```

    ##           Psi HerbTranspiration HydraulicInput HydraulicOutput PlantExtraction
    ## 1 -0.03591934                 0              0    0.7199065403    0.7199065403
    ## 2 -0.03318611                 0              0    0.1384719421    0.1384719421
    ## 3 -0.03301617                 0              0    0.0078299284    0.0078299284
    ## 4 -0.03300440                 0              0    0.0008525875    0.0008525875

### Soil and canopy energy balance

Element `EnergyBalance` contains subdaily variation in atmosphere,
canopy and soil temperatures, as well as canopy and soil energy balance
components.

``` r

names(sd1$EnergyBalance)
```

    ## [1] "Temperature"         "SoilTemperature"     "CanopyEnergyBalance"
    ## [4] "SoilEnergyBalance"   "TemperatureLayers"   "VaporPressureLayers"

Package `medfate` provides a `plot` function for objects of class
`spwb_day` that can be used to inspect the results of the simulation. We
use this function to display subdaily dynamics in plant, soil and canopy
variables. For example, we can use it to display temperature variations
(only the temperature of the topmost soil layer is drawn):

``` r

plot(sd1, type = "Temperature")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-26-1.png)

``` r

plot(sd1, type = "CanopyEnergyBalance")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-26-2.png)

``` r

plot(sd1, type = "SoilEnergyBalance")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-26-3.png)

### Plant output

Element `Plants` contains output values by plant cohort. Several output
variables can be inspected in this element.

``` r

sd1$Plants
```

    ##               LAI    LAIlive     FPAR Extraction Transpiration
    ## T1_148 1.20935408 1.20935408 88.57701 0.62989100    0.62989100
    ## T2_168 0.56790878 0.56790878 69.05317 0.20563718    0.20563718
    ## S1_165 0.04877848 0.04877848 39.89406 0.03153282    0.03153282
    ##        MistletoeTranspiration GrossPhotosynthesis NetPhotosynthesis    RootPsi
    ## T1_148                      0          2.89337986        2.73690248 -0.3982251
    ## T2_168                      0          1.38524378        1.29983524 -0.3202356
    ## S1_165                      0          0.09581454        0.09031333 -0.5013209
    ##           StemPsi LeafPLC      StemPLC LeafPsiMin  LeafPsiMax      dEdP
    ## T1_148 -1.2513865       0 5.363683e-09 -1.5379962 -0.03981994 0.5124954
    ## T2_168 -0.5179479       0 3.645308e-05 -0.8080977 -0.03734174 0.8780433
    ## S1_165 -0.7210844       0 3.091533e-03 -0.9636222 -0.04140432 1.1873386
    ##                DDS   StemRWC   LeafRWC      LFMC  WaterBalance
    ## T1_148 0.001393765 0.9972162 0.9563039 118.86092 -6.288373e-17
    ## T2_168 0.049387480 0.9983946 0.9920787 108.49764  1.672382e-17
    ## S1_165 0.002168005 0.9937055 0.9852409  97.28107  1.355253e-19

While `Plants` contains one value per cohort and variable that
summarizes the whole simulated day, information by disaggregated by time
step can be accessed in `PlantsInst`. Moreover, we can use function
[`plot.spwb_day()`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md)
to draw plots of sub-daily variation per species of plant transpiration
per ground area (L·m$`^{-2}`$), transpiration per leaf area (also in
L·m$`^{-2}`$), plant net photosynthesis (in g C·m$`^{-2}`$), and plant
water potential (in MPa):

``` r

plot(sd1, type = "PlantTranspiration", bySpecies = T)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-28-1.png)

``` r

plot(sd1, type = "TranspirationPerLeaf", bySpecies = T)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-28-2.png)

``` r

plot(sd1, type = "NetPhotosynthesis", bySpecies = T)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-28-3.png)

``` r

plot(sd1, type = "LeafPsiAverage", bySpecies = T)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-28-4.png)

### Output for sunlit and shade leaves

The model distinguishes between sunlit and shade leaves for stomatal
regulation. Static properties of sunlit and shade leaves, for each
cohort, can be accessed via:

``` r

sd1$SunlitLeaves
```

    ##        LeafPsiMin  LeafPsiMax      GSWMin     GSWMax  TempMin  TempMax
    ## T1_148  -1.588172 -0.03981994 0.002240670 0.28134434 1.355752 11.52658
    ## T2_168  -1.217468 -0.03842744 0.003206471 0.08780195 1.353462 17.73773
    ## S1_165  -1.477458 -0.04560648 0.007894569 0.16149927 1.348698 16.88081

``` r

sd1$ShadeLeaves
```

    ##        LeafPsiMin  LeafPsiMax      GSWMin     GSWMax   TempMin  TempMax
    ## T1_148 -1.4377825 -0.03981994 0.002219031 0.28033228 1.0272810 9.882972
    ## T2_168 -0.5023063 -0.03734174 0.003228388 0.09948119 0.6445816 9.783669
    ## S1_165 -0.7945205 -0.04140432 0.007641817 0.17912181 0.7388108 9.488250

Instantaneous values are also stored for sunlit and shade leaves. We can
also use the `plot` function for objects of class `spwb_day` to draw
instantaneous variations in temperature for sunlit and shade leaves:

``` r

plot(sd1, type = "LeafTemperature", bySpecies=TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-30-1.png)

Note that sunlit leaves of some species reach temperatures higher than
the canopy. We can also plot variations in instantaneous gross and net
photosynthesis rates:

``` r

plot(sd1, type = "LeafGrossPhotosynthesis", bySpecies=TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-31-1.png)

``` r

plot(sd1, type = "LeafNetPhotosynthesis", bySpecies=TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-32-1.png)

Or variations in stomatal conductance:

``` r

plot(sd1, type = "LeafStomatalConductance", bySpecies=TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-33-1.png)

Or variations in vapour pressure deficit:

``` r

plot(sd1, type = "LeafVPD", bySpecies=TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-34-1.png)

Or variations in leaf water potential:

``` r

plot(sd1, type = "LeafPsi", bySpecies=TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-35-1.png)

``` r

plot(sd1, type = "LeafCi", bySpecies=TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-36-1.png)

``` r

plot(sd1, type = "LeafIntrinsicWUE", bySpecies=TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-37-1.png)

## Water balance for multiple days

### Running the model

Users will often use function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md) to run
the soil water balance model for several days. This function requires
the `spwbInput` object, the `soil` object and the meteorological data
frame. However, running
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
modified the input objects. In particular, the soil moisture at the end
of the simulation was:

``` r

x$soil$W
```

    ## [1] 0.9821059 0.9988130 0.9998966 0.9999719

And the temperature of soil layers:

``` r

x$soil$Temp
```

    ## [1] 7.649311 3.060544 2.380712 2.363290

We can also see the current state of canopy variables:

``` r

x$canopy
```

    ##    zlow zmid  zup    LAIlive LAIexpanded LAIdead LAImistletoe    Tair Cair
    ## 1     0   50  100 0.04877848  0.04877848       0            0 5.80097  386
    ## 2   100  150  200 0.00000000  0.00000000       0            0 5.80097  386
    ## 3   200  250  300 0.07498750  0.07498750       0            0 5.80097  386
    ## 4   300  350  400 0.32143491  0.32143491       0            0 5.80097  386
    ## 5   400  450  500 0.46098893  0.46098893       0            0 5.80097  386
    ## 6   500  550  600 0.45798803  0.45798803       0            0 5.80097  386
    ## 7   600  650  700 0.30777120  0.30777120       0            0 5.80097  386
    ## 8   700  750  800 0.15409228  0.15409228       0            0 5.80097  386
    ## 9   800  850  900 0.00000000  0.00000000       0            0 5.80097  386
    ## 10  900  950 1000 0.00000000  0.00000000       0            0 5.80097  386
    ## 11 1000 1050 1100 0.00000000  0.00000000       0            0 5.80097  386
    ## 12 1100 1150 1200 0.00000000  0.00000000       0            0 5.80097  386
    ## 13 1200 1250 1300 0.00000000  0.00000000       0            0 5.80097  386
    ## 14 1300 1350 1400 0.00000000  0.00000000       0            0 5.80097  386
    ## 15 1400 1450 1500 0.00000000  0.00000000       0            0 5.80097  386
    ## 16 1500 1550 1600 0.00000000  0.00000000       0            0 5.80097  386
    ## 17 1600 1650 1700 0.00000000  0.00000000       0            0 5.80097  386
    ## 18 1700 1750 1800 0.00000000  0.00000000       0            0 5.80097  386
    ## 19 1800 1850 1900 0.00000000  0.00000000       0            0 5.80097  386
    ## 20 1900 1950 2000 0.00000000  0.00000000       0            0 5.80097  386
    ## 21 2000 2050 2100 0.00000000  0.00000000       0            0 5.80097  386
    ## 22 2100 2150 2200 0.00000000  0.00000000       0            0 5.80097  386
    ## 23 2200 2250 2300 0.00000000  0.00000000       0            0 5.80097  386
    ## 24 2300 2350 2400 0.00000000  0.00000000       0            0 5.80097  386
    ## 25 2400 2450 2500 0.00000000  0.00000000       0            0 5.80097  386
    ## 26 2500 2550 2600 0.00000000  0.00000000       0            0 5.80097  386
    ## 27 2600 2650 2700 0.00000000  0.00000000       0            0 5.80097  386
    ## 28 2700 2750 2800 0.00000000  0.00000000       0            0 5.80097  386
    ##        VPair
    ## 1  0.5170718
    ## 2  0.5170718
    ## 3  0.5170718
    ## 4  0.5170718
    ## 5  0.5170718
    ## 6  0.5170718
    ## 7  0.5170718
    ## 8  0.5170718
    ## 9  0.5170718
    ## 10 0.5170718
    ## 11 0.5170718
    ## 12 0.5170718
    ## 13 0.5170718
    ## 14 0.5170718
    ## 15 0.5170718
    ## 16 0.5170718
    ## 17 0.5170718
    ## 18 0.5170718
    ## 19 0.5170718
    ## 20 0.5170718
    ## 21 0.5170718
    ## 22 0.5170718
    ## 23 0.5170718
    ## 24 0.5170718
    ## 25 0.5170718
    ## 26 0.5170718
    ## 27 0.5170718
    ## 28 0.5170718

We simply use function
[`resetInputs()`](https://emf-creaf.github.io/medfate/reference/resetInputs.md)
to reset state variables to their default values, so that the new
simulation is not affected by the end state of the previous simulation:

``` r

resetInputs(x)
x$soil$W
```

    ## [1] 1 1 1 1

``` r

x$soil$Temp
```

    ## [1] NA NA NA NA

``` r

x$canopy
```

    ##    zlow zmid  zup    LAIlive LAIexpanded LAIdead LAImistletoe Tair Cair VPair
    ## 1     0   50  100 0.04877848  0.04877848       0            0   NA   NA    NA
    ## 2   100  150  200 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 3   200  250  300 0.07498750  0.07498750       0            0   NA   NA    NA
    ## 4   300  350  400 0.32143491  0.32143491       0            0   NA   NA    NA
    ## 5   400  450  500 0.46098893  0.46098893       0            0   NA   NA    NA
    ## 6   500  550  600 0.45798803  0.45798803       0            0   NA   NA    NA
    ## 7   600  650  700 0.30777120  0.30777120       0            0   NA   NA    NA
    ## 8   700  750  800 0.15409228  0.15409228       0            0   NA   NA    NA
    ## 9   800  850  900 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 10  900  950 1000 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 11 1000 1050 1100 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 12 1100 1150 1200 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 13 1200 1250 1300 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 14 1300 1350 1400 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 15 1400 1450 1500 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 16 1500 1550 1600 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 17 1600 1650 1700 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 18 1700 1750 1800 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 19 1800 1850 1900 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 20 1900 1950 2000 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 21 2000 2050 2100 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 22 2100 2150 2200 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 23 2200 2250 2300 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 24 2300 2350 2400 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 25 2400 2450 2500 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 26 2500 2550 2600 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 27 2600 2650 2700 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 28 2700 2750 2800 0.00000000  0.00000000       0            0   NA   NA    NA

Now we are ready to call function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md):

``` r

S <- spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
```

    ## Initial plant water content (mm): 9.19483
    ## Initial soil water content (mm): 290.875
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2001]:............
    ## 
    ## Final plant water content (mm): 9.18601
    ## Final soil water content (mm): 269.068
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.00881253
    ## Plant water balance result (mm): -7.98572e-16
    ## Change in soil water content (mm): -21.8067
    ## Soil water balance result (mm): -21.8067
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): -7.10543e-15
    ## Water balance components:
    ##   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
    ##   Interception (mm) 98 Net rainfall (mm) 364
    ##   Infiltration (mm) 401 Infiltration excess (mm) 14 Saturation excess (mm) 0 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 17  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 277  Mistletoe transpiration (mm) 0
    ##   Plant extraction from soil (mm) 277  Plant water balance (mm) -0 Hydraulic redistribution (mm) 7
    ##   Runoff (mm) 14 Deep drainage (mm) 129

Function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md)
returns an object of class *spwb*. If we inspect its elements, we
realize that the output is arranged differently than in
[`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md):

``` r

names(S)
```

    ##  [1] "latitude"      "topography"    "weather"       "spwbInput"    
    ##  [5] "spwbOutput"    "WaterBalance"  "EnergyBalance" "Temperature"  
    ##  [9] "Soil"          "Snow"          "Stand"         "Plants"       
    ## [13] "SunlitLeaves"  "ShadeLeaves"   "subdaily"

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
    ## 2001-01-01    3.2972315                  0                0      0    2.8504389
    ## 2001-01-02    0.9721473                  0                0      0    0.4317983
    ## 2001-01-03    0.0000000                  0                0      0    0.0000000
    ## 2001-01-04    4.2353558                  0                0      0    3.1100106
    ## 2001-01-05    0.7332674                  0                0      0    0.1830706
    ## 2001-01-06   11.5940646                  0                0      0    4.1214138
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2001-01-01               0          2.0186700     1.571877       0.4387811
    ## 2001-01-02               0          2.0664933     1.526144       0.5000000
    ## 2001-01-03               0          0.9837776     0.000000       0.5000000
    ## 2001-01-04               0          1.7031847     1.561617       0.1319572
    ## 2001-01-05               0          1.7785187     1.151134       0.5000000
    ## 2001-01-06               0          2.2029054     1.765736       0.4309166
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01                 0     0.008011580   0.008011580
    ## 2001-01-02                 0     0.040349062   0.040349062
    ## 2001-01-03                 0     0.483777643   0.483777643
    ## 2001-01-04                 0     0.009610367   0.009610367
    ## 2001-01-05                 0     0.127384653   0.127384653
    ## 2001-01-06                 0     0.006252746   0.006252746
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2001-01-01                      0            0.000000e+00
    ## 2001-01-02                      0            0.000000e+00
    ## 2001-01-03                      0            0.000000e+00
    ## 2001-01-04                      0            0.000000e+00
    ## 2001-01-05                      0            0.000000e+00
    ## 2001-01-06                      0            1.789931e-05

Elements `Plants` is itself a list with several elements that contain
daily output results by plant cohorts, for example leaf minimum (midday)
water potentials are:

``` r

head(S$Plants$LeafPsiMin)
```

    ##               T1_148     T2_168     S1_165
    ## 2001-01-01 -1.202493 -0.4091069 -0.4884300
    ## 2001-01-02 -1.314571 -0.4159543 -0.6086803
    ## 2001-01-03 -1.277215 -0.5248507 -0.5884302
    ## 2001-01-04 -1.170085 -0.3684948 -0.4256515
    ## 2001-01-05 -1.345484 -0.4050873 -0.6048143
    ## 2001-01-06 -1.296085 -0.4054554 -0.4748413

### Plotting and summarizing results

Package `medfate` also provides a `plot` function for objects of class
`spwb`. It can be used to show the meteorological input. Additionally,
it can also be used to draw soil and plant variables. In the code below
we draw water fluxes, soil water potentials, plant transpiration and
plant (mid-day) water potential:

``` r

plot(S, type="Evapotranspiration")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-47-1.png)

``` r

plot(S, type="SoilPsi", bySpecies = TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-47-2.png)

``` r

plot(S, type="PlantTranspiration", bySpecies = TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-47-3.png)

``` r

plot(S, type="LeafPsiMin", bySpecies = TRUE)
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-47-4.png)

Alternatively, one can interactively create plots using function
`shinyplot`, e.g.:

``` r

shinyplot(S)
```

While the simulation model uses daily steps, users may be interested in
outputs at larger time scales. The package provides a `summary` for
objects of class `spwb`. This function can be used to summarize the
model’s output at different temporal steps (i.e. weekly, annual, …). For
example, to obtain the water balance by months one can use:

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
    ## 2001-07-01   2.57275281         0.20763417                0 0.20763417
    ## 2001-08-01   1.51362812         0.01204280                0 0.01204280
    ## 2001-09-01   1.03132466         0.00000000                0 0.00000000
    ## 2001-10-01   0.99133649         0.04440967                0 0.04440967
    ## 2001-11-01   1.84023816         0.20849338                0 0.20849338
    ## 2001-12-01   0.01963594         0.00000000                0 0.00000000
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2001-01-01  1.388467000               0          1.0136946   0.56270036
    ## 2001-02-01  0.027420692               0          0.6557828   0.05362304
    ## 2001-03-01  1.131270628               0          1.1881489   0.50458808
    ## 2001-04-01  0.028306087               0          0.9350160   0.16341945
    ## 2001-05-01  0.134832563               0          1.1705981   0.19275992
    ## 2001-06-01  0.000000000               0          0.9901624   0.06642325
    ## 2001-07-01  0.009788784               0          1.7506756   0.49771893
    ## 2001-08-01  0.065130310               0          1.7999564   0.39655800
    ## 2001-09-01  0.101262813               0          1.1964697   0.27518837
    ## 2001-10-01  0.218000234               0          1.0415246   0.30023559
    ## 2001-11-01  1.088762604               0          0.5686714   0.15693127
    ## 2001-12-01  0.000000000               0          0.5215887   0.03082587
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01     0.165948968                 0       0.2850452     0.2850452
    ## 2001-02-01     0.035170682                 0       0.5669891     0.5669891
    ## 2001-03-01     0.104890944                 0       0.5786698     0.5786698
    ## 2001-04-01     0.017867468                 0       0.7537291     0.7537291
    ## 2001-05-01     0.030107198                 0       0.9477310     0.9477310
    ## 2001-06-01     0.004151822                 0       0.9195874     0.9195874
    ## 2001-07-01     0.051227419                 0       1.2017293     1.2017293
    ## 2001-08-01     0.023978023                 0       1.3794203     1.3794203
    ## 2001-09-01     0.025292032                 0       0.8959893     0.8959893
    ## 2001-10-01     0.029876046                 0       0.7114130     0.7114130
    ## 2001-11-01     0.055803710                 0       0.3559365     0.3559365
    ## 2001-12-01     0.014760176                 0       0.4760026     0.4760026
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2001-01-01                      0            9.667711e-07
    ## 2001-02-01                      0            2.171704e-04
    ## 2001-03-01                      0            8.932041e-04
    ## 2001-04-01                      0            6.962279e-03
    ## 2001-05-01                      0            2.260135e-02
    ## 2001-06-01                      0            1.492290e-01
    ## 2001-07-01                      0            1.900298e-02
    ## 2001-08-01                      0            3.986888e-03
    ## 2001-09-01                      0            1.154380e-02
    ## 2001-10-01                      0            1.268837e-02
    ## 2001-11-01                      0            5.596708e-03
    ## 2001-12-01                      0            1.671150e-03

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of plant cohorts by months:

``` r

summary(S, freq="months",FUN=mean, output="PlantStress")
```

    ##                   T1_148     T2_168       S1_165
    ## 2001-01-01  0.0002064360 0.01712581 0.0005151787
    ## 2001-02-01  0.0007221917 0.03096680 0.0016919680
    ## 2001-03-01  0.0009343046 0.03579608 0.0018792397
    ## 2001-04-01  0.0021278910 0.04433516 0.0031890216
    ## 2001-05-01  0.0171716204 0.06234330 0.0088845504
    ## 2001-06-01  0.1974574478 0.10755189 0.1343658062
    ## 2001-07-01  0.0224613800 0.08631355 0.0430472187
    ## 2001-08-01  0.0066600310 0.08365915 0.0254609850
    ## 2001-09-01  0.0049048454 0.05738261 0.0178187014
    ## 2001-10-01 -0.0297394746 0.04543140 0.0089859504
    ## 2001-11-01 -0.0141330959 0.02643787 0.0030882595
    ## 2001-12-01  0.0007683433 0.02611106 0.0023722229

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
    ## 2001-01-01     7.427610e-05      0.0002435350   0.01628465
    ## 2001-01-02     2.061980e-04      0.0004090454   0.01581267
    ## 2001-01-03     2.275427e-04      0.0005729462   0.02163518
    ## 2001-01-04     6.350982e-05      0.0003199936   0.01416743
    ## 2001-01-05     2.722371e-04      0.0003704242   0.01512584
    ## 2001-01-06     1.443194e-04      0.0003449075   0.01325645

Or we can combine the aggregation by species with a temporal aggregation
(here monthly averages):

``` r

summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01     0.0002064360      0.0005151787   0.01712581
    ## 2001-02-01     0.0007221917      0.0016919680   0.03096680
    ## 2001-03-01     0.0009343046      0.0018792397   0.03579608
    ## 2001-04-01     0.0021278910      0.0031890216   0.04433516
    ## 2001-05-01     0.0171716204      0.0088845504   0.06234330
    ## 2001-06-01     0.1974574478      0.1343658062   0.10755189
    ## 2001-07-01     0.0224613800      0.0430472187   0.08631355
    ## 2001-08-01     0.0066600310      0.0254609850   0.08365915
    ## 2001-09-01     0.0049048454      0.0178187014   0.05738261
    ## 2001-10-01    -0.0297394746      0.0089859504   0.04543140
    ## 2001-11-01    -0.0141330959      0.0030882595   0.02643787
    ## 2001-12-01     0.0007683433      0.0023722229   0.02611106

## References

- De Cáceres M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L,
  Poyatos R, Cabon A, Granda V, Forner A, Valladares F, Martínez-Vilalta
  J (2021) Unravelling the effect of species mixing on water use and
  drought stress in holm oak forests: a modelling approach. Agricultural
  and Forest Meteorology 296
  (<https://doi.org/10.1016/j.agrformet.2020.108233>).

- Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022)
  SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations
  of plant water status and drought-induced mortality at the ecosystem
  level. Geoscientific Model Development 15, 5593-5626
  (<https://doi.org/10.5194/gmd-15-5593-2022>).
