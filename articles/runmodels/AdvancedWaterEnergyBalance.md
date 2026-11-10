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
    ## T1_158 158  Pinus halepensis
    ## T2_179 179      Quercus ilex
    ## S1_176 176 Quercus coccifera

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
    ## T1_158 800 0.6534132 0.82389823   0.82389823        0             0  NA  <NA>
    ## T2_179 660 0.6359169 0.62107792   0.62107792        0             0  NA  <NA>
    ## S1_176  80 0.8032817 0.05274013   0.05274013        0             0  NA  <NA>

Belowground parameters can be seen in `below`:

``` r

x$below
```

    ##        Z50  Z95 Z100 poolProportions
    ## T1_158 100  300   NA       0.5501030
    ## T2_179 300 1000   NA       0.4146833
    ## S1_176 200 1000   NA       0.0352137

and in `belowLayers`:

``` r

x$belowLayers
```

    ## $V
    ##                1          2           3            4
    ## T1_158 0.9498377 0.04811006 0.001774047 0.0002781442
    ## T2_179 0.5008953 0.45059411 0.040648313 0.0078622840
    ## S1_176 0.6799879 0.27379114 0.035676316 0.0105446776
    ## 
    ## $L
    ##                1        2        3        4
    ## T1_158  3578.714 1564.102 2134.426 4084.160
    ## T2_179 12204.805 9631.822 5217.298 6237.250
    ## S1_176  1907.898 1779.971 2349.397 4300.342
    ## 
    ## $VGrhizo_kmax
    ##               1          2           3           4
    ## T1_158  2036332   103141.9    3803.334    596.3061
    ## T2_179 46380599 41722941.0 3763846.755 728011.3127
    ## S1_176 34771884 14000593.4 1824345.335 539213.0079
    ## 
    ## $VCroot_kmax
    ##               1         2           3            4
    ## T1_158 3.097495 0.3589707 0.009700014 0.0007947957
    ## T2_179 1.606734 1.8314915 0.305017448 0.0493495758
    ## S1_176 2.361274 1.0190767 0.100606036 0.0162454141
    ## 
    ## $Wpool
    ##        1 2 3 4
    ## T1_158 1 1 1 1
    ## T2_179 1 1 1 1
    ## S1_176 1 1 1 1
    ## 
    ## $RhizoPsi
    ##             1      2      3      4
    ## T1_158 -0.033 -0.033 -0.033 -0.033
    ## T2_179 -0.033 -0.033 -0.033 -0.033
    ## S1_176 -0.033 -0.033 -0.033 -0.033
    ## 
    ## $RHOP
    ## $RHOP$T1_158
    ##                1         2         3         4
    ## T1_158 0.5501030 0.5501030 0.5501030 0.5501030
    ## T2_179 0.4146833 0.4146833 0.4146833 0.4146833
    ## S1_176 0.0352137 0.0352137 0.0352137 0.0352137
    ## 
    ## $RHOP$T2_179
    ##                1         2         3         4
    ## T1_158 0.5501030 0.5501030 0.5501030 0.5501030
    ## T2_179 0.4146833 0.4146833 0.4146833 0.4146833
    ## S1_176 0.0352137 0.0352137 0.0352137 0.0352137
    ## 
    ## $RHOP$S1_176
    ##                1         2         3         4
    ## T1_158 0.5501030 0.5501030 0.5501030 0.5501030
    ## T2_179 0.4146833 0.4146833 0.4146833 0.4146833
    ## S1_176 0.0352137 0.0352137 0.0352137 0.0352137

The `spwbInput`object also includes cohort parameter values for several
kinds of traits. For example, plant anatomy parameters are described in
`paramsAnatomy`:

``` r

x$paramsAnatomy
```

    ##        Hmed    Al2As      SLA LeafWidth LeafDensity WoodDensity FineRootDensity
    ## T1_158  880 1982.660 4.990020      0.10   0.2881200        0.60       0.3671652
    ## T2_179  550 1108.262 5.580754      1.63   0.6650000        0.89       0.4276170
    ## S1_176   90 2436.475 7.152702      1.32   0.3669009        0.73       0.4276170
    ##        conduit2sapwood       SRL RLD     r635
    ## T1_158       0.9258273 2115.0874  10 1.964226
    ## T2_179       0.6263370  735.7025  10 1.805872
    ## S1_176       0.6263370  996.2760  10 2.289452

Parameters related to plant transpiration and photosynthesis can be seen
in `paramsTranspiration`:

``` r

x$paramsTranspiration
```

    ##              Gswmin     Gswmax  Vmax298   Jmax298 Kmax_stemxylem Kmax_rootxylem
    ## T1_158 0.0006655459 0.04554533 69.30208 119.93960      0.4152667       1.661067
    ## T2_179 0.0028095533 0.21575000 39.44837  93.30728      0.8268000       3.307200
    ## S1_176 0.0104552469 0.17847825 58.79590 104.34879      0.2900000       1.160000
    ##        VCleaf_kmax VCleaf_slope VCleaf_P50 VCleaf_c  VCleaf_d VCleafapo_kmax
    ## T1_158    1.261346     74.72896  -1.939080 4.412593 -2.107019             NA
    ## T2_179    5.334302     41.40155  -3.500000 4.412593 -3.803127             NA
    ## S1_176    2.191344     47.93469  -3.022976 4.412593 -3.284789             NA
    ##        VCleafapo_slope VCleafapo_P50 VCleafapo_c VCleafapo_d kleaf_symp
    ## T1_158        43.65232         -5.14    5.760425   -5.477666         NA
    ## T2_179        21.11111         -6.88    8.087529   -7.198963         NA
    ## S1_176        12.20615         -7.00    3.071843   -7.887063         NA
    ##        VCstem_kmax VCstem_slope VCstem_P50 VCstem_c  VCstem_d VCroot_kmax
    ## T1_158    2.423697     43.65232      -5.14 5.760425 -5.477666    3.466961
    ## T2_179   11.288633     21.11111      -6.88 8.087529 -7.198963    3.792592
    ## S1_176    7.473393     12.20615      -7.00 3.071843 -7.887063    3.497202
    ##        VCroot_slope VCroot_P50 VCroot_c   VCroot_d VGrhizo_kmax Plant_kmax
    ## T1_158    102.98940     -0.880 3.193046 -0.9870361      2143874  0.6694187
    ## T2_179     11.52311     -2.390 1.036296 -3.4040609     92595398  1.8528039
    ## S1_176     16.46463     -4.389 2.469094 -5.0913431     51136035  1.1414327
    ##          FR_leaf   FR_stem   FR_root
    ## T1_158 0.5307175 0.2761973 0.1930852
    ## T2_179 0.3473376 0.1641301 0.4885323
    ## S1_176 0.5208826 0.1527329 0.3263846

Parameters related to pressure-volume curves and water storage capacity
of leaf and stem organs are in `paramsWaterStorage`:

``` r

x$paramsWaterStorage
```

    ##           maxFMC maxMCleaf maxMCstem   LeafPI0  LeafEPS    LeafAF     Vleaf
    ## T1_158 122.16441  142.2755 101.30719 -1.236151 10.52923 0.3447917 0.5654136
    ## T2_179 109.23042  159.3801  47.00007 -2.070000 18.27000 0.4580000 0.1530992
    ## S1_176  98.39177  132.9039  71.62682 -3.224891 17.23000 0.2400000 0.2902653
    ##        StemPI0  StemEPS    StemAF Vsapwood
    ## T1_158 -1.9760 12.82234 0.9258273 7.421896
    ## T2_179 -3.1824 44.33022 0.6263370 8.648152
    ## S1_176 -2.5168 22.41056 0.6263370 3.649952

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
    ##               3.9023342               0.0000000               0.0000000 
    ##                 NetRain                Snowmelt                   Runon 
    ##               0.0000000               0.0000000               0.0000000 
    ##            Infiltration      InfiltrationExcess        SaturationExcess 
    ##               0.0000000               0.0000000               0.0000000 
    ##                  Runoff            DeepDrainage         CapillarityRise 
    ##               0.0000000               0.0000000               0.0000000 
    ##         SoilEvaporation       HerbTranspiration         PlantExtraction 
    ##               0.5000000               0.0000000               0.5647084 
    ##           Transpiration  MistletoeTranspiration HydraulicRedistribution 
    ##               0.5647084               0.0000000               0.0000000

And `Soil` contains water evaporated from each soil layer, water
transpired from each soil layer and the final soil water potential:

``` r

sd1$Soil
```

    ##           Psi HerbTranspiration HydraulicInput HydraulicOutput PlantExtraction
    ## 1 -0.03492670                 0              0     0.316639704     0.316639704
    ## 2 -0.03328340                 0              0     0.210553087     0.210553087
    ## 3 -0.03306677                 0              0     0.032318515     0.032318515
    ## 4 -0.03302683                 0              0     0.005197109     0.005197109

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

    ##               LAI    LAIlive Phenology     FPAR Extraction Transpiration
    ## T1_158 0.82389823 0.82389823         1 92.06855 0.14834387    0.14834387
    ## T2_179 0.62107792 0.62107792         1 75.49574 0.38574189    0.38574189
    ## S1_176 0.05274013 0.05274013         1 46.97224 0.03062266    0.03062266
    ##        MistletoeTranspiration GrossPhotosynthesis NetPhotosynthesis    RootPsi
    ## T1_158                      0           1.8248077         1.7267856 -0.1744177
    ## T2_179                      0           1.2351319         1.1814757 -0.5408952
    ## S1_176                      0           0.1143856         0.1080716 -0.4358687
    ##           StemPsi      LeafPLC      StemPLC LeafPsiMin  LeafPsiMax      dEdP
    ## T1_158 -0.3763572 6.038646e-06 1.012553e-07 -0.7659753 -0.03521198 0.4508832
    ## T2_179 -0.6989382 7.742359e-08 3.263293e-09 -1.0339301 -0.03460471 1.2181082
    ## S1_176 -0.6243011 2.067046e-03 2.347408e-04 -1.2702230 -0.04532272 0.7679645
    ##                DDS   StemRWC   LeafRWC      LFMC  WaterBalance
    ## T1_158 0.002768963 0.9992345 0.9839824 120.96613  2.439455e-18
    ## T2_179 0.026612066 0.9980299 0.9904052 108.34229 -4.201283e-19
    ## S1_176 0.003860019 0.9958905 0.9804248  97.08963  7.318365e-19

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

    ##        LeafPsiMin  LeafPsiMax       GSWMin     GSWMax  TempMin  TempMax
    ## T1_158 -0.8162142 -0.03521198 0.0004840174 0.04527644 1.213131 11.41738
    ## T2_179 -1.3855247 -0.03540713 0.0020104459 0.12209879 1.211720 17.89513
    ## S1_176 -1.6726845 -0.05245694 0.0076069976 0.15292880 1.206226 17.05683

``` r

sd1$ShadeLeaves
```

    ##        LeafPsiMin  LeafPsiMax       GSWMin     GSWMax   TempMin   TempMax
    ## T1_158 -0.7157860 -0.03521198 0.0005038862 0.04551875 0.4836943  9.778602
    ## T2_179 -0.6513587 -0.03460471 0.0020294114 0.12470652 0.1038145 10.045358
    ## S1_176 -1.0696775 -0.04532272 0.0074749994 0.17002951 0.6442540 10.069963

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

    ## [1] 0.9880212 0.9981950 0.9995733 0.9998285

And the temperature of soil layers:

``` r

x$soil$Temp
```

    ## [1]  8.441068  5.216682  7.119441 11.721633

We can also see the current state of canopy variables:

``` r

x$canopy
```

    ##    zlow zmid  zup    LAIlive LAIexpanded LAIdead LAImistletoe     Tair Cair
    ## 1     0   50  100 0.05274013  0.05274013       0            0 5.477382  386
    ## 2   100  150  200 0.00000000  0.00000000       0            0 5.094773  386
    ## 3   200  250  300 0.07129816  0.07129816       0            0 4.712165  386
    ## 4   300  350  400 0.27802215  0.27802215       0            0 4.329556  386
    ## 5   400  450  500 0.38952510  0.38952510       0            0 3.946947  386
    ## 6   500  550  600 0.37113716  0.37113716       0            0 3.564338  386
    ## 7   600  650  700 0.23001492  0.23001492       0            0 3.181730  386
    ## 8   700  750  800 0.10497865  0.10497865       0            0 2.799121  386
    ## 9   800  850  900 0.00000000  0.00000000       0            0 2.799121  386
    ## 10  900  950 1000 0.00000000  0.00000000       0            0 2.799121  386
    ## 11 1000 1050 1100 0.00000000  0.00000000       0            0 2.799121  386
    ## 12 1100 1150 1200 0.00000000  0.00000000       0            0 2.799121  386
    ## 13 1200 1250 1300 0.00000000  0.00000000       0            0 2.799121  386
    ## 14 1300 1350 1400 0.00000000  0.00000000       0            0 2.799121  386
    ## 15 1400 1450 1500 0.00000000  0.00000000       0            0 2.799121  386
    ## 16 1500 1550 1600 0.00000000  0.00000000       0            0 2.799121  386
    ## 17 1600 1650 1700 0.00000000  0.00000000       0            0 2.799121  386
    ## 18 1700 1750 1800 0.00000000  0.00000000       0            0 2.799121  386
    ## 19 1800 1850 1900 0.00000000  0.00000000       0            0 2.799121  386
    ## 20 1900 1950 2000 0.00000000  0.00000000       0            0 2.799121  386
    ## 21 2000 2050 2100 0.00000000  0.00000000       0            0 2.799121  386
    ## 22 2100 2150 2200 0.00000000  0.00000000       0            0 2.799121  386
    ## 23 2200 2250 2300 0.00000000  0.00000000       0            0 2.799121  386
    ## 24 2300 2350 2400 0.00000000  0.00000000       0            0 2.799121  386
    ## 25 2400 2450 2500 0.00000000  0.00000000       0            0 2.799121  386
    ## 26 2500 2550 2600 0.00000000  0.00000000       0            0 2.799121  386
    ## 27 2600 2650 2700 0.00000000  0.00000000       0            0 2.799121  386
    ## 28 2700 2750 2800 0.00000000  0.00000000       0            0 2.799121  386
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
    ## 1     0   50  100 0.05274013  0.05274013       0            0   NA   NA    NA
    ## 2   100  150  200 0.00000000  0.00000000       0            0   NA   NA    NA
    ## 3   200  250  300 0.07129816  0.07129816       0            0   NA   NA    NA
    ## 4   300  350  400 0.27802215  0.27802215       0            0   NA   NA    NA
    ## 5   400  450  500 0.38952510  0.38952510       0            0   NA   NA    NA
    ## 6   500  550  600 0.37113716  0.37113716       0            0   NA   NA    NA
    ## 7   600  650  700 0.23001492  0.23001492       0            0   NA   NA    NA
    ## 8   700  750  800 0.10497865  0.10497865       0            0   NA   NA    NA
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

    ## Initial plant water content (mm): 12.2513
    ## Initial soil water content (mm): 290.875
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2001]:............
    ## 
    ## Final plant water content (mm): 12.2492
    ## Final soil water content (mm): 270.014
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.00215035
    ## Plant water balance result (mm): -1.06282e-15
    ## Change in soil water content (mm): -20.8606
    ## Soil water balance result (mm): -20.8606
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): 0
    ## Water balance components:
    ##   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
    ##   Interception (mm) 79 Net rainfall (mm) 383 Runon (mm) 0
    ##   Infiltration (mm) 412 Infiltration excess (mm) 22 Saturation excess (mm) 0 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 32  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 204  Mistletoe transpiration (mm) 0
    ##   Plant extraction from soil (mm) 204  Plant water balance (mm) -0 Hydraulic redistribution (mm) 4
    ##   Runoff (mm) 22 Deep drainage (mm) 197

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

    ##                  PET Precipitation      Rain Snow    NetRain Snowmelt Runon
    ## 2001-01-01 0.8828475      4.869109  4.869109    0  3.6578430        0     0
    ## 2001-01-02 1.6375337      2.498292  2.498292    0  1.3026155        0     0
    ## 2001-01-03 1.3017026      0.000000  0.000000    0  0.0000000        0     0
    ## 2001-01-04 0.5690790      5.796973  5.796973    0  4.5956523        0     0
    ## 2001-01-05 1.6760567      1.884401  1.884401    0  0.8616218        0     0
    ## 2001-01-06 1.2077028     13.359801 13.359801    0 11.9754283        0     0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2001-01-01    3.6578430                  0                0      0    3.1616586
    ## 2001-01-02    1.3026155                  0                0      0    0.7490317
    ## 2001-01-03    0.0000000                  0                0      0    0.0000000
    ## 2001-01-04    4.5956523                  0                0      0    3.6924806
    ## 2001-01-05    0.8616218                  0                0      0    0.2781413
    ## 2001-01-06   11.9754283                  0                0      0    4.1214138
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2001-01-01               0          1.7074502     1.211266       0.4944700
    ## 2001-01-02               0          1.7492599     1.195676       0.5000000
    ## 2001-01-03               0          0.7256467     0.000000       0.5000000
    ## 2001-01-04               0          1.3788456     1.201321       0.1750146
    ## 2001-01-05               0          1.6062602     1.022780       0.5000000
    ## 2001-01-06               0          1.8849710     1.384372       0.5000000
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01                 0    0.0017143936  0.0017143936
    ## 2001-01-02                 0    0.0535837943  0.0535837943
    ## 2001-01-03                 0    0.2256467448  0.2256467448
    ## 2001-01-04                 0    0.0025102369  0.0025102369
    ## 2001-01-05                 0    0.0834805865  0.0834805865
    ## 2001-01-06                 0    0.0005986964  0.0005986964
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2001-01-01                      0                       0
    ## 2001-01-02                      0                       0
    ## 2001-01-03                      0                       0
    ## 2001-01-04                      0                       0
    ## 2001-01-05                      0                       0
    ## 2001-01-06                      0                       0

Elements `Plants` is itself a list with several elements that contain
daily output results by plant cohorts, for example leaf minimum (midday)
water potentials are:

``` r

head(S$Plants$LeafPsiMin)
```

    ##                T1_158     T2_179     S1_176
    ## 2001-01-01 -0.4321934 -0.6075626 -0.7797777
    ## 2001-01-02 -0.4458637 -0.4985138 -0.9374211
    ## 2001-01-03 -0.4602478 -0.6237631 -0.8457173
    ## 2001-01-04 -0.3511891 -0.4589976 -0.6905754
    ## 2001-01-05 -0.4241627 -0.5886844 -0.9246433
    ## 2001-01-06 -0.3746708 -0.5034911 -0.7748204

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
    ##            Runon Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2001-01-01     0   1.84350049         0.00000000                0 0.00000000
    ## 2001-02-01     0   0.23845253         0.00000000                0 0.00000000
    ## 2001-03-01     0   2.01794132         0.00000000                0 0.00000000
    ## 2001-04-01     0   0.49146494         0.00000000                0 0.00000000
    ## 2001-05-01     0   0.60312887         0.00000000                0 0.00000000
    ## 2001-06-01     0   0.16640364         0.00000000                0 0.00000000
    ## 2001-07-01     0   2.57606293         0.31340525                0 0.31340525
    ## 2001-08-01     0   1.56188967         0.04384247                0 0.04384247
    ## 2001-09-01     0   1.08679369         0.00000000                0 0.00000000
    ## 2001-10-01     0   0.98664682         0.10817787                0 0.10817787
    ## 2001-11-01     0   1.83052339         0.25442985                0 0.25442985
    ## 2001-12-01     0   0.02307311         0.00000000                0 0.00000000
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2001-01-01   1.57168289               0          0.7773293   0.45301063
    ## 2001-02-01   0.04696911               0          0.4560936   0.04764393
    ## 2001-03-01   1.46135211               0          0.9021566   0.41885712
    ## 2001-04-01   0.13624961               0          0.6251960   0.13909570
    ## 2001-05-01   0.25738910               0          0.9171789   0.16024458
    ## 2001-06-01   0.00000000               0          1.1811524   0.05319145
    ## 2001-07-01   0.80599183               0          1.4215261   0.38863773
    ## 2001-08-01   0.50076769               0          1.4612622   0.31649677
    ## 2001-09-01   0.34688048               0          0.9457869   0.21971934
    ## 2001-10-01   0.40972350               0          0.8444113   0.24115707
    ## 2001-11-01   0.85052200               0          0.4709287   0.12070957
    ## 2001-12-01   0.00000000               0          0.3298811   0.02738870
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01     0.201047772                 0       0.1232709     0.1232709
    ## 2001-02-01     0.080827451                 0       0.3276222     0.3276222
    ## 2001-03-01     0.150082164                 0       0.3332173     0.3332173
    ## 2001-04-01     0.052137042                 0       0.4339633     0.4339633
    ## 2001-05-01     0.089982145                 0       0.6669522     0.6669522
    ## 2001-06-01     0.007583687                 0       1.1203773     1.1203773
    ## 2001-07-01     0.114408807                 0       0.9184796     0.9184796
    ## 2001-08-01     0.081013939                 0       1.0637515     1.0637515
    ## 2001-09-01     0.092071589                 0       0.6339960     0.6339960
    ## 2001-10-01     0.089197418                 0       0.5140568     0.5140568
    ## 2001-11-01     0.071398785                 0       0.2788203     0.2788203
    ## 2001-12-01     0.024626536                 0       0.2778658     0.2778658
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2001-01-01                      0            3.992254e-06
    ## 2001-02-01                      0            1.765603e-04
    ## 2001-03-01                      0            5.569043e-04
    ## 2001-04-01                      0            1.141370e-03
    ## 2001-05-01                      0            3.188577e-03
    ## 2001-06-01                      0            4.570761e-02
    ## 2001-07-01                      0            2.122542e-02
    ## 2001-08-01                      0            1.597728e-02
    ## 2001-09-01                      0            2.316096e-02
    ## 2001-10-01                      0            2.528167e-02
    ## 2001-11-01                      0            6.476684e-03
    ## 2001-12-01                      0            2.488413e-03

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of plant cohorts by months:

``` r

summary(S, freq="months",FUN=mean, output="PlantStress")
```

    ##                 T1_158     T2_179      S1_176
    ## 2001-01-01 0.000368473 0.01183003 0.000778966
    ## 2001-02-01 0.001988338 0.01778528 0.003554299
    ## 2001-03-01 0.001465891 0.01980082 0.004309948
    ## 2001-04-01 0.001685023 0.02307740 0.005536437
    ## 2001-05-01 0.004562799 0.03112796 0.011621784
    ## 2001-06-01 0.020176510 0.05852703 0.048695029
    ## 2001-07-01 0.008080109 0.04843653 0.028706490
    ## 2001-08-01 0.011486104 0.05022502 0.028594945
    ## 2001-09-01 0.004619614 0.03719911 0.017533303
    ## 2001-10-01 0.002715559 0.03176043 0.010120853
    ## 2001-11-01 0.001165185 0.01887562 0.004987009
    ## 2001-12-01 0.001900284 0.01680754 0.002640939

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
    ## 2001-01-01     0.0002749076      0.0045775452   0.01302928
    ## 2001-01-02     0.0002859634      0.0019255623   0.01130166
    ## 2001-01-03     0.0003361710      0.0006590488   0.01389147
    ## 2001-01-04     0.0001908709      0.0002522494   0.01135001
    ## 2001-01-05     0.0002600548      0.0006319488   0.01237614
    ## 2001-01-06     0.0001847598      0.0003874138   0.01120844

Or we can combine the aggregation by species with a temporal aggregation
(here monthly averages):

``` r

summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01      0.000368473       0.000778966   0.01183003
    ## 2001-02-01      0.001988338       0.003554299   0.01778528
    ## 2001-03-01      0.001465891       0.004309948   0.01980082
    ## 2001-04-01      0.001685023       0.005536437   0.02307740
    ## 2001-05-01      0.004562799       0.011621784   0.03112796
    ## 2001-06-01      0.020176510       0.048695029   0.05852703
    ## 2001-07-01      0.008080109       0.028706490   0.04843653
    ## 2001-08-01      0.011486104       0.028594945   0.05022502
    ## 2001-09-01      0.004619614       0.017533303   0.03719911
    ## 2001-10-01      0.002715559       0.010120853   0.03176043
    ## 2001-11-01      0.001165185       0.004987009   0.01887562
    ## 2001-12-01      0.001900284       0.002640939   0.01680754

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
