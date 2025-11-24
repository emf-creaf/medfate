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

    ##   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha
    ## 1    300   25   25 Silt loam NA       NA 1.5  25 0.0485 5401.471 89.16112
    ## 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112
    ## 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112
    ## 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112
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

    ##    zlow zmid  zup LAIlive LAIexpanded LAIdead Tair Cair VPair
    ## 1     0   50  100      NA          NA      NA   NA   NA    NA
    ## 2   100  150  200      NA          NA      NA   NA   NA    NA
    ## 3   200  250  300      NA          NA      NA   NA   NA    NA
    ## 4   300  350  400      NA          NA      NA   NA   NA    NA
    ## 5   400  450  500      NA          NA      NA   NA   NA    NA
    ## 6   500  550  600      NA          NA      NA   NA   NA    NA
    ## 7   600  650  700      NA          NA      NA   NA   NA    NA
    ## 8   700  750  800      NA          NA      NA   NA   NA    NA
    ## 9   800  850  900      NA          NA      NA   NA   NA    NA
    ## 10  900  950 1000      NA          NA      NA   NA   NA    NA
    ## 11 1000 1050 1100      NA          NA      NA   NA   NA    NA
    ## 12 1100 1150 1200      NA          NA      NA   NA   NA    NA
    ## 13 1200 1250 1300      NA          NA      NA   NA   NA    NA
    ## 14 1300 1350 1400      NA          NA      NA   NA   NA    NA
    ## 15 1400 1450 1500      NA          NA      NA   NA   NA    NA
    ## 16 1500 1550 1600      NA          NA      NA   NA   NA    NA
    ## 17 1600 1650 1700      NA          NA      NA   NA   NA    NA
    ## 18 1700 1750 1800      NA          NA      NA   NA   NA    NA
    ## 19 1800 1850 1900      NA          NA      NA   NA   NA    NA
    ## 20 1900 1950 2000      NA          NA      NA   NA   NA    NA
    ## 21 2000 2050 2100      NA          NA      NA   NA   NA    NA
    ## 22 2100 2150 2200      NA          NA      NA   NA   NA    NA
    ## 23 2200 2250 2300      NA          NA      NA   NA   NA    NA
    ## 24 2300 2350 2400      NA          NA      NA   NA   NA    NA
    ## 25 2400 2450 2500      NA          NA      NA   NA   NA    NA
    ## 26 2500 2550 2600      NA          NA      NA   NA   NA    NA
    ## 27 2600 2650 2700      NA          NA      NA   NA   NA    NA
    ## 28 2700 2750 2800      NA          NA      NA   NA   NA    NA

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

    ##          H        CR   LAI_live LAI_expanded LAI_dead Age ObsID
    ## T1_148 800 0.6605196 0.84874773   0.84874773        0  NA  <NA>
    ## T2_168 660 0.6055642 0.70557382   0.70557382        0  NA  <NA>
    ## S1_165  80 0.8032817 0.03062604   0.03062604        0  NA  <NA>

Belowground parameters can be seen in `below`:

``` r
x$below
```

    ##        Z50  Z95 Z100
    ## T1_148 100  600   NA
    ## T2_168 300 1000   NA
    ## S1_165 200 1000   NA

and in `belowLayers`:

``` r
x$belowLayers
```

    ## $V
    ##                1         2          3           4
    ## T1_148 0.8604899 0.1194556 0.01511005 0.004944476
    ## T2_168 0.5008953 0.4505941 0.04064831 0.007862284
    ## S1_165 0.6799879 0.2737911 0.03567632 0.010544678
    ## 
    ## $L
    ##               1        2        3        4
    ## T1_148 2289.062 1566.552 2250.052 4226.166
    ## T2_168 1817.571 2100.346 2410.127 4285.194
    ## S1_165 1085.030 1380.808 2170.587 4146.637
    ## 
    ## $VGrhizo_kmax
    ##               1        2         3         4
    ## T1_148 11478618  1593494  201562.5  65957.49
    ## T2_168 36121247 32493859 2931286.8 566975.78
    ## S1_165 10941459  4405482  574055.7 169670.90
    ## 
    ## $VCroot_kmax
    ##               1         2          3           4
    ## T1_148 2.382795 0.4833484 0.04256689 0.007416044
    ## T2_168 1.568929 1.2213562 0.09601747 0.010445417
    ## S1_165 2.407779 0.7618041 0.06314806 0.009770000
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

The `spwbInput`object also includes cohort parameter values for several
kinds of traits. For example, plant anatomy parameters are described in
`paramsAnatomy`:

``` r
x$paramsAnatomy
```

    ##        Hmed    Al2As      SLA LeafWidth LeafDensity WoodDensity FineRootDensity
    ## T1_148  850 1317.523 5.140523 0.1384772   0.2982842   0.6077016       0.2982842
    ## T2_168  500 3908.823 6.340000 1.7674359   0.4893392   0.9008264       0.4893392
    ## S1_165   80 4189.325 4.980084 1.3761085   0.3709679   0.4389106       0.3709679
    ##        conduit2sapwood      SRL RLD     r635
    ## T1_148       0.9236406 3172.572  10 1.964226
    ## T2_168       0.6238125 4398.812  10 1.805872
    ## S1_165       0.6238125 4398.812  10 2.289452

Parameters related to plant transpiration and photosynthesis can be seen
in `paramsTranspiration`:

``` r
x$paramsTranspiration
```

    ##             Gswmin    Gswmax  Vmax298  Jmax298 Kmax_stemxylem Kmax_rootxylem
    ## T1_148 0.003086667 0.2850000 72.19617 124.1687           0.15           0.60
    ## T2_168 0.004473333 0.2007222 68.51600 118.7863           0.40           1.60
    ## S1_165 0.010455247 0.2830167 62.78100 118.4486           0.29           1.16
    ##        VCleaf_kmax VCleafapo_kmax VCleaf_slope VCleaf_P50  VCleaf_c  VCleaf_d
    ## T1_148    4.000000        8.00000    133.86620  -2.303772 11.137050 -2.380849
    ## T2_168    4.000000        8.00000     19.14428  -1.964085  1.339370 -2.582279
    ## S1_165    9.579077       19.15815     25.47382  -2.663333  2.254991 -3.133381
    ##        kleaf_symp VCstem_kmax VCstem_slope VCstem_P50  VCstem_c  VCstem_d
    ## T1_148    8.00000    1.339563     68.30291  -5.139633 12.709999 -5.290000
    ## T2_168    8.00000    1.620936     14.60786  -6.964747  3.560000 -7.720000
    ## S1_165   19.15815    4.599269     12.31134  -6.980000  3.095442 -7.857378
    ##        VCroot_kmax VCroot_slope VCroot_P50  VCroot_c  VCroot_d VGrhizo_kmax
    ## T1_148    2.916127    103.96607  -2.966325 11.137050 -3.065569     13339632
    ## T2_168    2.896748     22.32794  -1.684034  1.339370 -2.214081     72113368
    ## S1_165    3.242501     31.37100  -1.173000  1.402489 -1.523324     16090667
    ##        Plant_kmax   FR_leaf   FR_stem   FR_root
    ## T1_148  0.7465846 0.1866462 0.5573346 0.2560193
    ## T2_168  0.8249857 0.2062464 0.5089563 0.2847972
    ## S1_165  1.5867376 0.1656462 0.3449978 0.4893561

Parameters related to pressure-volume curves and water storage capacity
of leaf and stem organs are in `paramsWaterStorage`:

``` r
x$paramsWaterStorage
```

    ##           maxFMC maxMCleaf maxMCstem   LeafPI0   LeafEPS LeafAF     Vleaf
    ## T1_148 126.03063  151.9063  99.19498 -1.591429  8.918571 0.3525 0.5258525
    ## T2_168  93.15304  131.4346  45.64970 -1.483333 19.260000 0.1700 0.2199087
    ## S1_165  96.53441   47.3984 134.64052 -2.370000 17.230000 0.2400 0.4108968
    ##          StemPI0   StemEPS    StemAF Vsapwood
    ## T1_148 -2.008039 13.256355 0.9236406 6.174277
    ## T2_168 -3.227438 46.420610 0.6238125 1.278142
    ## S1_165 -1.305868  6.297155 0.6238125 1.064511

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

### Stomatal regulation and photosynthesis

The soil water balance model determines stomatal conductance and
transpiration separately for sunlit and shade leaves. Stomatal
conductance is determined after building a photosynthesis function
corresponding to the supply function and finding the value of stomatal
conductance that maximizes carbon revenue while avoiding hydraulic
damage (a profit-maximization approach). Given a meteorological and soil
inputs and a chosen day and timestep, function
[`transp_stomatalRegulationPlot()`](https://emf-creaf.github.io/medfate/reference/transp_stomatalregulation.md)
allows displaying the supply and photosynthesis curves for sunlit and
shade leaves, along with an indication of the values corresponding to
the chosen stomatal aperture:

``` r
d <- 100
transp_stomatalRegulationPlot(x, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="E")
```

    ## Package 'meteoland' [ver. 2.2.4]

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-19-1.png)

``` r
transp_stomatalRegulationPlot(x, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="An")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-19-2.png)

``` r
transp_stomatalRegulationPlot(x, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="Gsw")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-19-3.png)

``` r
transp_stomatalRegulationPlot(x, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="T")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-19-4.png)

``` r
transp_stomatalRegulationPlot(x, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="VPD")
```

![](AdvancedWaterEnergyBalance_files/figure-html/unnamed-chunk-19-5.png)

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
    ##              3.90233421              0.00000000              0.00000000 
    ##                 NetRain                Snowmelt                   Runon 
    ##              0.00000000              0.00000000              0.00000000 
    ##            Infiltration      InfiltrationExcess        SaturationExcess 
    ##              0.00000000              0.00000000              0.00000000 
    ##                  Runoff            DeepDrainage         CapillarityRise 
    ##              0.00000000              0.00000000              0.00000000 
    ##         SoilEvaporation       HerbTranspiration         PlantExtraction 
    ##              0.50000000              0.04872542              0.68347525 
    ##           Transpiration HydraulicRedistribution 
    ##              0.68347525              0.00000000

And `Soil` contains water evaporated from each soil layer, water
transpired from each soil layer and the final soil water potential:

``` r
sd1$Soil
```

    ##           Psi HerbTranspiration HydraulicInput HydraulicOutput PlantExtraction
    ## 1 -0.03546870      0.0444001775              0     0.493843783     0.493843783
    ## 2 -0.03323777      0.0034620610              0     0.173313552     0.173313552
    ## 3 -0.03303084      0.0006078123              0     0.014326790     0.014326790
    ## 4 -0.03301159      0.0002553696              0     0.001991124     0.001991124

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
    ## T1_148 0.84874773 0.84874773 92.18285 0.43983754    0.43983754
    ## T2_168 0.70557382 0.70557382 72.36365 0.23125117    0.23125117
    ## S1_165 0.03062604 0.03062604 44.32407 0.01238653    0.01238653
    ##        GrossPhotosynthesis NetPhotosynthesis    RootPsi    StemPsi LeafPLC
    ## T1_148           2.1295147        2.01694450 -0.4190276 -1.2593181       0
    ## T2_168           1.7448236        1.63767365 -0.3353606 -0.8559196       0
    ## S1_165           0.0666783        0.06278927 -0.3759137 -0.6029081       0
    ##             StemPLC LeafPsiMin  LeafPsiMax      dEdP        DDS   StemRWC
    ## T1_148 5.822231e-09 -1.5416393 -0.03994108 0.5034372 0.00162010 0.9971770
    ## T2_168 2.174009e-04 -1.1346709 -0.04019180 0.5264951 0.05511836 0.9973200
    ## S1_165 1.941646e-04 -0.7151737 -0.04188680 1.0418297 0.02787710 0.9880834
    ##          LeafRWC      LFMC WaterBalance
    ## T1_148 0.9581789 122.65887 4.895173e-17
    ## T2_168 0.9827563  91.84342 7.128629e-18
    ## S1_165 0.9891112  95.40533 1.162129e-18

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
    ## T1_148  -1.642069 -0.03994108 0.002242631 0.27455938 1.274602 12.14597
    ## T2_168  -1.593485 -0.04199022 0.003269028 0.07398477 1.272483 17.43958
    ## S1_165  -1.122060 -0.04703362 0.007561239 0.09776418 1.267959 17.16641

``` r
sd1$ShadeLeaves
```

    ##        LeafPsiMin  LeafPsiMax      GSWMin     GSWMax   TempMin  TempMax
    ## T1_148 -1.4414211 -0.03994108 0.002249194 0.27901525 0.9445492 10.29826
    ## T2_168 -0.7791464 -0.04019180 0.003253667 0.08234382 0.5424548 10.38280
    ## S1_165 -0.5152537 -0.04188680 0.007596148 0.10706902 0.6786286 10.16417

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

    ## [1] 0.9847706 0.9984846 0.9998028 0.9999259

And the temperature of soil layers:

``` r
x$soil$Temp
```

    ## [1] 8.088756 3.159986 2.383625 2.363290

We can also see the current state of canopy variables:

``` r
x$canopy
```

    ##    zlow zmid  zup    LAIlive LAIexpanded LAIdead     Tair Cair     VPair
    ## 1     0   50  100 0.03062604  0.03062604       0 5.684231  386 0.5170718
    ## 2   100  150  200 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 3   200  250  300 0.06200693  0.06200693       0 5.684231  386 0.5170718
    ## 4   300  350  400 0.29933459  0.29933459       0 5.684231  386 0.5170718
    ## 5   400  450  500 0.43266592  0.43266592       0 5.684231  386 0.5170718
    ## 6   500  550  600 0.41005056  0.41005056       0 5.684231  386 0.5170718
    ## 7   600  650  700 0.24368578  0.24368578       0 5.684231  386 0.5170718
    ## 8   700  750  800 0.10657777  0.10657777       0 5.684231  386 0.5170718
    ## 9   800  850  900 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 10  900  950 1000 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 11 1000 1050 1100 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 12 1100 1150 1200 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 13 1200 1250 1300 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 14 1300 1350 1400 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 15 1400 1450 1500 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 16 1500 1550 1600 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 17 1600 1650 1700 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 18 1700 1750 1800 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 19 1800 1850 1900 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 20 1900 1950 2000 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 21 2000 2050 2100 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 22 2100 2150 2200 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 23 2200 2250 2300 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 24 2300 2350 2400 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 25 2400 2450 2500 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 26 2500 2550 2600 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 27 2600 2650 2700 0.00000000  0.00000000       0 5.684231  386 0.5170718
    ## 28 2700 2750 2800 0.00000000  0.00000000       0 5.684231  386 0.5170718

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

    ##    zlow zmid  zup    LAIlive LAIexpanded LAIdead Tair Cair VPair
    ## 1     0   50  100 0.03062604  0.03062604       0   NA   NA    NA
    ## 2   100  150  200 0.00000000  0.00000000       0   NA   NA    NA
    ## 3   200  250  300 0.06200693  0.06200693       0   NA   NA    NA
    ## 4   300  350  400 0.29933459  0.29933459       0   NA   NA    NA
    ## 5   400  450  500 0.43266592  0.43266592       0   NA   NA    NA
    ## 6   500  550  600 0.41005056  0.41005056       0   NA   NA    NA
    ## 7   600  650  700 0.24368578  0.24368578       0   NA   NA    NA
    ## 8   700  750  800 0.10657777  0.10657777       0   NA   NA    NA
    ## 9   800  850  900 0.00000000  0.00000000       0   NA   NA    NA
    ## 10  900  950 1000 0.00000000  0.00000000       0   NA   NA    NA
    ## 11 1000 1050 1100 0.00000000  0.00000000       0   NA   NA    NA
    ## 12 1100 1150 1200 0.00000000  0.00000000       0   NA   NA    NA
    ## 13 1200 1250 1300 0.00000000  0.00000000       0   NA   NA    NA
    ## 14 1300 1350 1400 0.00000000  0.00000000       0   NA   NA    NA
    ## 15 1400 1450 1500 0.00000000  0.00000000       0   NA   NA    NA
    ## 16 1500 1550 1600 0.00000000  0.00000000       0   NA   NA    NA
    ## 17 1600 1650 1700 0.00000000  0.00000000       0   NA   NA    NA
    ## 18 1700 1750 1800 0.00000000  0.00000000       0   NA   NA    NA
    ## 19 1800 1850 1900 0.00000000  0.00000000       0   NA   NA    NA
    ## 20 1900 1950 2000 0.00000000  0.00000000       0   NA   NA    NA
    ## 21 2000 2050 2100 0.00000000  0.00000000       0   NA   NA    NA
    ## 22 2100 2150 2200 0.00000000  0.00000000       0   NA   NA    NA
    ## 23 2200 2250 2300 0.00000000  0.00000000       0   NA   NA    NA
    ## 24 2300 2350 2400 0.00000000  0.00000000       0   NA   NA    NA
    ## 25 2400 2450 2500 0.00000000  0.00000000       0   NA   NA    NA
    ## 26 2500 2550 2600 0.00000000  0.00000000       0   NA   NA    NA
    ## 27 2600 2650 2700 0.00000000  0.00000000       0   NA   NA    NA
    ## 28 2700 2750 2800 0.00000000  0.00000000       0   NA   NA    NA

Now we are ready to call function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md):

``` r
S <- spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)
```

    ## Initial plant water content (mm): 6.78662
    ## Initial soil water content (mm): 290.875
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2001]:............
    ## 
    ## Final plant water content (mm): 6.78165
    ## Final soil water content (mm): 273.985
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.00497469
    ## Plant water balance result (mm): 2.04696e-16
    ## Change in soil water content (mm): -16.8895
    ## Soil water balance result (mm): -16.8895
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): -7.10543e-15
    ## Water balance components:
    ##   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
    ##   Interception (mm) 92 Net rainfall (mm) 370
    ##   Infiltration (mm) 402 Infiltration excess (mm) 19 Saturation excess (mm) 0 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 22  Herbaceous transpiration (mm) 13 Woody plant transpiration (mm) 242
    ##   Plant extraction from soil (mm) 242  Plant water balance (mm) 0 Hydraulic redistribution (mm) 4
    ##   Runoff (mm) 19 Deep drainage (mm) 142

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
    ## 2001-01-01 0.8828475      4.869109  4.869109    0  3.4241795        0
    ## 2001-01-02 1.6375337      2.498292  2.498292    0  1.0717469        0
    ## 2001-01-03 1.3017026      0.000000  0.000000    0  0.0000000        0
    ## 2001-01-04 0.5690790      5.796973  5.796973    0  4.3625616        0
    ## 2001-01-05 1.6760567      1.884401  1.884401    0  0.7539027        0
    ## 2001-01-06 1.2077028     13.359801 13.359801    0 11.7240275        0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2001-01-01    3.4241795                  0                0      0    2.9587033
    ## 2001-01-02    1.0717469                  0                0      0    0.4987713
    ## 2001-01-03    0.0000000                  0                0      0    0.0000000
    ## 2001-01-04    4.3625616                  0                0      0    3.3039457
    ## 2001-01-05    0.7539027                  0                0      0    0.1775195
    ## 2001-01-06   11.7240275                  0                0      0    4.1214138
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2001-01-01               0          1.9104055     1.444929       0.4478948
    ## 2001-01-02               0          1.9995203     1.426545       0.5000000
    ## 2001-01-03               0          0.8999863     0.000000       0.5000000
    ## 2001-01-04               0          1.5930409     1.434411       0.1435328
    ## 2001-01-05               0          1.7579878     1.130499       0.5000000
    ## 2001-01-06               0          2.1071751     1.635773       0.4514992
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01       0.011023432     0.006557998   0.006557998
    ## 2001-01-02       0.020446613     0.052528994   0.052528994
    ## 2001-01-03       0.016253351     0.383732929   0.383732929
    ## 2001-01-04       0.007105368     0.007991412   0.007991412
    ## 2001-01-05       0.020927620     0.106561357   0.106561357
    ## 2001-01-06       0.015079613     0.004823168   0.004823168
    ##            HydraulicRedistribution
    ## 2001-01-01            0.000000e+00
    ## 2001-01-02            0.000000e+00
    ## 2001-01-03            0.000000e+00
    ## 2001-01-04            0.000000e+00
    ## 2001-01-05            0.000000e+00
    ## 2001-01-06            1.318394e-05

Elements `Plants` is itself a list with several elements that contain
daily output results by plant cohorts, for example leaf minimum (midday)
water potentials are:

``` r
head(S$Plants$LeafPsiMin)
```

    ##               T1_148     T2_168     S1_165
    ## 2001-01-01 -1.256406 -0.7425001 -0.3902716
    ## 2001-01-02 -1.368721 -0.6525553 -0.3961108
    ## 2001-01-03 -1.331580 -0.8165961 -0.4527407
    ## 2001-01-04 -1.223824 -0.6063900 -0.3488515
    ## 2001-01-05 -1.349259 -0.6851557 -0.4326192
    ## 2001-01-06 -1.299799 -0.6350143 -0.3832306

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
    ## 2001-01-01 1.011397    2.41127383 1.87415609 0.5371177 1.34613589 0.42235503
    ## 2001-02-01 2.278646    0.17855109 0.08778069 0.0907704 0.03511889 0.19831578
    ## 2001-03-01 2.368035    2.41917349 2.41917349 0.0000000 1.93933438 0.01762496
    ## 2001-04-01 3.086567    0.63056064 0.29195973 0.3386009 0.13472156 0.33860091
    ## 2001-05-01 3.662604    0.76337345 0.76337345 0.0000000 0.57881566 0.00000000
    ## 2001-06-01 5.265359    0.21959509 0.21959509 0.0000000 0.15746695 0.00000000
    ## 2001-07-01 4.443053    3.27810591 3.27810591 0.0000000 2.81775241 0.00000000
    ## 2001-08-01 4.463242    1.92222891 1.92222891 0.0000000 1.55188110 0.00000000
    ## 2001-09-01 3.453891    1.30651303 1.30651303 0.0000000 1.04946567 0.00000000
    ## 2001-10-01 2.405506    1.33598175 1.33598175 0.0000000 1.05395924 0.00000000
    ## 2001-11-01 1.716591    2.20566281 1.47764599 0.7280168 1.33350571 0.72801682
    ## 2001-12-01 1.608082    0.05046181 0.05046181 0.0000000 0.02018853 0.00000000
    ##            Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2001-01-01   1.76849092         0.00000000                0 0.00000000
    ## 2001-02-01   0.23343467         0.00000000                0 0.00000000
    ## 2001-03-01   1.95695933         0.00000000                0 0.00000000
    ## 2001-04-01   0.47332247         0.00000000                0 0.00000000
    ## 2001-05-01   0.57881566         0.00000000                0 0.00000000
    ## 2001-06-01   0.15746695         0.00000000                0 0.00000000
    ## 2001-07-01   2.57550745         0.24224496                0 0.24224496
    ## 2001-08-01   1.51548112         0.03639998                0 0.03639998
    ## 2001-09-01   1.04946567         0.00000000                0 0.00000000
    ## 2001-10-01   0.95308781         0.10087142                0 0.10087142
    ## 2001-11-01   1.82721466         0.23430787                0 0.23430787
    ## 2001-12-01   0.02018853         0.00000000                0 0.00000000
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2001-01-01   1.44758966               0          0.9340615   0.52802019
    ## 2001-02-01   0.01127971               0          0.5698550   0.05266179
    ## 2001-03-01   1.27019065               0          1.0903517   0.47983911
    ## 2001-04-01   0.00000000               0          0.8124561   0.15723817
    ## 2001-05-01   0.06486495               0          1.1244635   0.18455779
    ## 2001-06-01   0.00000000               0          1.2236277   0.06212814
    ## 2001-07-01   0.00000000               0          1.6049742   0.46035350
    ## 2001-08-01   0.10052299               0          1.5957976   0.37034782
    ## 2001-09-01   0.05147845               0          1.0891619   0.25704736
    ## 2001-10-01   0.44595880               0          0.9896279   0.28202252
    ## 2001-11-01   1.21970539               0          0.6280551   0.14414028
    ## 2001-12-01   0.00000000               0          0.4463346   0.03027328
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2001-01-01     0.170070965        0.01262819       0.2233421     0.2233421
    ## 2001-02-01     0.039758590        0.02844309       0.4489915     0.4489915
    ## 2001-03-01     0.120166740        0.02955994       0.4607859     0.4607859
    ## 2001-04-01     0.012059561        0.03849435       0.6046640     0.6046640
    ## 2001-05-01     0.073971446        0.04566666       0.8202676     0.8202676
    ## 2001-06-01     0.004705277        0.06308762       1.0937066     1.0937066
    ## 2001-07-01     0.090027076        0.05492135       0.9996723     0.9996723
    ## 2001-08-01     0.043429912        0.05567784       1.1263420     1.1263420
    ## 2001-09-01     0.044398488        0.04310265       0.7446134     0.7446134
    ## 2001-10-01     0.052609415        0.03002612       0.6249698     0.6249698
    ## 2001-11-01     0.058375882        0.02142490       0.4041141     0.4041141
    ## 2001-12-01     0.016188020        0.02006189       0.3798114     0.3798114
    ##            HydraulicRedistribution
    ## 2001-01-01            6.918653e-07
    ## 2001-02-01            2.635511e-04
    ## 2001-03-01            9.990381e-04
    ## 2001-04-01            5.625717e-03
    ## 2001-05-01            4.493326e-03
    ## 2001-06-01            9.805231e-02
    ## 2001-07-01            2.316109e-02
    ## 2001-08-01            2.129799e-03
    ## 2001-09-01            1.742284e-03
    ## 2001-10-01            5.064231e-04
    ## 2001-11-01            1.579466e-03
    ## 2001-12-01            2.412461e-03

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of plant cohorts by months:

``` r
summary(S, freq="months",FUN=mean, output="PlantStress")
```

    ##                  T1_148     T2_168     S1_165
    ## 2001-01-01 0.0002425468 0.02014629 0.01093263
    ## 2001-02-01 0.0007935805 0.03565554 0.01957627
    ## 2001-03-01 0.0010323840 0.04068858 0.02225634
    ## 2001-04-01 0.0013955918 0.04998857 0.03007254
    ## 2001-05-01 0.0029743016 0.06859507 0.04096868
    ## 2001-06-01 0.0209510693 0.12143927 0.12707895
    ## 2001-07-01 0.0083441325 0.09544334 0.06788208
    ## 2001-08-01 0.0063815217 0.09315796 0.05681723
    ## 2001-09-01 0.0028046184 0.06572578 0.03782227
    ## 2001-10-01 0.0019498386 0.05308549 0.02940347
    ## 2001-11-01 0.0006230769 0.03228246 0.01844611
    ## 2001-12-01 0.0006651983 0.03284306 0.01996970

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
    ## 2001-01-01     1.239166e-04       0.009799371   0.02111355
    ## 2001-01-02     2.690256e-04       0.009886504   0.01865970
    ## 2001-01-03     2.951949e-04       0.012992165   0.02712720
    ## 2001-01-04     8.704417e-05       0.009361148   0.01819566
    ## 2001-01-05     2.837136e-04       0.010119591   0.01877508
    ## 2001-01-06     1.669594e-04       0.008493509   0.01632105

Or we can combine the aggregation by species with a temporal aggregation
(here monthly averages):

``` r
summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Pinus halepensis Quercus coccifera Quercus ilex
    ## 2001-01-01     0.0002425468        0.01093263   0.02014629
    ## 2001-02-01     0.0007935805        0.01957627   0.03565554
    ## 2001-03-01     0.0010323840        0.02225634   0.04068858
    ## 2001-04-01     0.0013955918        0.03007254   0.04998857
    ## 2001-05-01     0.0029743016        0.04096868   0.06859507
    ## 2001-06-01     0.0209510693        0.12707895   0.12143927
    ## 2001-07-01     0.0083441325        0.06788208   0.09544334
    ## 2001-08-01     0.0063815217        0.05681723   0.09315796
    ## 2001-09-01     0.0028046184        0.03782227   0.06572578
    ## 2001-10-01     0.0019498386        0.02940347   0.05308549
    ## 2001-11-01     0.0006230769        0.01844611   0.03228246
    ## 2001-12-01     0.0006651983        0.01996970   0.03284306

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
