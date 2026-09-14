# Forest dynamics

## About this vignette

This document describes how to run the forest dynamics model of
`medfate`, described in De Cáceres et al. (2023) and implemented in
function
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md).
This document is meant to teach users to run the simulation model with
function
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md).
Details of the model design and formulation can be found at the
corresponding chapters of the [medfate
book](https://emf-creaf.github.io/medfatebook/index.html).

Because the model builds on the growth and water balance models, the
reader is assumed here to be familiarized with
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md) and
[`growth()`](https://emf-creaf.github.io/medfate/reference/growth.md)
(otherwise read vignettes [*Basic water
balance*](https://emf-creaf.github.io/medfate/articles/runmodels/BasicWaterBalance.html)
and [*Forest
growth*](https://emf-creaf.github.io/medfate/articles/runmodels/ForestGrowth.html)).

## Preparing model inputs

Any forest dynamics model needs information on climate, vegetation and
soils of the forest stand to be simulated. Moreover, since models in
`medfate` differentiate between species, information on species-specific
model parameters is also needed. In this subsection we explain the
different steps to prepare the data needed to run function
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

Model inputs are explained in greater detail in vignettes
[*Understanding model
inputs*](https://emf-creaf.github.io/medfate/articles/intro/UnderstandingInputs.html)
and [*Preparing model
inputs*](https://emf-creaf.github.io/medfate/articles/intro/PreparingInputs.html).
Here we only review the different steps required to run function
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

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

We can keep track of cohort age if we define a column called `Age` in
tree or shrub data, for example let us assume we know the age of the two
tree cohorts:

``` r

exampleforest$treeData$Age <- c(40, 24)
```

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

Here we will run simulations of forest dynamics using the basic water
balance model (i.e. `transpirationMode = "Granier"`). The complexity of
the soil water balance calculations can be changed by using `"Sperry"`
as input to
[`defaultControl()`](https://emf-creaf.github.io/medfate/reference/defaultControl.md).
However, when running
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
sub-daily output will never be stored (i.e. setting
`subdailyResults = TRUE` is useless).

## Executing the forest dynamics model

In this vignette we will fake a ten-year weather input by repeating the
example weather data frame ten times.

``` r

meteo <- rbind(examplemeteo, examplemeteo, examplemeteo, examplemeteo,
                    examplemeteo, examplemeteo, examplemeteo, examplemeteo,
                    examplemeteo, examplemeteo)
meteo$dates <- as.character(seq(as.Date("2001-01-01"), 
                                as.Date("2010-12-29"), by="day"))
```

Now we run the forest dynamics model using all inputs (note that no
intermediate input object is needed, as in
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
[`growth()`](https://emf-creaf.github.io/medfate/reference/growth.md)):

``` r

fd<-fordyn(exampleforest, examplesoil, SpParamsMED, meteo, control, 
           latitude = 41.82592, elevation = 100)
```

    ## Simulating year 2001 (1/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2002 (2/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2003 (3/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2004 (4/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2005 (5/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2006 (6/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2007 (7/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2008 (8/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2009 (9/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1
    ## Simulating year 2010 (10/10):  (a) Growth/mortality, (b) Regeneration nT = 2 nS = 1

It is worth noting that, while
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
calls function
[`growth()`](https://emf-creaf.github.io/medfate/reference/growth.md)
internally for each simulated year, the `verbose` option of the control
parameters only affects function
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
(i.e. all console output from
[`growth()`](https://emf-creaf.github.io/medfate/reference/growth.md) is
hidden). Recruitment and summaries are done only once a year at the
level of function
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

## Inspecting model outputs

### Stand, species and cohort summaries and plots

Among other outputs, function
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
calculates standard summary statistics that describe the structural and
compositional state of the forest at each time step. For example, we can
access stand-level statistics using:

``` r

fd$StandSummary
```

    ##    Step NumTreeSpecies NumTreeCohorts NumShrubSpecies NumShrubCohorts
    ## 1     0              2              2               1               1
    ## 2     1              2              2               1               1
    ## 3     2              2              2               1               1
    ## 4     3              2              2               1               1
    ## 5     4              2              2               1               1
    ## 6     5              2              2               1               1
    ## 7     6              2              2               1               1
    ## 8     7              2              2               1               1
    ## 9     8              2              2               1               1
    ## 10    9              2              2               1               1
    ## 11   10              2              2               1               1
    ##    TreeDensityLive TreeBasalAreaLive DominantTreeHeight DominantTreeDiameter
    ## 1         552.0000          25.03330           800.0000             37.55000
    ## 2         551.3673          25.14751           802.6528             37.59931
    ## 3         550.7312          25.25992           805.2975             37.64859
    ## 4         550.0917          25.37262           807.9395             37.69792
    ## 5         549.4470          25.48540           810.5746             37.74723
    ## 6         548.8007          25.59857           813.2082             37.79662
    ## 7         548.1510          25.71189           815.8383             37.84606
    ## 8         547.4979          25.82532           818.4644             37.89553
    ## 9         546.8395          25.93872           821.0864             37.94503
    ## 10        546.1795          26.05234           823.7061             37.99460
    ## 11        545.5197          26.16623           826.3204             38.04418
    ##    QuadraticMeanTreeDiameter HartBeckingIndex ShrubCoverLive BasalAreaDead
    ## 1                   24.02949         53.20353       3.750000    0.00000000
    ## 2                   24.09806         53.05811       3.859594    0.03899665
    ## 3                   24.16580         52.91440       3.921379    0.03938876
    ## 4                   24.23373         52.77201       3.989613    0.03978145
    ## 5                   24.30177         52.63130       4.061031    0.04028792
    ## 6                   24.37000         52.49173       4.132921    0.04057763
    ## 7                   24.43835         52.35351       4.205669    0.04098032
    ## 8                   24.50680         52.21665       4.279827    0.04138608
    ## 9                   24.57533         52.08123       4.355399    0.04190964
    ## 10                  24.64397         51.94694       4.430820    0.04220704
    ## 11                  24.71271         51.81390       4.507320    0.04238723
    ##    ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1     0.000000000            0             0
    ## 2     0.005832415            0             0
    ## 3     0.005974952            0             0
    ## 4     0.006073624            0             0
    ## 5     0.006197510            0             0
    ## 6     0.006290721            0             0
    ## 7     0.006401805            0             0
    ## 8     0.006514596            0             0
    ## 9     0.006647820            0             0
    ## 10    0.006745796            0             0
    ## 11    0.006824435            0             0

Species-level analogous statistics are shown using:

``` r

fd$SpeciesSummary
```

    ##    Step           Species NumCohorts TreeDensityLive TreeBasalAreaLive
    ## 1     0  Pinus halepensis          1        168.0000         18.604547
    ## 2     0 Quercus coccifera          1              NA                NA
    ## 3     0      Quercus ilex          1        384.0000          6.428755
    ## 4     1  Pinus halepensis          1        167.6997         18.620105
    ## 5     1 Quercus coccifera          1              NA                NA
    ## 6     1      Quercus ilex          1        383.6676          6.527402
    ## 7     2  Pinus halepensis          1        167.3978         18.635324
    ## 8     2 Quercus coccifera          1              NA                NA
    ## 9     2      Quercus ilex          1        383.3334          6.624600
    ## 10    3  Pinus halepensis          1        167.0942         18.650307
    ## 11    3 Quercus coccifera          1              NA                NA
    ## 12    3      Quercus ilex          1        382.9975          6.722314
    ## 13    4  Pinus halepensis          1        166.7881         18.664877
    ## 14    4 Quercus coccifera          1              NA                NA
    ## 15    4      Quercus ilex          1        382.6589          6.820526
    ## 16    5  Pinus halepensis          1        166.4812         18.679319
    ## 17    5 Quercus coccifera          1              NA                NA
    ## 18    5      Quercus ilex          1        382.3195          6.919250
    ## 19    6  Pinus halepensis          1        166.1726         18.693502
    ## 20    6 Quercus coccifera          1              NA                NA
    ## 21    6      Quercus ilex          1        381.9784          7.018383
    ## 22    7  Pinus halepensis          1        165.8624         18.707414
    ## 23    7 Quercus coccifera          1              NA                NA
    ## 24    7      Quercus ilex          1        381.6355          7.117903
    ## 25    8  Pinus halepensis          1        165.5496         18.720950
    ## 26    8 Quercus coccifera          1              NA                NA
    ## 27    8      Quercus ilex          1        381.2899          7.217769
    ## 28    9  Pinus halepensis          1        165.2361         18.734343
    ## 29    9 Quercus coccifera          1              NA                NA
    ## 30    9      Quercus ilex          1        380.9434          7.318000
    ## 31   10  Pinus halepensis          1        164.9226         18.747629
    ## 32   10 Quercus coccifera          1              NA                NA
    ## 33   10      Quercus ilex          1        380.5971          7.418597
    ##    ShrubCoverLive BasalAreaDead ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1              NA   0.000000000             NA            0            NA
    ## 2        3.750000            NA    0.000000000           NA             0
    ## 3              NA   0.000000000             NA            0            NA
    ## 4              NA   0.033341011             NA            0            NA
    ## 5        3.859594            NA    0.005832415           NA             0
    ## 6              NA   0.005655639             NA            0            NA
    ## 7              NA   0.033613747             NA            0            NA
    ## 8        3.921379            NA    0.005974952           NA             0
    ## 9              NA   0.005775013             NA            0            NA
    ## 10             NA   0.033885722             NA            0            NA
    ## 11       3.989613            NA    0.006073624           NA             0
    ## 12             NA   0.005895723             NA            0            NA
    ## 13             NA   0.034253231             NA            0            NA
    ## 14       4.061031            NA    0.006197510           NA             0
    ## 15             NA   0.006034685             NA            0            NA
    ## 16             NA   0.034435265             NA            0            NA
    ## 17       4.132921            NA    0.006290721           NA             0
    ## 18             NA   0.006142369             NA            0            NA
    ## 19             NA   0.034712191             NA            0            NA
    ## 20       4.205669            NA    0.006401805           NA             0
    ## 21             NA   0.006268125             NA            0            NA
    ## 22             NA   0.034990587             NA            0            NA
    ## 23       4.279827            NA    0.006514596           NA             0
    ## 24             NA   0.006395493             NA            0            NA
    ## 25             NA   0.035367260             NA            0            NA
    ## 26       4.355399            NA    0.006647820           NA             0
    ## 27             NA   0.006542380             NA            0            NA
    ## 28             NA   0.035551952             NA            0            NA
    ## 29       4.430820            NA    0.006745796           NA             0
    ## 30             NA   0.006655085             NA            0            NA
    ## 31             NA   0.035637330             NA            0            NA
    ## 32       4.507320            NA    0.006824435           NA             0
    ## 33             NA   0.006749904             NA            0            NA

Package `medfate` provides a simple `plot` function for objects of class
`fordyn`. For example, we can show the interannual variation in
stand-level basal area using:

``` r

plot(fd, type = "StandBasalArea")
```

![Stand basal area over
time](ForestDynamics_files/figure-html/unnamed-chunk-11-1.png)

### Tree/shrub tables

Another useful output of
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
are tables in long format with cohort structural information (i.e. DBH,
height, density, etc) for each time step:

``` r

fd$TreeTable
```

    ##    Step Year Cohort          Species      DBH   Height        N Z50  Z95 Z100
    ## 1     0   NA T1_158 Pinus halepensis 37.55000 800.0000 168.0000 100  300   NA
    ## 2     0   NA T2_179     Quercus ilex 14.60000 660.0000 384.0000 300 1000   NA
    ## 3     1 2001 T1_158 Pinus halepensis 37.59931 802.6528 167.6997 100  300   NA
    ## 4     1 2001 T2_179     Quercus ilex 14.71796 661.8168 383.6676 300 1000   NA
    ## 5     2 2002 T1_158 Pinus halepensis 37.64859 805.2975 167.3978 100  300   NA
    ## 6     2 2002 T2_179     Quercus ilex 14.83360 663.5916 383.3334 300 1000   NA
    ## 7     3 2003 T1_158 Pinus halepensis 37.69792 807.9395 167.0942 100  300   NA
    ## 8     3 2003 T2_179     Quercus ilex 14.94915 665.3627 382.9975 300 1000   NA
    ## 9     4 2004 T1_158 Pinus halepensis 37.74723 810.5746 166.7881 100  300   NA
    ## 10    4 2004 T2_179     Quercus ilex 15.06462 667.1305 382.6589 300 1000   NA
    ## 11    5 2005 T1_158 Pinus halepensis 37.79662 813.2082 166.4812 100  300   NA
    ## 12    5 2005 T2_179     Quercus ilex 15.17998 668.8950 382.3195 300 1000   NA
    ## 13    6 2006 T1_158 Pinus halepensis 37.84606 815.8383 166.1726 100  300   NA
    ## 14    6 2006 T2_179     Quercus ilex 15.29517 670.6544 381.9784 300 1000   NA
    ## 15    7 2007 T1_158 Pinus halepensis 37.89553 818.4644 165.8624 100  300   NA
    ## 16    7 2007 T2_179     Quercus ilex 15.41014 672.4086 381.6355 300 1000   NA
    ## 17    8 2008 T1_158 Pinus halepensis 37.94503 821.0864 165.5496 100  300   NA
    ## 18    8 2008 T2_179     Quercus ilex 15.52490 674.1577 381.2899 300 1000   NA
    ## 19    9 2009 T1_158 Pinus halepensis 37.99460 823.7061 165.2361 100  300   NA
    ## 20    9 2009 T2_179     Quercus ilex 15.63943 675.9009 380.9434 300 1000   NA
    ## 21   10 2010 T1_158 Pinus halepensis 38.04418 826.3204 164.9226 100  300   NA
    ## 22   10 2010 T2_179     Quercus ilex 15.75372 677.6379 380.5971 300 1000   NA
    ##    Age ObsID
    ## 1   40  <NA>
    ## 2   24  <NA>
    ## 3   40    NA
    ## 4   24    NA
    ## 5   41    NA
    ## 6   25    NA
    ## 7   42    NA
    ## 8   26    NA
    ## 9   43    NA
    ## 10  27    NA
    ## 11  44    NA
    ## 12  28    NA
    ## 13  45    NA
    ## 14  29    NA
    ## 15  46    NA
    ## 16  30    NA
    ## 17  47    NA
    ## 18  31    NA
    ## 19  48    NA
    ## 20  32    NA
    ## 21  49    NA
    ## 22  33    NA

The same can be shown for dead trees:

``` r

fd$DeadTreeTable
```

    ##    Step Year Cohort          Species      DBH   Height         N N_starvation
    ## 1     1 2001 T1_158 Pinus halepensis 37.59931 802.6528 0.3002818            0
    ## 2     1 2001 T2_179     Quercus ilex 14.71796 661.8168 0.3324271            0
    ## 3     2 2002 T1_158 Pinus halepensis 37.64859 805.2975 0.3019463            0
    ## 4     2 2002 T2_179     Quercus ilex 14.83360 663.5916 0.3341719            0
    ## 5     3 2003 T1_158 Pinus halepensis 37.69792 807.9395 0.3035932            0
    ## 6     3 2003 T2_179     Quercus ilex 14.94915 665.3627 0.3359033            0
    ## 7     4 2004 T1_158 Pinus halepensis 37.74723 810.5746 0.3060846            0
    ## 8     4 2004 T2_179     Quercus ilex 15.06462 667.1305 0.3385701            0
    ## 9     5 2005 T1_158 Pinus halepensis 37.79662 813.2082 0.3069075            0
    ## 10    5 2005 T2_179     Quercus ilex 15.17998 668.8950 0.3393934            0
    ## 11    6 2006 T1_158 Pinus halepensis 37.84606 815.8383 0.3085680            0
    ## 12    6 2006 T2_179     Quercus ilex 15.29517 670.6544 0.3411453            0
    ## 13    7 2007 T1_158 Pinus halepensis 37.89553 818.4644 0.3102311            0
    ## 14    7 2007 T2_179     Quercus ilex 15.41014 672.4086 0.3429026            0
    ## 15    8 2008 T1_158 Pinus halepensis 37.94503 821.0864 0.3127532            0
    ## 16    8 2008 T2_179     Quercus ilex 15.52490 674.1577 0.3456114            0
    ## 17    9 2009 T1_158 Pinus halepensis 37.99460 823.7061 0.3135666            0
    ## 18    9 2009 T2_179     Quercus ilex 15.63943 675.9009 0.3464350            0
    ## 19   10 2010 T1_158 Pinus halepensis 38.04418 826.3204 0.3135010            0
    ## 20   10 2010 T2_179     Quercus ilex 15.75372 677.6379 0.3462911            0
    ##    N_dessication N_burnt N_resprouting_stumps Z50  Z95 Z100 Age ObsID
    ## 1              0       0                    0 100  300   NA  40    NA
    ## 2              0       0                    0 300 1000   NA  24    NA
    ## 3              0       0                    0 100  300   NA  40    NA
    ## 4              0       0                    0 300 1000   NA  24    NA
    ## 5              0       0                    0 100  300   NA  41    NA
    ## 6              0       0                    0 300 1000   NA  25    NA
    ## 7              0       0                    0 100  300   NA  42    NA
    ## 8              0       0                    0 300 1000   NA  26    NA
    ## 9              0       0                    0 100  300   NA  43    NA
    ## 10             0       0                    0 300 1000   NA  27    NA
    ## 11             0       0                    0 100  300   NA  44    NA
    ## 12             0       0                    0 300 1000   NA  28    NA
    ## 13             0       0                    0 100  300   NA  45    NA
    ## 14             0       0                    0 300 1000   NA  29    NA
    ## 15             0       0                    0 100  300   NA  46    NA
    ## 16             0       0                    0 300 1000   NA  30    NA
    ## 17             0       0                    0 100  300   NA  47    NA
    ## 18             0       0                    0 300 1000   NA  31    NA
    ## 19             0       0                    0 100  300   NA  48    NA
    ## 20             0       0                    0 300 1000   NA  32    NA

### Accessing the output from function growth()

Since function
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
makes internal calls to function
[`growth()`](https://emf-creaf.github.io/medfate/reference/growth.md),
it stores the result in a vector called `GrowthResults`, which we can
use to inspect intra-annual patterns of desired variables. For example,
the following shows the leaf area for individuals of the three cohorts
during the second year:

``` r

plot(fd$GrowthResults[[2]], "LeafArea", bySpecies = T)
```

![Leaf area variation over one
year](ForestDynamics_files/figure-html/unnamed-chunk-14-1.png) Instead
of examining year by year, it is possible to plot the whole series of
results by passing a `fordyn` object to the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) function:

``` r

plot(fd, "LeafArea")
```

![Leaf area variation for multiple
years](ForestDynamics_files/figure-html/unnamed-chunk-15-1.png)

We can also create interactive plots for particular steps using function
[`shinyplot()`](https://emf-creaf.github.io/medfate/reference/shinyplot.md),
e.g.:

``` r

shinyplot(fd$GrowthResults[[1]])
```

Finally, calling function
[`extract()`](https://emf-creaf.github.io/medfate/reference/extract.md)
will extract and bind outputs for all the internal calls to function
[`growth()`](https://emf-creaf.github.io/medfate/reference/growth.md):

``` r

medfate::extract(fd, "forest", addunits = TRUE) |>
  tibble::as_tibble()
```

    ## # A tibble: 3,650 × 53
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
    ## # ℹ 3,640 more rows
    ## # ℹ 45 more variables: InfiltrationExcess [L/m^2], SaturationExcess [L/m^2],
    ## #   Runoff [L/m^2], DeepDrainage [L/m^2], CapillarityRise [L/m^2],
    ## #   Evapotranspiration [L/m^2], Interception [L/m^2], SoilEvaporation [L/m^2],
    ## #   HerbTranspiration [L/m^2], PlantExtraction [L/m^2], Transpiration [L/m^2],
    ## #   MistletoeTranspiration [L/m^2], HydraulicRedistribution [L/m^2],
    ## #   LAI [m^2/m^2], LAIherb [m^2/m^2], LAIlive [m^2/m^2], …

## Forest dynamics including management

The package allows including forest management in simulations of forest
dynamics. This is done in a very flexible manner, in the sense that
[`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
allows the user to supply an arbitrary function implementing a desired
management strategy for the stand whose dynamics are to be simulated.
The package includes, however, an in-built default function called
[`defaultManagementFunction()`](https://emf-creaf.github.io/medfate/reference/defaultManagementFunction.md)
along with a flexible parameterization, a list with defaults provided by
function
[`defaultManagementArguments()`](https://emf-creaf.github.io/medfate/reference/defaultManagementFunction.md).

Here we provide an example of simulations including forest management:

``` r

# Default arguments
args <- defaultManagementArguments()
# Here one can modify defaults before calling fordyn()
#
# Simulation
fd<-fordyn(exampleforest, examplesoil, SpParamsMED, meteo, control, 
           latitude = 41.82592, elevation = 100,
           management_function = defaultManagementFunction,
           management_args = args)
```

    ## Simulating year 2001 (1/10):  (a) Growth/mortality & management [thinning], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2002 (2/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2003 (3/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2004 (4/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2005 (5/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2006 (6/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2007 (7/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2008 (8/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2009 (9/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2
    ## Simulating year 2010 (10/10):  (a) Growth/mortality & management [none], (b) Regeneration nT = 2 nS = 2

When management is included in simulations, two additional tables are
produced, corresponding to the trees and shrubs that were cut, e.g.:

``` r

fd$CutTreeTable
```

    ##   Step Year Cohort          Species      DBH   Height          N Z50  Z95 Z100
    ## 1    1 2001 T1_158 Pinus halepensis 37.59931 802.6528   9.158143 100  300   NA
    ## 2    1 2001 T2_179     Quercus ilex 14.71796 661.8168 383.667573 300 1000   NA
    ##   Age ObsID
    ## 1  40    NA
    ## 2  24    NA

Management parameters were those of an irregular model with thinning
interventions from ‘below’, indicating that smaller trees were to be cut
earlier:

``` r

args$type
```

    ## [1] "irregular"

``` r

args$thinning
```

    ## [1] "below"

Note that in this example, there is resprouting of *Quercus ilex* after
the thinning intervention, evidenced by the new cohort (T3_168)
appearing in year 2001:

``` r

fd$TreeTable
```

    ##    Step Year Cohort          Species      DBH    Height         N Z50  Z95 Z100
    ## 1     0   NA T1_158 Pinus halepensis 37.55000 800.00000  168.0000 100  300   NA
    ## 2     0   NA T2_179     Quercus ilex 14.60000 660.00000  384.0000 300 1000   NA
    ## 3     1 2001 T1_158 Pinus halepensis 37.59931 802.65283  158.5416 100  300   NA
    ## 4     1 2001 T3_179     Quercus ilex  1.00000  47.23629 3000.0000 300 1000   NA
    ## 5     2 2002 T1_158 Pinus halepensis 37.71656 808.86594  158.3707 100  300   NA
    ## 6     2 2002 T3_179     Quercus ilex  1.00000  47.23629 2998.2956 300 1000   NA
    ## 7     3 2003 T1_158 Pinus halepensis 37.80703 813.63837  158.1989 100  300   NA
    ## 8     3 2003 T3_179     Quercus ilex  1.00000  47.23629 2996.5831 300 1000   NA
    ## 9     4 2004 T1_158 Pinus halepensis 37.85620 816.22405  158.0262 100  300   NA
    ## 10    4 2004 T3_179     Quercus ilex  1.00000  47.23629 2994.8619 300 1000   NA
    ## 11    5 2005 T1_158 Pinus halepensis 37.90537 818.80343  157.8538 100  300   NA
    ## 12    5 2005 T3_179     Quercus ilex  1.00000  47.23629 2993.1435 300 1000   NA
    ## 13    6 2006 T1_158 Pinus halepensis 37.95455 821.37809  157.6813 100  300   NA
    ## 14    6 2006 T3_179     Quercus ilex  1.00000  47.23629 2991.4231 300 1000   NA
    ## 15    7 2007 T1_158 Pinus halepensis 38.00377 823.94880  157.5086 100  300   NA
    ## 16    7 2007 T3_179     Quercus ilex  1.00000  47.23629 2989.7009 300 1000   NA
    ## 17    8 2008 T1_158 Pinus halepensis 38.05290 826.50895  157.3352 100  300   NA
    ## 18    8 2008 T3_179     Quercus ilex  1.00000  47.23629 2987.9719 300 1000   NA
    ## 19    9 2009 T1_158 Pinus halepensis 38.10207 829.06562  157.1622 100  300   NA
    ## 20    9 2009 T3_179     Quercus ilex  1.00000  47.23629 2986.2459 300 1000   NA
    ## 21   10 2010 T1_158 Pinus halepensis 38.15125 831.61694  156.9900 100  300   NA
    ## 22   10 2010 T3_179     Quercus ilex  1.00000  47.23629 2984.5274 300 1000   NA
    ##    Age ObsID
    ## 1   40  <NA>
    ## 2   24  <NA>
    ## 3   40    NA
    ## 4   24  <NA>
    ## 5   41    NA
    ## 6   24    NA
    ## 7   42    NA
    ## 8   25    NA
    ## 9   43    NA
    ## 10  26    NA
    ## 11  44    NA
    ## 12  27    NA
    ## 13  45    NA
    ## 14  28    NA
    ## 15  46    NA
    ## 16  29    NA
    ## 17  47    NA
    ## 18  30    NA
    ## 19  48    NA
    ## 20  31    NA
    ## 21  49    NA
    ## 22  32    NA

## References

- De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini
  M, García-Valdés R, Nadal-Sala D, Sabaté S, Martin-StPaul N, Morin X,
  D’Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A
  trait-enabled model to simulate Mediterranean forest function and
  dynamics at regional scales. Geoscientific Model Development 16:
  3165-3201 (<https://doi.org/10.5194/gmd-16-3165-2023>).
