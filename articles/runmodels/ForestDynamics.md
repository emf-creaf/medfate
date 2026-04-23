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
    ## Simulating year 2010 (10/10):  (a) Growth/mortality, (b) Regeneration nT = 3 nS = 1

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
    ## 11   10              2              3               1               1
    ##    TreeDensityLive TreeBasalAreaLive DominantTreeHeight DominantTreeDiameter
    ## 1         552.0000          25.03330           800.0000             37.55000
    ## 2         551.3665          25.19812           806.0398             37.66254
    ## 3         550.7277          25.36587           812.1650             37.77725
    ## 4         550.0833          25.53495           818.3192             37.89311
    ## 5         549.4316          25.70414           824.4592             38.00930
    ## 6         548.7762          25.87336           830.5710             38.12557
    ## 7         548.1151          26.04233           836.6497             38.24182
    ## 8         547.4485          26.21084           842.6868             38.35788
    ## 9         546.7745          26.37863           848.6811             38.47372
    ## 10        546.0968          26.54569           854.6312             38.58931
    ## 11       3545.4173          26.94778           860.5372             38.70464
    ##    QuadraticMeanTreeDiameter HartBeckingIndex ShrubCoverLive BasalAreaDead
    ## 1                   24.02949         53.20353       3.750000    0.00000000
    ## 2                   24.12231         52.83519       3.110315    0.03914829
    ## 3                   24.21650         52.46713       3.203175    0.03978289
    ## 4                   24.31131         52.10303       3.278440    0.04043208
    ## 5                   24.40618         51.74567       3.375165    0.04120402
    ## 6                   24.50100         51.39556       3.455862    0.04175824
    ## 7                   24.59570         51.05290       3.553698    0.04243147
    ## 8                   24.69016         50.71800       3.637475    0.04311056
    ## 9                   24.78432         50.39081       3.740642    0.04391554
    ## 10                  24.87810         50.07102       3.827436    0.04448481
    ## 11                  24.97153         49.75834       3.933696    0.04492895
    ##    ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1     0.000000000            0             0
    ## 2     0.005315610            0             0
    ## 3     0.004837298            0             0
    ## 4     0.004969985            0             0
    ## 5     0.005111929            0             0
    ## 6     0.005237396            0             0
    ## 7     0.005371235            0             0
    ## 8     0.005513560            0             0
    ## 9     0.005669192            0             0
    ## 10    0.005802661            0             0
    ## 11    0.005914416            0             0

Species-level analogous statistics are shown using:

``` r
fd$SpeciesSummary
```

    ##    Step           Species NumCohorts TreeDensityLive TreeBasalAreaLive
    ## 1     0  Pinus halepensis          1        168.0000         18.604547
    ## 2     0 Quercus coccifera          1              NA                NA
    ## 3     0      Quercus ilex          1        384.0000          6.428755
    ## 4     1  Pinus halepensis          1        167.6993         18.682730
    ## 5     1 Quercus coccifera          1              NA                NA
    ## 6     1      Quercus ilex          1        383.6672          6.515389
    ## 7     2  Pinus halepensis          1        167.3960         18.762710
    ## 8     2 Quercus coccifera          1              NA                NA
    ## 9     2      Quercus ilex          1        383.3317          6.603157
    ## 10    3  Pinus halepensis          1        167.0898         18.843449
    ## 11    3 Quercus coccifera          1              NA                NA
    ## 12    3      Quercus ilex          1        382.9935          6.691504
    ## 13    4  Pinus halepensis          1        166.7801         18.924044
    ## 14    4 Quercus coccifera          1              NA                NA
    ## 15    4      Quercus ilex          1        382.6515          6.780091
    ## 16    5  Pinus halepensis          1        166.4684         19.004416
    ## 17    5 Quercus coccifera          1              NA                NA
    ## 18    5      Quercus ilex          1        382.3078          6.868948
    ## 19    6  Pinus halepensis          1        166.1540         19.084367
    ## 20    6 Quercus coccifera          1              NA                NA
    ## 21    6      Quercus ilex          1        381.9612          6.957966
    ## 22    7  Pinus halepensis          1        165.8367         19.163720
    ## 23    7 Quercus coccifera          1              NA                NA
    ## 24    7      Quercus ilex          1        381.6118          7.047123
    ## 25    8  Pinus halepensis          1        165.5158         19.242337
    ## 26    8 Quercus coccifera          1              NA                NA
    ## 27    8      Quercus ilex          1        381.2587          7.136290
    ## 28    9  Pinus halepensis          1        165.1931         19.320381
    ## 29    9 Quercus coccifera          1              NA                NA
    ## 30    9      Quercus ilex          1        380.9038          7.225309
    ## 31   10  Pinus halepensis          1        164.8693         19.397952
    ## 32   10 Quercus coccifera          1              NA                NA
    ## 33   10      Quercus ilex          2       3380.5480          7.549830
    ##    ShrubCoverLive BasalAreaDead ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1              NA   0.000000000             NA            0            NA
    ## 2        3.750000            NA    0.000000000           NA             0
    ## 3              NA   0.000000000             NA            0            NA
    ## 4              NA   0.033496910             NA            0            NA
    ## 5        3.110315            NA    0.005315610           NA             0
    ## 6              NA   0.005651379             NA            0            NA
    ## 7              NA   0.034003864             NA            0            NA
    ## 8        3.203175            NA    0.004837298           NA             0
    ## 9              NA   0.005779028             NA            0            NA
    ## 10             NA   0.034522571             NA            0            NA
    ## 11       3.278440            NA    0.004969985           NA             0
    ## 12             NA   0.005909504             NA            0            NA
    ## 13             NA   0.035145249             NA            0            NA
    ## 14       3.375165            NA    0.005111929           NA             0
    ## 15             NA   0.006058773             NA            0            NA
    ## 16             NA   0.035581340             NA            0            NA
    ## 17       3.455862            NA    0.005237396           NA             0
    ## 18             NA   0.006176899             NA            0            NA
    ## 19             NA   0.036118094             NA            0            NA
    ## 20       3.553698            NA    0.005371235           NA             0
    ## 21             NA   0.006313376             NA            0            NA
    ## 22             NA   0.036658945             NA            0            NA
    ## 23       3.637475            NA    0.005513560           NA             0
    ## 24             NA   0.006451618             NA            0            NA
    ## 25             NA   0.037305919             NA            0            NA
    ## 26       3.740642            NA    0.005669192           NA             0
    ## 27             NA   0.006609626             NA            0            NA
    ## 28             NA   0.037751942             NA            0            NA
    ## 29       3.827436            NA    0.005802661           NA             0
    ## 30             NA   0.006732865             NA            0            NA
    ## 31             NA   0.038091353             NA            0            NA
    ## 32       3.933696            NA    0.005914416           NA             0
    ## 33             NA   0.006837601             NA            0            NA

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

    ##    Step Year Cohort          Species      DBH    Height         N      Z50  Z95
    ## 1     0   NA T1_148 Pinus halepensis 37.55000 800.00000  168.0000 100.0000  300
    ## 2     0   NA T2_168     Quercus ilex 14.60000 660.00000  384.0000 300.0000 1000
    ## 3     1 2001 T1_148 Pinus halepensis 37.66254 806.03980  167.6993 100.0000  300
    ## 4     1 2001 T2_168     Quercus ilex 14.70442 663.24378  383.6672 300.0000 1000
    ## 5     2 2002 T1_148 Pinus halepensis 37.77725 812.16500  167.3960 100.0000  300
    ## 6     2 2002 T2_168     Quercus ilex 14.80961 666.50255  383.3317 300.0000 1000
    ## 7     3 2003 T1_148 Pinus halepensis 37.89311 818.31921  167.0898 100.0000  300
    ## 8     3 2003 T2_168     Quercus ilex 14.91493 669.75755  382.9935 300.0000 1000
    ## 9     4 2004 T1_148 Pinus halepensis 38.00930 824.45916  166.7801 100.0000  300
    ## 10    4 2004 T2_168     Quercus ilex 15.02004 672.99770  382.6515 300.0000 1000
    ## 11    5 2005 T1_148 Pinus halepensis 38.12557 830.57098  166.4684 100.0000  300
    ## 12    5 2005 T2_168     Quercus ilex 15.12494 676.22324  382.3078 300.0000 1000
    ## 13    6 2006 T1_148 Pinus halepensis 38.24182 836.64974  166.1540 100.0000  300
    ## 14    6 2006 T2_168     Quercus ilex 15.22953 679.43110  381.9612 300.0000 1000
    ## 15    7 2007 T1_148 Pinus halepensis 38.35788 842.68680  165.8367 100.0000  300
    ## 16    7 2007 T2_168     Quercus ilex 15.33381 682.62080  381.6118 300.0000 1000
    ## 17    8 2008 T1_148 Pinus halepensis 38.47372 848.68107  165.5158 100.0000  300
    ## 18    8 2008 T2_168     Quercus ilex 15.43766 685.78896  381.2587 300.0000 1000
    ## 19    9 2009 T1_148 Pinus halepensis 38.58931 854.63116  165.1931 100.0000  300
    ## 20    9 2009 T2_168     Quercus ilex 15.54088 688.92955  380.9038 300.0000 1000
    ## 21   10 2010 T1_148 Pinus halepensis 38.70464 860.53725  164.8693 100.0000  300
    ## 22   10 2010 T2_168     Quercus ilex 15.64350 692.04340  380.5480 300.0000 1000
    ## 23   10 2010 T3_168     Quercus ilex  1.00000  47.23629 3000.0000 439.9059 5020
    ##    Z100 Age ObsID
    ## 1    NA  40  <NA>
    ## 2    NA  24  <NA>
    ## 3    NA  40    NA
    ## 4    NA  24    NA
    ## 5    NA  41    NA
    ## 6    NA  25    NA
    ## 7    NA  42    NA
    ## 8    NA  26    NA
    ## 9    NA  43    NA
    ## 10   NA  27    NA
    ## 11   NA  44    NA
    ## 12   NA  28    NA
    ## 13   NA  45    NA
    ## 14   NA  29    NA
    ## 15   NA  46    NA
    ## 16   NA  30    NA
    ## 17   NA  47    NA
    ## 18   NA  31    NA
    ## 19   NA  48    NA
    ## 20   NA  32    NA
    ## 21   NA  49    NA
    ## 22   NA  33    NA
    ## 23   NA   5  <NA>

The same can be shown for dead trees:

``` r
fd$DeadTreeTable
```

    ##    Step Year Cohort          Species      DBH   Height         N N_starvation
    ## 1     1 2001 T1_148 Pinus halepensis 37.66254 806.0398 0.3006739            0
    ## 2     1 2001 T2_168     Quercus ilex 14.70442 663.2438 0.3327888            0
    ## 3     2 2002 T1_148 Pinus halepensis 37.77725 812.1650 0.3033735            0
    ## 4     2 2002 T2_168     Quercus ilex 14.80961 666.5026 0.3354887            0
    ## 5     3 2003 T1_148 Pinus halepensis 37.89311 818.3192 0.3061207            0
    ## 6     3 2003 T2_168     Quercus ilex 14.91493 669.7575 0.3382351            0
    ## 7     4 2004 T1_148 Pinus halepensis 38.00930 824.4592 0.3097397            0
    ## 8     4 2004 T2_168     Quercus ilex 15.02004 672.9977 0.3419421            0
    ## 9     5 2005 T1_148 Pinus halepensis 38.12557 830.5710 0.3116733            0
    ## 10    5 2005 T2_168     Quercus ilex 15.12494 676.2232 0.3437901            0
    ## 11    6 2006 T1_148 Pinus halepensis 38.24182 836.6497 0.3144545            0
    ## 12    6 2006 T2_168     Quercus ilex 15.22953 679.4311 0.3465761            0
    ## 13    7 2007 T1_148 Pinus halepensis 38.35788 842.6868 0.3172348            0
    ## 14    7 2007 T2_168     Quercus ilex 15.33381 682.6208 0.3493644            0
    ## 15    8 2008 T1_148 Pinus halepensis 38.47372 848.6811 0.3208924            0
    ## 16    8 2008 T2_168     Quercus ilex 15.43766 685.7890 0.3531215            0
    ## 17    9 2009 T1_148 Pinus halepensis 38.58931 854.6312 0.3227865            0
    ## 18    9 2009 T2_168     Quercus ilex 15.54088 688.9295 0.3549431            0
    ## 19   10 2010 T1_148 Pinus halepensis 38.70464 860.5372 0.3237504            0
    ## 20   10 2010 T2_168     Quercus ilex 15.64350 692.0434 0.3557506            0
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
extract(fd, "forest", addunits = TRUE) |>
  tibble::as_tibble()
```

    ## # A tibble: 3,650 × 51
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
    ## # ℹ 3,640 more rows
    ## # ℹ 43 more variables: InfiltrationExcess [L/m^2], SaturationExcess [L/m^2],
    ## #   Runoff [L/m^2], DeepDrainage [L/m^2], CapillarityRise [L/m^2],
    ## #   Evapotranspiration [L/m^2], Interception [L/m^2], SoilEvaporation [L/m^2],
    ## #   HerbTranspiration [L/m^2], PlantExtraction [L/m^2], Transpiration [L/m^2],
    ## #   HydraulicRedistribution [L/m^2], LAI [m^2/m^2], LAIherb [m^2/m^2],
    ## #   LAIlive [m^2/m^2], LAIexpanded [m^2/m^2], LAIdead [m^2/m^2], Cm [L/m^2], …

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
    ## 1    1 2001 T1_148 Pinus halepensis 37.66254 806.0398   9.371545 100  300   NA
    ## 2    1 2001 T2_168     Quercus ilex 14.70442 663.2438 383.667211 300 1000   NA
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

    ##    Step Year Cohort          Species       DBH    Height         N Z50  Z95
    ## 1     0   NA T1_148 Pinus halepensis 37.550000 800.00000  168.0000 100  300
    ## 2     0   NA T2_168     Quercus ilex 14.600000 660.00000  384.0000 300 1000
    ## 3     1 2001 T1_148 Pinus halepensis 37.662535 806.03980  158.3278 100  300
    ## 4     1 2001 T3_168     Quercus ilex  1.000000  47.23629 3000.0000 300 1000
    ## 5     2 2002 T1_148 Pinus halepensis 37.778885 812.18630  158.1566 100  300
    ## 6     2 2002 T3_168     Quercus ilex  1.018272  48.34510 2942.7965 300 1000
    ## 7     3 2003 T1_148 Pinus halepensis 37.898441 818.47581  157.9844 100  300
    ## 8     3 2003 T3_168     Quercus ilex  1.034429  49.32910 2893.9568 300 1000
    ## 9     4 2004 T1_148 Pinus halepensis 38.018473 824.75661  157.8106 100  300
    ## 10    4 2004 T3_168     Quercus ilex  1.050784  50.32563 2846.0925 300 1000
    ## 11    5 2005 T1_148 Pinus halepensis 38.138562 831.00666  157.6362 100  300
    ## 12    5 2005 T3_168     Quercus ilex  1.067465  51.34262 2798.8362 300 1000
    ## 13    6 2006 T1_148 Pinus halepensis 38.258585 837.21947  157.4607 100  300
    ## 14    6 2006 T3_168     Quercus ilex  1.084522  52.38306 2752.0629 300 1000
    ## 15    7 2007 T1_148 Pinus halepensis 38.378552 843.39589  157.2841 100  300
    ## 16    7 2007 T3_168     Quercus ilex  1.101961  53.44735 2705.7874 300 1000
    ## 17    8 2008 T1_148 Pinus halepensis 38.498342 849.52993  157.1059 100  300
    ## 18    8 2008 T3_168     Quercus ilex  1.119785  54.53577 2660.0203 300 1000
    ## 19    9 2009 T1_148 Pinus halepensis 38.617930 855.62066  156.9271 100  300
    ## 20    9 2009 T3_168     Quercus ilex  1.138003  55.64882 2614.7745 300 1000
    ## 21   10 2010 T1_148 Pinus halepensis 38.737298 861.66752  156.7482 100  300
    ## 22   10 2010 T3_168     Quercus ilex  1.156625  56.78717 2570.0513 300 1000
    ##    Z100 Age ObsID
    ## 1    NA  40  <NA>
    ## 2    NA  24  <NA>
    ## 3    NA  40    NA
    ## 4    NA  24  <NA>
    ## 5    NA  41    NA
    ## 6    NA  24    NA
    ## 7    NA  42    NA
    ## 8    NA  25    NA
    ## 9    NA  43    NA
    ## 10   NA  26    NA
    ## 11   NA  44    NA
    ## 12   NA  27    NA
    ## 13   NA  45    NA
    ## 14   NA  28    NA
    ## 15   NA  46    NA
    ## 16   NA  29    NA
    ## 17   NA  47    NA
    ## 18   NA  30    NA
    ## 19   NA  48    NA
    ## 20   NA  31    NA
    ## 21   NA  49    NA
    ## 22   NA  32    NA

## References

- De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini
  M, García-Valdés R, Nadal-Sala D, Sabaté S, Martin-StPaul N, Morin X,
  D’Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A
  trait-enabled model to simulate Mediterranean forest function and
  dynamics at regional scales. Geoscientific Model Development 16:
  3165-3201 (<https://doi.org/10.5194/gmd-16-3165-2023>).
