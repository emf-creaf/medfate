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
    ## 2         551.3673          25.14663           802.6468             37.59922
    ## 3         550.7312          25.26037           805.2904             37.64849
    ## 4         550.0916          25.37443           807.9294             37.69779
    ## 5         549.4469          25.48875           810.5718             37.74726
    ## 6         548.8004          25.60370           813.2023             37.79659
    ## 7         548.1505          25.71886           815.8317             37.84602
    ## 8         547.4970          25.83417           818.4594             37.89553
    ## 9         546.8383          25.94945           821.0822             37.94505
    ## 10        546.1778          26.06533           823.6969             37.99451
    ## 11        545.5175          26.18162           826.3135             38.04412
    ##    QuadraticMeanTreeDiameter HartBeckingIndex ShrubCoverLive BasalAreaDead
    ## 1                   24.02949         53.20353       3.750000    0.00000000
    ## 2                   24.09764         53.05850       3.980119    0.03899540
    ## 3                   24.16602         52.91486       4.229441    0.03938881
    ## 4                   24.23459         52.77267       4.494970    0.03978613
    ## 5                   24.30337         52.63149       4.777099    0.04029770
    ## 6                   24.37245         52.49213       5.078686    0.04059290
    ## 7                   24.44168         52.35396       5.398639    0.04100207
    ## 8                   24.51102         52.21701       5.737831    0.04141446
    ## 9                   24.58044         52.08155       6.097254    0.04194485
    ## 10                  24.65015         51.94760       6.480797    0.04224967
    ## 11                  24.72003         51.81444       6.887360    0.04243836
    ##    ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1     0.000000000            0             0
    ## 2     0.005923307            0             0
    ## 3     0.006293002            0             0
    ## 4     0.006687784            0             0
    ## 5     0.007127848            0             0
    ## 6     0.007556200            0             0
    ## 7     0.008032917            0             0
    ## 8     0.008538383            0             0
    ## 9     0.009099863            0             0
    ## 10    0.009645135            0             0
    ## 11    0.010193119            0             0

Species-level analogous statistics are shown using:

``` r

fd$SpeciesSummary
```

    ##    Step           Species NumCohorts TreeDensityLive TreeBasalAreaLive
    ## 1     0  Pinus halepensis          1        168.0000         18.604547
    ## 2     0 Quercus coccifera          1              NA                NA
    ## 3     0      Quercus ilex          1        384.0000          6.428755
    ## 4     1  Pinus halepensis          1        167.6997         18.620014
    ## 5     1 Quercus coccifera          1              NA                NA
    ## 6     1      Quercus ilex          1        383.6676          6.526616
    ## 7     2  Pinus halepensis          1        167.3978         18.635232
    ## 8     2 Quercus coccifera          1              NA                NA
    ## 9     2      Quercus ilex          1        383.3334          6.625137
    ## 10    3  Pinus halepensis          1        167.0942         18.650177
    ## 11    3 Quercus coccifera          1              NA                NA
    ## 12    3      Quercus ilex          1        382.9975          6.724253
    ## 13    4  Pinus halepensis          1        166.7880         18.664895
    ## 14    4 Quercus coccifera          1              NA                NA
    ## 15    4      Quercus ilex          1        382.6589          6.823857
    ## 16    5  Pinus halepensis          1        166.4810         18.679274
    ## 17    5 Quercus coccifera          1              NA                NA
    ## 18    5      Quercus ilex          1        382.3194          6.924428
    ## 19    6  Pinus halepensis          1        166.1723         18.693435
    ## 20    6 Quercus coccifera          1              NA                NA
    ## 21    6      Quercus ilex          1        381.9781          7.025424
    ## 22    7  Pinus halepensis          1        165.8620         18.707363
    ## 23    7 Quercus coccifera          1              NA                NA
    ## 24    7      Quercus ilex          1        381.6351          7.126809
    ## 25    8  Pinus halepensis          1        165.5490         18.720899
    ## 26    8 Quercus coccifera          1              NA                NA
    ## 27    8      Quercus ilex          1        381.2893          7.228549
    ## 28    9  Pinus halepensis          1        165.2352         18.734161
    ## 29    9 Quercus coccifera          1              NA                NA
    ## 30    9      Quercus ilex          1        380.9426          7.331167
    ## 31   10  Pinus halepensis          1        164.9214         18.747449
    ## 32   10 Quercus coccifera          1              NA                NA
    ## 33   10      Quercus ilex          1        380.5961          7.434170
    ##    ShrubCoverLive BasalAreaDead ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1              NA   0.000000000             NA            0            NA
    ## 2        3.750000            NA    0.000000000           NA             0
    ## 3              NA   0.000000000             NA            0            NA
    ## 4              NA   0.033340488             NA            0            NA
    ## 5        3.980119            NA    0.005923307           NA             0
    ## 6              NA   0.005654907             NA            0            NA
    ## 7              NA   0.033613364             NA            0            NA
    ## 8        4.229441            NA    0.006293002           NA             0
    ## 9              NA   0.005775450             NA            0            NA
    ## 10             NA   0.033888300             NA            0            NA
    ## 11       4.494970            NA    0.006687784           NA             0
    ## 12             NA   0.005897831             NA            0            NA
    ## 13             NA   0.034259197             NA            0            NA
    ## 14       4.777099            NA    0.007127848           NA             0
    ## 15             NA   0.006038503             NA            0            NA
    ## 16             NA   0.034444547             NA            0            NA
    ## 17       5.078686            NA    0.007556200           NA             0
    ## 18             NA   0.006148358             NA            0            NA
    ## 19             NA   0.034725615             NA            0            NA
    ## 20       5.398639            NA    0.008032917           NA             0
    ## 21             NA   0.006276453             NA            0            NA
    ## 22             NA   0.035008254             NA            0            NA
    ## 23       5.737831            NA    0.008538383           NA             0
    ## 24             NA   0.006406204             NA            0            NA
    ## 25             NA   0.035389283             NA            0            NA
    ## 26       6.097254            NA    0.009099863           NA             0
    ## 27             NA   0.006555566             NA            0            NA
    ## 28             NA   0.035578422             NA            0            NA
    ## 29       6.480797            NA    0.009645135           NA             0
    ## 30             NA   0.006671250             NA            0            NA
    ## 31             NA   0.035669190             NA            0            NA
    ## 32       6.887360            NA    0.010193119           NA             0
    ## 33             NA   0.006769166             NA            0            NA

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
    ## 3     1 2001 T1_158 Pinus halepensis 37.59922 802.6468 167.6997 100  300   NA
    ## 4     1 2001 T2_179     Quercus ilex 14.71708 661.7987 383.6676 300 1000   NA
    ## 5     2 2002 T1_158 Pinus halepensis 37.64849 805.2904 167.3978 100  300   NA
    ## 6     2 2002 T2_179     Quercus ilex 14.83420 663.5957 383.3334 300 1000   NA
    ## 7     3 2003 T1_158 Pinus halepensis 37.69779 807.9294 167.0942 100  300   NA
    ## 8     3 2003 T2_179     Quercus ilex 14.95130 665.3902 382.9975 300 1000   NA
    ## 9     4 2004 T1_158 Pinus halepensis 37.74726 810.5718 166.7880 100  300   NA
    ## 10    4 2004 T2_179     Quercus ilex 15.06830 667.1809 382.6589 300 1000   NA
    ## 11    5 2005 T1_158 Pinus halepensis 37.79659 813.2023 166.4810 100  300   NA
    ## 12    5 2005 T2_179     Quercus ilex 15.18566 668.9785 382.3194 300 1000   NA
    ## 13    6 2006 T1_158 Pinus halepensis 37.84602 815.8317 166.1723 100  300   NA
    ## 14    6 2006 T2_179     Quercus ilex 15.30284 670.7708 381.9781 300 1000   NA
    ## 15    7 2007 T1_158 Pinus halepensis 37.89553 818.4594 165.8620 100  300   NA
    ## 16    7 2007 T2_179     Quercus ilex 15.41979 672.5574 381.6351 300 1000   NA
    ## 17    8 2008 T1_158 Pinus halepensis 37.94505 821.0822 165.5490 100  300   NA
    ## 18    8 2008 T2_179     Quercus ilex 15.53650 674.3381 381.2893 300 1000   NA
    ## 19    9 2009 T1_158 Pinus halepensis 37.99451 823.6969 165.2352 100  300   NA
    ## 20    9 2009 T2_179     Quercus ilex 15.65351 676.1243 380.9426 300 1000   NA
    ## 21   10 2010 T1_158 Pinus halepensis 38.04412 826.3135 164.9214 100  300   NA
    ## 22   10 2010 T2_179     Quercus ilex 15.77027 677.9035 380.5961 300 1000   NA
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
    ## 1     1 2001 T1_158 Pinus halepensis 37.59922 802.6468 0.3002785            0
    ## 2     1 2001 T2_179     Quercus ilex 14.71708 661.7987 0.3324241            0
    ## 3     2 2002 T1_158 Pinus halepensis 37.64849 805.2904 0.3019443            0
    ## 4     2 2002 T2_179     Quercus ilex 14.83420 663.5957 0.3341701            0
    ## 5     3 2003 T1_158 Pinus halepensis 37.69779 807.9294 0.3036184            0
    ## 6     3 2003 T2_179     Quercus ilex 14.95130 665.3902 0.3359265            0
    ## 7     4 2004 T1_158 Pinus halepensis 37.74726 810.5718 0.3061375            0
    ## 8     4 2004 T2_179     Quercus ilex 15.06830 667.1809 0.3386189            0
    ## 9     5 2005 T1_158 Pinus halepensis 37.79659 813.2023 0.3069907            0
    ## 10    5 2005 T2_179     Quercus ilex 15.18566 668.9785 0.3394701            0
    ## 11    6 2006 T1_158 Pinus halepensis 37.84602 815.8317 0.3086879            0
    ## 12    6 2006 T2_179     Quercus ilex 15.30284 670.7708 0.3412560            0
    ## 13    7 2007 T1_158 Pinus halepensis 37.89553 818.4594 0.3103878            0
    ## 14    7 2007 T2_179     Quercus ilex 15.41979 672.5574 0.3430473            0
    ## 15    8 2008 T1_158 Pinus halepensis 37.94505 821.0822 0.3129476            0
    ## 16    8 2008 T2_179     Quercus ilex 15.53650 674.3381 0.3457910            0
    ## 17    9 2009 T1_158 Pinus halepensis 37.99451 823.6969 0.3138015            0
    ## 18    9 2009 T2_179     Quercus ilex 15.65351 676.1243 0.3466520            0
    ## 19   10 2010 T1_158 Pinus halepensis 38.04412 826.3135 0.3137821            0
    ## 20   10 2010 T2_179     Quercus ilex 15.77027 677.9035 0.3465509            0
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

    ## # A tibble: 3,650 × 54
    ##    date           PET Precipitation    Rain    Snow NetRain Snowmelt   Runon
    ##    <date>     [L/m^2]       [L/m^2] [L/m^2] [L/m^2] [L/m^2]  [L/m^2] [L/m^2]
    ##  1 2001-01-01   0.883          4.87    4.87    0      3.66      0          0
    ##  2 2001-01-02   1.64           2.50    2.50    0      1.30      0          0
    ##  3 2001-01-03   1.30           0       0       0      0         0          0
    ##  4 2001-01-04   0.569          5.80    5.80    0      4.60      0          0
    ##  5 2001-01-05   1.68           1.88    1.88    0      0.862     0          0
    ##  6 2001-01-06   1.21          13.4    13.4     0     12.0       0          0
    ##  7 2001-01-07   0.637          5.38    0       5.38   0         0          0
    ##  8 2001-01-08   0.832          0       0       0      0         0          0
    ##  9 2001-01-09   1.98           0       0       0      0         0          0
    ## 10 2001-01-10   0.829          5.12    5.12    0      3.91      5.38       0
    ## # ℹ 3,640 more rows
    ## # ℹ 46 more variables: Infiltration [L/m^2], InfiltrationExcess [L/m^2],
    ## #   SaturationExcess [L/m^2], Runoff [L/m^2], DeepDrainage [L/m^2],
    ## #   CapillarityRise [L/m^2], Evapotranspiration [L/m^2], Interception [L/m^2],
    ## #   SoilEvaporation [L/m^2], HerbTranspiration [L/m^2],
    ## #   PlantExtraction [L/m^2], Transpiration [L/m^2],
    ## #   MistletoeTranspiration [L/m^2], HydraulicRedistribution [L/m^2], …

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
    ## 1    1 2001 T1_158 Pinus halepensis 37.59922 802.6468   9.162901 100  300   NA
    ## 2    1 2001 T2_179     Quercus ilex 14.71708 661.7987 383.667576 300 1000   NA
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
    ## 3     1 2001 T1_158 Pinus halepensis 37.59922 802.64683  158.5368 100  300   NA
    ## 4     1 2001 T3_179     Quercus ilex  1.00000  47.23629 3000.0000 300 1000   NA
    ## 5     2 2002 T1_158 Pinus halepensis 37.71758 808.91651  158.3659 100  300   NA
    ## 6     2 2002 T3_179     Quercus ilex  1.00000  47.23629 2998.2956 300 1000   NA
    ## 7     3 2003 T1_158 Pinus halepensis 37.80710 813.63671  158.1941 100  300   NA
    ## 8     3 2003 T3_179     Quercus ilex  1.00000  47.23629 2996.5829 300 1000   NA
    ## 9     4 2004 T1_158 Pinus halepensis 37.85627 816.22107  158.0214 100  300   NA
    ## 10    4 2004 T3_179     Quercus ilex  1.00000  47.23629 2994.8618 300 1000   NA
    ## 11    5 2005 T1_158 Pinus halepensis 37.90545 818.80107  157.8491 100  300   NA
    ## 12    5 2005 T3_179     Quercus ilex  1.00000  47.23629 2993.1434 300 1000   NA
    ## 13    6 2006 T1_158 Pinus halepensis 37.95467 821.37705  157.6765 100  300   NA
    ## 14    6 2006 T3_179     Quercus ilex  1.00000  47.23629 2991.4231 300 1000   NA
    ## 15    7 2007 T1_158 Pinus halepensis 38.00381 823.94367  157.5038 100  300   NA
    ## 16    7 2007 T3_179     Quercus ilex  1.00000  47.23629 2989.7009 300 1000   NA
    ## 17    8 2008 T1_158 Pinus halepensis 38.05303 826.50823  157.3305 100  300   NA
    ## 18    8 2008 T3_179     Quercus ilex  1.00000  47.23629 2987.9720 300 1000   NA
    ## 19    9 2009 T1_158 Pinus halepensis 38.10222 829.06645  157.1575 100  300   NA
    ## 20    9 2009 T3_179     Quercus ilex  1.00000  47.23629 2986.2460 300 1000   NA
    ## 21   10 2010 T1_158 Pinus halepensis 38.15135 831.61555  156.9853 100  300   NA
    ## 22   10 2010 T3_179     Quercus ilex  1.00000  47.23629 2984.5275 300 1000   NA
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
