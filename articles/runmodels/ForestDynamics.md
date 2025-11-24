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

    ## Simulating year 2001 (1/10):  (a) Growth/mortality

    ## Package 'meteoland' [ver. 2.2.4]

    ## , (b) Regeneration nT = 2 nS = 1
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
    ## 2         551.3664          25.20608           806.1349             37.66427
    ## 3         550.7271          25.37991           812.2899             37.77951
    ## 4         550.0821          25.55420           818.4330             37.89512
    ## 5         549.4296          25.72846           824.5461             38.01077
    ## 6         548.7731          25.90244           830.6197             38.12627
    ## 7         548.1108          26.07567           836.6486             38.24151
    ## 8         547.4429          26.24831           842.6313             38.35647
    ## 9         546.7673          26.42020           848.5662             38.47111
    ## 10        546.0879          26.59156           854.4526             38.58539
    ## 11        545.4065          26.76259           860.2914             38.69934
    ##    QuadraticMeanTreeDiameter HartBeckingIndex ShrubCoverLive BasalAreaDead
    ## 1                   24.02949         53.20353       3.750000    0.00000000
    ## 2                   24.12613         52.82897       3.092002    0.03916800
    ## 3                   24.22322         52.45909       3.139798    0.03982580
    ## 4                   24.32050         52.09585       3.188170    0.04049396
    ## 5                   24.41777         51.74031       3.237219    0.04128410
    ## 6                   24.51483         51.39269       3.286758    0.04185535
    ## 7                   24.61153         51.05317       3.336766    0.04254437
    ## 8                   24.70793         50.72161       3.387265    0.04323866
    ## 9                   24.80401         50.39796       3.438831    0.04405921
    ## 10                  24.89979         50.08190       3.489708    0.04464393
    ## 11                  24.99534         49.77305       3.541115    0.04510400
    ##    ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1     0.000000000            0             0
    ## 2     0.005308865            0             0
    ## 3     0.004784383            0             0
    ## 4     0.004858246            0             0
    ## 5     0.004946676            0             0
    ## 6     0.005008905            0             0
    ## 7     0.005085388            0             0
    ## 8     0.005162607            0             0
    ## 9     0.005255279            0             0
    ## 10    0.005319864            0             0
    ## 11    0.005368564            0             0

Species-level analogous statistics are shown using:

``` r
fd$SpeciesSummary
```

    ##    Step           Species NumCohorts TreeDensityLive TreeBasalAreaLive
    ## 1     0  Pinus halepensis          1        168.0000         18.604547
    ## 2     0 Quercus coccifera          1              NA                NA
    ## 3     0      Quercus ilex          1        384.0000          6.428755
    ## 4     1  Pinus halepensis          1        167.6992         18.684443
    ## 5     1 Quercus coccifera          1              NA                NA
    ## 6     1      Quercus ilex          1        383.6671          6.521640
    ## 7     2  Pinus halepensis          1        167.3956         18.764917
    ## 8     2 Quercus coccifera          1              NA                NA
    ## 9     2      Quercus ilex          1        383.3314          6.614992
    ## 10    3  Pinus halepensis          1        167.0892         18.845380
    ## 11    3 Quercus coccifera          1              NA                NA
    ## 12    3      Quercus ilex          1        382.9929          6.708824
    ## 13    4  Pinus halepensis          1        166.7790         18.925382
    ## 14    4 Quercus coccifera          1              NA                NA
    ## 15    4      Quercus ilex          1        382.6505          6.803075
    ## 16    5  Pinus halepensis          1        166.4668         19.004926
    ## 17    5 Quercus coccifera          1              NA                NA
    ## 18    5      Quercus ilex          1        382.3063          6.897510
    ## 19    6  Pinus halepensis          1        166.1517         19.083803
    ## 20    6 Quercus coccifera          1              NA                NA
    ## 21    6      Quercus ilex          1        381.9591          6.991865
    ## 22    7  Pinus halepensis          1        165.8338         19.161973
    ## 23    7 Quercus coccifera          1              NA                NA
    ## 24    7      Quercus ilex          1        381.6091          7.086334
    ## 25    8  Pinus halepensis          1        165.5121         19.239286
    ## 26    8 Quercus coccifera          1              NA                NA
    ## 27    8      Quercus ilex          1        381.2552          7.180910
    ## 28    9  Pinus halepensis          1        165.1884         19.315913
    ## 29    9 Quercus coccifera          1              NA                NA
    ## 30    9      Quercus ilex          1        380.8995          7.275647
    ## 31   10  Pinus halepensis          1        164.8637         19.391970
    ## 32   10 Quercus coccifera          1              NA                NA
    ## 33   10      Quercus ilex          1        380.5428          7.370624
    ##    ShrubCoverLive BasalAreaDead ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1              NA   0.000000000             NA            0            NA
    ## 2        3.750000            NA    0.000000000           NA             0
    ## 3              NA   0.000000000             NA            0            NA
    ## 4              NA   0.033509817             NA            0            NA
    ## 5        3.092002            NA    0.005308865           NA             0
    ## 6              NA   0.005658184             NA            0            NA
    ## 7              NA   0.034032872             NA            0            NA
    ## 8        3.139798            NA    0.004784383           NA             0
    ## 9              NA   0.005792932             NA            0            NA
    ## 10             NA   0.034563773             NA            0            NA
    ## 11       3.188170            NA    0.004858246           NA             0
    ## 12             NA   0.005930183             NA            0            NA
    ## 13             NA   0.035197604             NA            0            NA
    ## 14       3.237219            NA    0.004946676           NA             0
    ## 15             NA   0.006086497             NA            0            NA
    ## 16             NA   0.035643831             NA            0            NA
    ## 17       3.286758            NA    0.005008905           NA             0
    ## 18             NA   0.006211519             NA            0            NA
    ## 19             NA   0.036189616             NA            0            NA
    ## 20       3.336766            NA    0.005085388           NA             0
    ## 21             NA   0.006354757             NA            0            NA
    ## 22             NA   0.036738863             NA            0            NA
    ## 23       3.387265            NA    0.005162607           NA             0
    ## 24             NA   0.006499793             NA            0            NA
    ## 25             NA   0.037394259             NA            0            NA
    ## 26       3.438831            NA    0.005255279           NA             0
    ## 27             NA   0.006664956             NA            0            NA
    ## 28             NA   0.037848416             NA            0            NA
    ## 29       3.489708            NA    0.005319864           NA             0
    ## 30             NA   0.006795519             NA            0            NA
    ## 31             NA   0.038196097             NA            0            NA
    ## 32       3.541115            NA    0.005368564           NA             0
    ## 33             NA   0.006907901             NA            0            NA

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
    ## 1     0   NA T1_148 Pinus halepensis 37.55000 800.0000 168.0000 100  600   NA
    ## 2     0   NA T2_168     Quercus ilex 14.60000 660.0000 384.0000 300 1000   NA
    ## 3     1 2001 T1_148 Pinus halepensis 37.66427 806.1349 167.6992 100  600   NA
    ## 4     1 2001 T2_168     Quercus ilex 14.71147 663.4696 383.6671 300 1000   NA
    ## 5     2 2002 T1_148 Pinus halepensis 37.77951 812.2899 167.3956 100  600   NA
    ## 6     2 2002 T2_168     Quercus ilex 14.82288 666.9271 383.3314 300 1000   NA
    ## 7     3 2003 T1_148 Pinus halepensis 37.89512 818.4330 167.0892 100  600   NA
    ## 8     3 2003 T2_168     Quercus ilex 14.93423 670.3732 382.9929 300 1000   NA
    ## 9     4 2004 T1_148 Pinus halepensis 38.01077 824.5461 166.7790 100  600   NA
    ## 10    4 2004 T2_168     Quercus ilex 15.04550 673.8066 382.6505 300 1000   NA
    ## 11    5 2005 T1_148 Pinus halepensis 38.12627 830.6197 166.4668 100  600   NA
    ## 12    5 2005 T2_168     Quercus ilex 15.15638 677.2185 382.3063 300 1000   NA
    ## 13    6 2006 T1_148 Pinus halepensis 38.24151 836.6486 166.1517 100  600   NA
    ## 14    6 2006 T2_168     Quercus ilex 15.26663 680.6006 381.9591 300 1000   NA
    ## 15    7 2007 T1_148 Pinus halepensis 38.35647 842.6313 165.8338 100  600   NA
    ## 16    7 2007 T2_168     Quercus ilex 15.37646 683.9602 381.6091 300 1000   NA
    ## 17    8 2008 T1_148 Pinus halepensis 38.47111 848.5662 165.5121 100  600   NA
    ## 18    8 2008 T2_168     Quercus ilex 15.48592 687.2977 381.2552 300 1000   NA
    ## 19    9 2009 T1_148 Pinus halepensis 38.58539 854.4526 165.1884 100  600   NA
    ## 20    9 2009 T2_168     Quercus ilex 15.59501 690.6143 380.8995 300 1000   NA
    ## 21   10 2010 T1_148 Pinus halepensis 38.69934 860.2914 164.8637 100  600   NA
    ## 22   10 2010 T2_168     Quercus ilex 15.70382 693.9120 380.5428 300 1000   NA
    ##    Age ObsID
    ## 1   40  <NA>
    ## 2   24  <NA>
    ## 3   40  <NA>
    ## 4   24  <NA>
    ## 5   41  <NA>
    ## 6   25  <NA>
    ## 7   42  <NA>
    ## 8   26  <NA>
    ## 9   43  <NA>
    ## 10  27  <NA>
    ## 11  44  <NA>
    ## 12  28  <NA>
    ## 13  45  <NA>
    ## 14  29  <NA>
    ## 15  46  <NA>
    ## 16  30  <NA>
    ## 17  47  <NA>
    ## 18  31  <NA>
    ## 19  48  <NA>
    ## 20  32  <NA>
    ## 21  49  <NA>
    ## 22  33  <NA>

The same can be shown for dead trees:

``` r
fd$DeadTreeTable
```

    ##    Step Year Cohort          Species      DBH   Height         N N_starvation
    ## 1     1 2001 T1_148 Pinus halepensis 37.66427 806.1349 0.3007620            0
    ## 2     1 2001 T2_168     Quercus ilex 14.71147 663.4696 0.3328701            0
    ## 3     2 2002 T1_148 Pinus halepensis 37.77951 812.2899 0.3035960            0
    ## 4     2 2002 T2_168     Quercus ilex 14.82288 666.9271 0.3356940            0
    ## 5     3 2003 T1_148 Pinus halepensis 37.89512 818.4330 0.3064535            0
    ## 6     3 2003 T2_168     Quercus ilex 14.93423 670.3732 0.3385419            0
    ## 7     4 2004 T1_148 Pinus halepensis 38.01077 824.5461 0.3101772            0
    ## 8     4 2004 T2_168     Quercus ilex 15.04550 673.8066 0.3423454            0
    ## 9     5 2005 T1_148 Pinus halepensis 38.12627 830.6197 0.3122093            0
    ## 10    5 2005 T2_168     Quercus ilex 15.15638 677.2185 0.3442841            0
    ## 11    6 2006 T1_148 Pinus halepensis 38.24151 836.6486 0.3150822            0
    ## 12    6 2006 T2_168     Quercus ilex 15.26663 680.6006 0.3471545            0
    ## 13    7 2007 T1_148 Pinus halepensis 38.35647 842.6313 0.3179497            0
    ## 14    7 2007 T2_168     Quercus ilex 15.37646 683.9602 0.3500230            0
    ## 15    8 2008 T1_148 Pinus halepensis 38.47111 848.5662 0.3216960            0
    ## 16    8 2008 T2_168     Quercus ilex 15.48592 687.2977 0.3538617            0
    ## 17    9 2009 T1_148 Pinus halepensis 38.58539 854.4526 0.3236771            0
    ## 18    9 2009 T2_168     Quercus ilex 15.59501 690.6143 0.3557634            0
    ## 19   10 2010 T1_148 Pinus halepensis 38.69934 860.2914 0.3247297            0
    ## 20   10 2010 T2_168     Quercus ilex 15.70382 693.9120 0.3566526            0
    ##    N_dessication N_burnt Z50  Z95 Z100 Age ObsID
    ## 1              0       0 100  600   NA  40  <NA>
    ## 2              0       0 300 1000   NA  24  <NA>
    ## 3              0       0 100  600   NA  40  <NA>
    ## 4              0       0 300 1000   NA  24  <NA>
    ## 5              0       0 100  600   NA  41  <NA>
    ## 6              0       0 300 1000   NA  25  <NA>
    ## 7              0       0 100  600   NA  42  <NA>
    ## 8              0       0 300 1000   NA  26  <NA>
    ## 9              0       0 100  600   NA  43  <NA>
    ## 10             0       0 300 1000   NA  27  <NA>
    ## 11             0       0 100  600   NA  44  <NA>
    ## 12             0       0 300 1000   NA  28  <NA>
    ## 13             0       0 100  600   NA  45  <NA>
    ## 14             0       0 300 1000   NA  29  <NA>
    ## 15             0       0 100  600   NA  46  <NA>
    ## 16             0       0 300 1000   NA  30  <NA>
    ## 17             0       0 100  600   NA  47  <NA>
    ## 18             0       0 300 1000   NA  31  <NA>
    ## 19             0       0 100  600   NA  48  <NA>
    ## 20             0       0 300 1000   NA  32  <NA>

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

Finally, we can create interactive plots for particular steps using
function
[`shinyplot()`](https://emf-creaf.github.io/medfate/reference/shinyplot.md),
e.g.:

``` r
shinyplot(fd$GrowthResults[[1]])
```

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
    ## 1    1 2001 T1_148 Pinus halepensis 37.66427 806.1349   9.336019 100  600   NA
    ## 2    1 2001 T2_168     Quercus ilex 14.71147 663.4696 383.667130 300 1000   NA
    ##   Age ObsID
    ## 1  40  <NA>
    ## 2  24  <NA>

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

    ##    Step Year Cohort          Species       DBH    Height         N      Z50
    ## 1     0   NA T1_148 Pinus halepensis 37.550000 800.00000  168.0000 100.0000
    ## 2     0   NA T2_168     Quercus ilex 14.600000 660.00000  384.0000 300.0000
    ## 3     1 2001 T1_148 Pinus halepensis 37.664271 806.13493  158.3632 100.0000
    ## 4     1 2001 T3_168     Quercus ilex  1.000000  47.23629 3000.0000 300.0000
    ## 5     2 2002 T1_148 Pinus halepensis 37.781094 812.32091  158.1918 100.0000
    ## 6     2 2002 T3_168     Quercus ilex  1.138025  55.69984 2614.7333 300.0000
    ## 7     3 2003 T1_148 Pinus halepensis 37.900792 818.62002  158.0190 100.0000
    ## 8     3 2003 T3_168     Quercus ilex  1.253656  62.74543 2359.1212 300.0000
    ## 9     4 2004 T1_148 Pinus halepensis 38.020670 824.89498  157.8443 100.0000
    ## 10    4 2004 T3_168     Quercus ilex  1.369321  69.79064 2147.8583 300.0000
    ## 11    5 2005 T1_148 Pinus halepensis 38.140537 831.13600  157.6687 100.0000
    ## 12    5 2005 T2_168     Quercus ilex  1.412760  74.29166 2389.7162 282.3466
    ## 13    6 2006 T1_148 Pinus halepensis 38.257950 837.21672  157.4920 100.0000
    ## 14    6 2006 T2_168     Quercus ilex  1.529022  81.37533 1910.1852 282.3466
    ## 15    7 2007 T1_148 Pinus halepensis 38.376098 843.30336  157.3139 100.0000
    ## 16    7 2007 T2_168     Quercus ilex  1.645689  88.48298 1766.5566 282.3466
    ## 17    8 2008 T1_148 Pinus halepensis 38.494595 849.37566  157.1340 100.0000
    ## 18    8 2008 T2_168     Quercus ilex  1.763580  95.67111 1641.2866 282.3466
    ## 19    9 2009 T1_148 Pinus halepensis 38.613214 855.42201  156.9531 100.0000
    ## 20    9 2009 T2_168     Quercus ilex  1.882208 102.90984 1531.5418 282.3466
    ## 21   10 2010 T1_148 Pinus halepensis 38.731818 861.43553  156.7719 100.0000
    ## 22   10 2010 T2_168     Quercus ilex  2.001386 110.18771 1434.7802 282.3466
    ##     Z95 Z100      Age ObsID
    ## 1   600   NA 40.00000  <NA>
    ## 2  1000   NA 24.00000  <NA>
    ## 3   600   NA 40.00000  <NA>
    ## 4  1000   NA 24.00000  <NA>
    ## 5   600   NA 41.00000  <NA>
    ## 6  1000   NA 24.00000  <NA>
    ## 7   600   NA 42.00000  <NA>
    ## 8  1000   NA 25.00000  <NA>
    ## 9   600   NA 43.00000  <NA>
    ## 10 1000   NA 26.00000  <NA>
    ## 11  600   NA 45.00000  <NA>
    ## 12 1000   NA 25.61679  <NA>
    ## 13  600   NA 45.00000  <NA>
    ## 14 1000   NA 25.61679  <NA>
    ## 15  600   NA 46.00000  <NA>
    ## 16 1000   NA 26.61679  <NA>
    ## 17  600   NA 47.00000  <NA>
    ## 18 1000   NA 27.61679  <NA>
    ## 19  600   NA 48.00000  <NA>
    ## 20 1000   NA 28.61679  <NA>
    ## 21  600   NA 49.00000  <NA>
    ## 22 1000   NA 29.61679  <NA>

## References

- De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini
  M, García-Valdés R, Nadal-Sala D, Sabaté S, Martin-StPaul N, Morin X,
  D’Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A
  trait-enabled model to simulate Mediterranean forest function and
  dynamics at regional scales. Geoscientific Model Development 16:
  3165-3201 (<https://doi.org/10.5194/gmd-16-3165-2023>).
