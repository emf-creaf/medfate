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

    ## Simulating year 2001 (1/10):  (a) Growth/mortality

    ## Package 'meteoland' [ver. 2.2.5]

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
    ## 2         551.3664          25.20512           806.0859             37.66335
    ## 3         550.7271          25.37888           812.2284             37.77835
    ## 4         550.0822          25.55365           818.3835             37.89419
    ## 5         549.4297          25.72873           824.5252             38.01037
    ## 6         548.7731          25.90403           830.6392             38.12664
    ## 7         548.1108          26.07891           836.7159             38.24280
    ## 8         547.4427          26.25338           842.7513             38.35879
    ## 9         546.7669          26.42725           848.7438             38.47455
    ## 10        546.0871          26.60074           854.6924             38.59006
    ## 11        545.4054          26.77396           860.5967             38.70531
    ##    QuadraticMeanTreeDiameter HartBeckingIndex ShrubCoverLive BasalAreaDead
    ## 1                   24.02949         53.20353       3.750000    0.00000000
    ## 2                   24.12566         52.83218       3.092078    0.03916503
    ## 3                   24.22273         52.46306       3.140005    0.03982124
    ## 4                   24.32023         52.09900       3.188580    0.04049045
    ## 5                   24.41789         51.74162       3.237908    0.04128341
    ## 6                   24.51559         51.39149       3.287826    0.04185936
    ## 7                   24.61306         51.04907       3.338318    0.04255468
    ## 8                   24.71032         50.71439       3.389394    0.04325634
    ## 9                   24.80732         50.38744       3.441633    0.04408521
    ## 10                  24.90410         50.06788       3.493286    0.04467893
    ## 11                  25.00067         49.75544       3.545553    0.04514843
    ##    ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1     0.000000000            0             0
    ## 2     0.005308909            0             0
    ## 3     0.004784578            0             0
    ## 4     0.004858683            0             0
    ## 5     0.004947469            0             0
    ## 6     0.005010182            0             0
    ## 7     0.005087305            0             0
    ## 8     0.005165317            0             0
    ## 9     0.005258939            0             0
    ## 10    0.005324601            0             0
    ## 11    0.005374499            0             0

Species-level analogous statistics are shown using:

``` r
fd$SpeciesSummary
```

    ##    Step           Species NumCohorts TreeDensityLive TreeBasalAreaLive
    ## 1     0  Pinus halepensis          1        168.0000         18.604547
    ## 2     0 Quercus coccifera          1              NA                NA
    ## 3     0      Quercus ilex          1        384.0000          6.428755
    ## 4     1  Pinus halepensis          1        167.6992         18.683535
    ## 5     1 Quercus coccifera          1              NA                NA
    ## 6     1      Quercus ilex          1        383.6671          6.521581
    ## 7     2  Pinus halepensis          1        167.3957         18.763774
    ## 8     2 Quercus coccifera          1              NA                NA
    ## 9     2      Quercus ilex          1        383.3315          6.615106
    ## 10    3  Pinus halepensis          1        167.0892         18.844455
    ## 11    3 Quercus coccifera          1              NA                NA
    ## 12    3      Quercus ilex          1        382.9929          6.709197
    ## 13    4  Pinus halepensis          1        166.7791         18.924994
    ## 14    4 Quercus coccifera          1              NA                NA
    ## 15    4      Quercus ilex          1        382.6506          6.803735
    ## 16    5  Pinus halepensis          1        166.4668         19.005299
    ## 17    5 Quercus coccifera          1              NA                NA
    ## 18    5      Quercus ilex          1        382.3063          6.898731
    ## 19    6  Pinus halepensis          1        166.1517         19.085089
    ## 20    6 Quercus coccifera          1              NA                NA
    ## 21    6      Quercus ilex          1        381.9591          6.993819
    ## 22    7  Pinus halepensis          1        165.8337         19.164273
    ## 23    7 Quercus coccifera          1              NA                NA
    ## 24    7      Quercus ilex          1        381.6090          7.089104
    ## 25    8  Pinus halepensis          1        165.5119         19.242701
    ## 26    8 Quercus coccifera          1              NA                NA
    ## 27    8      Quercus ilex          1        381.2550          7.184545
    ## 28    9  Pinus halepensis          1        165.1880         19.320544
    ## 29    9 Quercus coccifera          1              NA                NA
    ## 30    9      Quercus ilex          1        380.8991          7.280195
    ## 31   10  Pinus halepensis          1        164.8631         19.397890
    ## 32   10 Quercus coccifera          1              NA                NA
    ## 33   10      Quercus ilex          1        380.5423          7.376068
    ##    ShrubCoverLive BasalAreaDead ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1              NA   0.000000000             NA            0            NA
    ## 2        3.750000            NA    0.000000000           NA             0
    ## 3              NA   0.000000000             NA            0            NA
    ## 4              NA   0.033507057             NA            0            NA
    ## 5        3.092078            NA    0.005308909           NA             0
    ## 6              NA   0.005657973             NA            0            NA
    ## 7              NA   0.034028535             NA            0            NA
    ## 8        3.140005            NA    0.004784578           NA             0
    ## 9              NA   0.005792710             NA            0            NA
    ## 10             NA   0.034560204             NA            0            NA
    ## 11       3.188580            NA    0.004858683           NA             0
    ## 12             NA   0.005930246             NA            0            NA
    ## 13             NA   0.035196391             NA            0            NA
    ## 14       3.237908            NA    0.004947469           NA             0
    ## 15             NA   0.006087017             NA            0            NA
    ## 16             NA   0.035646461             NA            0            NA
    ## 17       3.287826            NA    0.005010182           NA             0
    ## 18             NA   0.006212899             NA            0            NA
    ## 19             NA   0.036197371             NA            0            NA
    ## 20       3.338318            NA    0.005087305           NA             0
    ## 21             NA   0.006357311             NA            0            NA
    ## 22             NA   0.036752630             NA            0            NA
    ## 23       3.389394            NA    0.005165317           NA             0
    ## 24             NA   0.006503714             NA            0            NA
    ## 25             NA   0.037414815             NA            0            NA
    ## 26       3.441633            NA    0.005258939           NA             0
    ## 27             NA   0.006670397             NA            0            NA
    ## 28             NA   0.037876344             NA            0            NA
    ## 29       3.493286            NA    0.005324601           NA             0
    ## 30             NA   0.006802588             NA            0            NA
    ## 31             NA   0.038231799             NA            0            NA
    ## 32       3.545553            NA    0.005374499           NA             0
    ## 33             NA   0.006916629             NA            0            NA

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
    ## 1     0   NA T1_148 Pinus halepensis 37.55000 800.0000 168.0000 100  300   NA
    ## 2     0   NA T2_168     Quercus ilex 14.60000 660.0000 384.0000 300 1000   NA
    ## 3     1 2001 T1_148 Pinus halepensis 37.66335 806.0859 167.6992 100  300   NA
    ## 4     1 2001 T2_168     Quercus ilex 14.71141 663.4675 383.6671 300 1000   NA
    ## 5     2 2002 T1_148 Pinus halepensis 37.77835 812.2284 167.3957 100  300   NA
    ## 6     2 2002 T2_168     Quercus ilex 14.82300 666.9308 383.3315 300 1000   NA
    ## 7     3 2003 T1_148 Pinus halepensis 37.89419 818.3835 167.0892 100  300   NA
    ## 8     3 2003 T2_168     Quercus ilex 14.93465 670.3856 382.9929 300 1000   NA
    ## 9     4 2004 T1_148 Pinus halepensis 38.01037 824.5252 166.7791 100  300   NA
    ## 10    4 2004 T2_168     Quercus ilex 15.04622 673.8285 382.6506 300 1000   NA
    ## 11    5 2005 T1_148 Pinus halepensis 38.12664 830.6392 166.4668 100  300   NA
    ## 12    5 2005 T2_168     Quercus ilex 15.15772 677.2589 382.3063 300 1000   NA
    ## 13    6 2006 T1_148 Pinus halepensis 38.24280 836.7159 166.1517 100  300   NA
    ## 14    6 2006 T2_168     Quercus ilex 15.26876 680.6652 381.9591 300 1000   NA
    ## 15    7 2007 T1_148 Pinus halepensis 38.35879 842.7513 165.8337 100  300   NA
    ## 16    7 2007 T2_168     Quercus ilex 15.37947 684.0511 381.6090 300 1000   NA
    ## 17    8 2008 T1_148 Pinus halepensis 38.47455 848.7438 165.5119 100  300   NA
    ## 18    8 2008 T2_168     Quercus ilex 15.48984 687.4162 381.2550 300 1000   NA
    ## 19    9 2009 T1_148 Pinus halepensis 38.59006 854.6924 165.1880 100  300   NA
    ## 20    9 2009 T2_168     Quercus ilex 15.59989 690.7614 380.8991 300 1000   NA
    ## 21   10 2010 T1_148 Pinus halepensis 38.70531 860.5967 164.8631 100  300   NA
    ## 22   10 2010 T2_168     Quercus ilex 15.70963 694.0867 380.5423 300 1000   NA
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
    ## 1     1 2001 T1_148 Pinus halepensis 37.66335 806.0859 0.3007519            0
    ## 2     1 2001 T2_168     Quercus ilex 14.71141 663.4675 0.3328608            0
    ## 3     2 2002 T1_148 Pinus halepensis 37.77835 812.2284 0.3035759            0
    ## 4     2 2002 T2_168     Quercus ilex 14.82300 666.9308 0.3356754            0
    ## 5     3 2003 T1_148 Pinus halepensis 37.89419 818.3835 0.3064370            0
    ## 6     3 2003 T2_168     Quercus ilex 14.93465 670.3856 0.3385267            0
    ## 7     4 2004 T1_148 Pinus halepensis 38.01037 824.5252 0.3101729            0
    ## 8     4 2004 T2_168     Quercus ilex 15.04622 673.8285 0.3423415            0
    ## 9     5 2005 T1_148 Pinus halepensis 38.12664 830.6392 0.3122263            0
    ## 10    5 2005 T2_168     Quercus ilex 15.15772 677.2589 0.3442996            0
    ## 11    6 2006 T1_148 Pinus halepensis 38.24280 836.7159 0.3151285            0
    ## 12    6 2006 T2_168     Quercus ilex 15.26876 680.6652 0.3471970            0
    ## 13    7 2007 T1_148 Pinus halepensis 38.35879 842.7513 0.3180305            0
    ## 14    7 2007 T2_168     Quercus ilex 15.37947 684.0511 0.3500972            0
    ## 15    8 2008 T1_148 Pinus halepensis 38.47455 848.7438 0.3218153            0
    ## 16    8 2008 T2_168     Quercus ilex 15.48984 687.4162 0.3539713            0
    ## 17    9 2009 T1_148 Pinus halepensis 38.59006 854.6924 0.3238376            0
    ## 18    9 2009 T2_168     Quercus ilex 15.59989 690.7614 0.3559108            0
    ## 19   10 2010 T1_148 Pinus halepensis 38.70531 860.5967 0.3249329            0
    ## 20   10 2010 T2_168     Quercus ilex 15.70963 694.0867 0.3568391            0
    ##    N_dessication N_burnt Z50  Z95 Z100 Age ObsID
    ## 1              0       0 100  300   NA  40  <NA>
    ## 2              0       0 300 1000   NA  24  <NA>
    ## 3              0       0 100  300   NA  40  <NA>
    ## 4              0       0 300 1000   NA  24  <NA>
    ## 5              0       0 100  300   NA  41  <NA>
    ## 6              0       0 300 1000   NA  25  <NA>
    ## 7              0       0 100  300   NA  42  <NA>
    ## 8              0       0 300 1000   NA  26  <NA>
    ## 9              0       0 100  300   NA  43  <NA>
    ## 10             0       0 300 1000   NA  27  <NA>
    ## 11             0       0 100  300   NA  44  <NA>
    ## 12             0       0 300 1000   NA  28  <NA>
    ## 13             0       0 100  300   NA  45  <NA>
    ## 14             0       0 300 1000   NA  29  <NA>
    ## 15             0       0 100  300   NA  46  <NA>
    ## 16             0       0 300 1000   NA  30  <NA>
    ## 17             0       0 100  300   NA  47  <NA>
    ## 18             0       0 300 1000   NA  31  <NA>
    ## 19             0       0 100  300   NA  48  <NA>
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

    ## # A tibble: 3,650 × 49
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
    ## # ℹ 41 more variables: InfiltrationExcess [L/m^2], SaturationExcess [L/m^2],
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
    ## 1    1 2001 T1_148 Pinus halepensis 37.66335 806.0859   9.334398 100  300   NA
    ## 2    1 2001 T2_168     Quercus ilex 14.71141 663.4675 383.667139 300 1000   NA
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

    ##    Step Year Cohort          Species       DBH    Height         N Z50  Z95
    ## 1     0   NA T1_148 Pinus halepensis 37.550000 800.00000  168.0000 100  300
    ## 2     0   NA T2_168     Quercus ilex 14.600000 660.00000  384.0000 300 1000
    ## 3     1 2001 T1_148 Pinus halepensis 37.663355 806.08587  158.3648 100  300
    ## 4     1 2001 T3_168     Quercus ilex  1.000000  47.23629 3000.0000 300 1000
    ## 5     2 2002 T1_148 Pinus halepensis 37.780051 812.26537  158.1934 100  300
    ## 6     2 2002 T3_168     Quercus ilex  1.138186  55.70973 2614.3391 300 1000
    ## 7     3 2003 T1_148 Pinus halepensis 37.899810 818.56798  158.0206 100  300
    ## 8     3 2003 T3_168     Quercus ilex  1.254005  62.76669 2358.4238 300 1000
    ## 9     4 2004 T1_148 Pinus halepensis 38.019865 824.85248  157.8460 100  300
    ## 10    4 2004 T3_168     Quercus ilex  1.369893  69.82547 2146.9044 300 1000
    ## 11    5 2005 T1_148 Pinus halepensis 38.139968 831.10600  157.6704 100  300
    ## 12    5 2005 T3_168     Quercus ilex  1.486753  76.94320 1967.9696 300 1000
    ## 13    6 2006 T1_148 Pinus halepensis 38.259999 837.32226  157.4934 100  300
    ## 14    6 2006 T3_168     Quercus ilex  1.604515  84.11653 1814.7860 300 1000
    ## 15    7 2007 T1_148 Pinus halepensis 38.379886 843.49782  157.3151 100  300
    ## 16    7 2007 T3_168     Quercus ilex  1.722913  91.33151 1682.5028 300 1000
    ## 17    8 2008 T1_148 Pinus halepensis 38.499589 849.63093  157.1349 100  300
    ## 18    8 2008 T3_168     Quercus ilex  1.841896  98.58787 1567.1993 300 1000
    ## 19    9 2009 T1_148 Pinus halepensis 38.619082 855.72055  156.9537 100  300
    ## 20    9 2009 T3_168     Quercus ilex  1.961392 105.88115 1465.8968 300 1000
    ## 21   10 2010 T1_148 Pinus halepensis 38.738350 861.76624  156.7722 100  300
    ## 22   10 2010 T3_168     Quercus ilex  2.081283 113.20418 1376.2995 300 1000
    ##    Z100 Age ObsID
    ## 1    NA  40  <NA>
    ## 2    NA  24  <NA>
    ## 3    NA  40  <NA>
    ## 4    NA  24  <NA>
    ## 5    NA  41  <NA>
    ## 6    NA  24  <NA>
    ## 7    NA  42  <NA>
    ## 8    NA  25  <NA>
    ## 9    NA  43  <NA>
    ## 10   NA  26  <NA>
    ## 11   NA  44  <NA>
    ## 12   NA  27  <NA>
    ## 13   NA  45  <NA>
    ## 14   NA  28  <NA>
    ## 15   NA  46  <NA>
    ## 16   NA  29  <NA>
    ## 17   NA  47  <NA>
    ## 18   NA  30  <NA>
    ## 19   NA  48  <NA>
    ## 20   NA  31  <NA>
    ## 21   NA  49  <NA>
    ## 22   NA  32  <NA>

## References

- De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini
  M, García-Valdés R, Nadal-Sala D, Sabaté S, Martin-StPaul N, Morin X,
  D’Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A
  trait-enabled model to simulate Mediterranean forest function and
  dynamics at regional scales. Geoscientific Model Development 16:
  3165-3201 (<https://doi.org/10.5194/gmd-16-3165-2023>).
