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
    ## 2         551.3673          25.14821           802.6546             37.59933
    ## 3         550.7311          25.26317           805.2995             37.64858
    ## 4         550.0914          25.37829           807.9405             37.69788
    ## 5         549.4465          25.49349           810.5784             37.74722
    ## 6         548.7999          25.60897           813.2125             37.79660
    ## 7         548.1497          25.72455           815.8400             37.84597
    ## 8         547.4961          25.84040           818.4711             37.89551
    ## 9         546.8371          25.95623           821.0971             37.94507
    ## 10        546.1764          26.07217           823.7118             37.99453
    ## 11        545.5158          26.18865           826.3341             38.04424
    ##    QuadraticMeanTreeDiameter HartBeckingIndex ShrubCoverLive BasalAreaDead
    ## 1                   24.02949         53.20353       3.750000    0.00000000
    ## 2                   24.09840         53.05799       3.987443    0.03899827
    ## 3                   24.16736         52.91426       4.242769    0.03939620
    ## 4                   24.23644         52.77196       4.513862    0.03979730
    ## 5                   24.30564         52.63108       4.801629    0.04031222
    ## 6                   24.37497         52.49150       5.107654    0.04061001
    ## 7                   24.44440         52.35346       5.432457    0.04102079
    ## 8                   24.51400         52.21631       5.777138    0.04143505
    ## 9                   24.58368         52.08066       6.142813    0.04196752
    ## 10                  24.65342         51.94673       6.531184    0.04227339
    ## 11                  24.72338         51.81323       6.943181    0.04246254
    ##    ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1     0.000000000            0             0
    ## 2     0.005927133            0             0
    ## 3     0.006307349            0             0
    ## 4     0.006710839            0             0
    ## 5     0.007159461            0             0
    ## 6     0.007595700            0             0
    ## 7     0.008079318            0             0
    ## 8     0.008592576            0             0
    ## 9     0.009163093            0             0
    ## 10    0.009716999            0             0
    ## 11    0.010272117            0             0

Species-level analogous statistics are shown using:

``` r

fd$SpeciesSummary
```

    ##    Step           Species NumCohorts TreeDensityLive TreeBasalAreaLive
    ## 1     0  Pinus halepensis          1        168.0000         18.604547
    ## 2     0 Quercus coccifera          1              NA                NA
    ## 3     0      Quercus ilex          1        384.0000          6.428755
    ## 4     1  Pinus halepensis          1        167.6997         18.620115
    ## 5     1 Quercus coccifera          1              NA                NA
    ## 6     1      Quercus ilex          1        383.6676          6.528097
    ## 7     2  Pinus halepensis          1        167.3977         18.635316
    ## 8     2 Quercus coccifera          1              NA                NA
    ## 9     2      Quercus ilex          1        383.3334          6.627851
    ## 10    3  Pinus halepensis          1        167.0940         18.650250
    ## 11    3 Quercus coccifera          1              NA                NA
    ## 12    3      Quercus ilex          1        382.9974          6.728043
    ## 13    4  Pinus halepensis          1        166.7878         18.664839
    ## 14    4 Quercus coccifera          1              NA                NA
    ## 15    4      Quercus ilex          1        382.6587          6.828650
    ## 16    5  Pinus halepensis          1        166.4807         18.679251
    ## 17    5 Quercus coccifera          1              NA                NA
    ## 18    5      Quercus ilex          1        382.3191          6.929722
    ## 19    6  Pinus halepensis          1        166.1720         18.693340
    ## 20    6 Quercus coccifera          1              NA                NA
    ## 21    6      Quercus ilex          1        381.9778          7.031209
    ## 22    7  Pinus halepensis          1        165.8614         18.707294
    ## 23    7 Quercus coccifera          1              NA                NA
    ## 24    7      Quercus ilex          1        381.6346          7.133104
    ## 25    8  Pinus halepensis          1        165.5484         18.720853
    ## 26    8 Quercus coccifera          1              NA                NA
    ## 27    8      Quercus ilex          1        381.2887          7.235378
    ## 28    9  Pinus halepensis          1        165.2344         18.734091
    ## 29    9 Quercus coccifera          1              NA                NA
    ## 30    9      Quercus ilex          1        380.9419          7.338075
    ## 31   10  Pinus halepensis          1        164.9205         18.747457
    ## 32   10 Quercus coccifera          1              NA                NA
    ## 33   10      Quercus ilex          1        380.5953          7.441190
    ##    ShrubCoverLive BasalAreaDead ShrubCoverDead BasalAreaCut ShrubCoverCut
    ## 1              NA   0.000000000             NA            0            NA
    ## 2        3.750000            NA    0.000000000           NA             0
    ## 3              NA   0.000000000             NA            0            NA
    ## 4              NA   0.033341905             NA            0            NA
    ## 5        3.987443            NA    0.005927133           NA             0
    ## 6              NA   0.005656365             NA            0            NA
    ## 7              NA   0.033617772             NA            0            NA
    ## 8        4.242769            NA    0.006307349           NA             0
    ## 9              NA   0.005778426             NA            0            NA
    ## 10             NA   0.033895166             NA            0            NA
    ## 11       4.513862            NA    0.006710839           NA             0
    ## 12             NA   0.005902133             NA            0            NA
    ## 13             NA   0.034268142             NA            0            NA
    ## 14       4.801629            NA    0.007159461           NA             0
    ## 15             NA   0.006044074             NA            0            NA
    ## 16             NA   0.034455340             NA            0            NA
    ## 17       5.107654            NA    0.007595700           NA             0
    ## 18             NA   0.006154671             NA            0            NA
    ## 19             NA   0.034737372             NA            0            NA
    ## 20       5.432457            NA    0.008079318           NA             0
    ## 21             NA   0.006283420             NA            0            NA
    ## 22             NA   0.035021194             NA            0            NA
    ## 23       5.777138            NA    0.008592576           NA             0
    ## 24             NA   0.006413857             NA            0            NA
    ## 25             NA   0.035403544             NA            0            NA
    ## 26       6.142813            NA    0.009163093           NA             0
    ## 27             NA   0.006563975             NA            0            NA
    ## 28             NA   0.035593479             NA            0            NA
    ## 29       6.531184            NA    0.009716999           NA             0
    ## 30             NA   0.006679912             NA            0            NA
    ## 31             NA   0.035684550             NA            0            NA
    ## 32       6.943181            NA    0.010272117           NA             0
    ## 33             NA   0.006777986             NA            0            NA

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
    ## 3     1 2001 T1_158 Pinus halepensis 37.59933 802.6546 167.6997 100  300   NA
    ## 4     1 2001 T2_179     Quercus ilex 14.71875 661.8336 383.6676 300 1000   NA
    ## 5     2 2002 T1_158 Pinus halepensis 37.64858 805.2995 167.3977 100  300   NA
    ## 6     2 2002 T2_179     Quercus ilex 14.83724 663.6604 383.3334 300 1000   NA
    ## 7     3 2003 T1_158 Pinus halepensis 37.69788 807.9405 167.0940 100  300   NA
    ## 8     3 2003 T2_179     Quercus ilex 14.95552 665.4814 382.9974 300 1000   NA
    ## 9     4 2004 T1_158 Pinus halepensis 37.74722 810.5784 166.7878 100  300   NA
    ## 10    4 2004 T2_179     Quercus ilex 15.07359 667.2968 382.6587 300 1000   NA
    ## 11    5 2005 T1_158 Pinus halepensis 37.79660 813.2125 166.4807 100  300   NA
    ## 12    5 2005 T2_179     Quercus ilex 15.19147 669.1074 382.3191 300 1000   NA
    ## 13    6 2006 T1_158 Pinus halepensis 37.84597 815.8400 166.1720 100  300   NA
    ## 14    6 2006 T2_179     Quercus ilex 15.30915 670.9123 381.9778 300 1000   NA
    ## 15    7 2007 T1_158 Pinus halepensis 37.89551 818.4711 165.8614 100  300   NA
    ## 16    7 2007 T2_179     Quercus ilex 15.42661 672.7116 381.6346 300 1000   NA
    ## 17    8 2008 T1_158 Pinus halepensis 37.94507 821.0971 165.5484 100  300   NA
    ## 18    8 2008 T2_179     Quercus ilex 15.54385 674.5051 381.2887 300 1000   NA
    ## 19    9 2009 T1_158 Pinus halepensis 37.99453 823.7118 165.2344 100  300   NA
    ## 20    9 2009 T2_179     Quercus ilex 15.66090 676.2935 380.9419 300 1000   NA
    ## 21   10 2010 T1_158 Pinus halepensis 38.04424 826.3341 164.9205 100  300   NA
    ## 22   10 2010 T2_179     Quercus ilex 15.77773 678.0761 380.5953 300 1000   NA
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
    ## 1     1 2001 T1_158 Pinus halepensis 37.59933 802.6546 0.3002897            0
    ## 2     1 2001 T2_179     Quercus ilex 14.71875 661.8336 0.3324344            0
    ## 3     2 2002 T1_158 Pinus halepensis 37.64858 805.2995 0.3019825            0
    ## 4     2 2002 T2_179     Quercus ilex 14.83724 663.6604 0.3342054            0
    ## 5     3 2003 T1_158 Pinus halepensis 37.69788 807.9405 0.3036785            0
    ## 6     3 2003 T2_179     Quercus ilex 14.95552 665.4814 0.3359820            0
    ## 7     4 2004 T1_158 Pinus halepensis 37.74722 810.5784 0.3062180            0
    ## 8     4 2004 T2_179     Quercus ilex 15.07359 667.2968 0.3386932            0
    ## 9     5 2005 T1_158 Pinus halepensis 37.79660 813.2125 0.3070868            0
    ## 10    5 2005 T2_179     Quercus ilex 15.19147 669.1074 0.3395589            0
    ## 11    6 2006 T1_158 Pinus halepensis 37.84597 815.8400 0.3087932            0
    ## 12    6 2006 T2_179     Quercus ilex 15.30915 670.9123 0.3413534            0
    ## 13    7 2007 T1_158 Pinus halepensis 37.89551 818.4711 0.3105027            0
    ## 14    7 2007 T2_179     Quercus ilex 15.42661 672.7116 0.3431535            0
    ## 15    8 2008 T1_158 Pinus halepensis 37.94507 821.0971 0.3130733            0
    ## 16    8 2008 T2_179     Quercus ilex 15.54385 674.5051 0.3459073            0
    ## 17    9 2009 T1_158 Pinus halepensis 37.99453 823.7118 0.3139340            0
    ## 18    9 2009 T2_179     Quercus ilex 15.66090 676.2935 0.3467747            0
    ## 19   10 2010 T1_158 Pinus halepensis 38.04424 826.3341 0.3139154            0
    ## 20   10 2010 T2_179     Quercus ilex 15.77773 678.0761 0.3466743            0
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
    ## 1    1 2001 T1_158 Pinus halepensis 37.59933 802.6546   9.153782 100  300   NA
    ## 2    1 2001 T2_179     Quercus ilex 14.71875 661.8336 383.667566 300 1000   NA
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
    ## 3     1 2001 T1_158 Pinus halepensis 37.59933 802.65462  158.5459 100  300   NA
    ## 4     1 2001 T3_179     Quercus ilex  1.00000  47.23629 3000.0000 300 1000   NA
    ## 5     2 2002 T1_158 Pinus halepensis 37.71708 808.89659  158.3750 100  300   NA
    ## 6     2 2002 T3_179     Quercus ilex  1.00000  47.23629 2998.2955 300 1000   NA
    ## 7     3 2003 T1_158 Pinus halepensis 37.80700 813.64114  158.2032 100  300   NA
    ## 8     3 2003 T3_179     Quercus ilex  1.00000  47.23629 2996.5828 300 1000   NA
    ## 9     4 2004 T1_158 Pinus halepensis 37.85616 816.22709  158.0305 100  300   NA
    ## 10    4 2004 T3_179     Quercus ilex  1.00000  47.23629 2994.8616 300 1000   NA
    ## 11    5 2005 T1_158 Pinus halepensis 37.90534 818.80798  157.8581 100  300   NA
    ## 12    5 2005 T3_179     Quercus ilex  1.00000  47.23629 2993.1431 300 1000   NA
    ## 13    6 2006 T1_158 Pinus halepensis 37.95449 821.38197  157.6855 100  300   NA
    ## 14    6 2006 T3_179     Quercus ilex  1.00000  47.23629 2991.4227 300 1000   NA
    ## 15    7 2007 T1_158 Pinus halepensis 38.00370 823.95269  157.5128 100  300   NA
    ## 16    7 2007 T3_179     Quercus ilex  1.00000  47.23629 2989.7004 300 1000   NA
    ## 17    8 2008 T1_158 Pinus halepensis 38.05286 826.51547  157.3395 100  300   NA
    ## 18    8 2008 T3_179     Quercus ilex  1.00000  47.23629 2987.9714 300 1000   NA
    ## 19    9 2009 T1_158 Pinus halepensis 38.10203 829.07273  157.1664 100  300   NA
    ## 20    9 2009 T3_179     Quercus ilex  1.00000  47.23629 2986.2453 300 1000   NA
    ## 21   10 2010 T1_158 Pinus halepensis 38.15120 831.62458  156.9942 100  300   NA
    ## 22   10 2010 T3_179     Quercus ilex  1.00000  47.23629 2984.5267 300 1000   NA
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
