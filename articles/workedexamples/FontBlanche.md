# Soil and plant water balances at Font-Blanche

## Introduction

### About this vignette

This document describes how to run the water balance model on a forest
plot at Font-Blanche (France), using the R function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md)
included in package `medfate`. The document indicates how to prepare the
model inputs, use the model simulation function, evaluate the
predictions against available observations and inspect the outputs.

### About the Font-Blanche research forest

The Font-Blanche research forest, located in southeastern France
(43º14′27″ N 5°40′45″ E, 420 m elevation), is composed of a top strata
of *Pinus halepensis* (Aleppo pine) reaching about 12 m, a lower strata
of *Quercus ilex* (holm oak), reaching about 6 m, and an understorey
strata dominated by *Quercus coccifera* but including other species such
as *Phillyrea latifolia*. It is spatially heterogeneous: not all trees
in each strata are contiguous, so trees from the lower stratas are
partially exposed to direct light. The forest grows on rocky and shallow
soils that have a low retention capacity and are of Jurassic limestone
origin. The climate is Mediterranean, with a water stress period in
summer, cold or mild winters and most precipitation occurring between
September and May. The experimental site, which is dedicated to study
forest carbon and water cycles, has an enclosed area of 80×80 m (Simioni
et al. 2013) but our specific plot is a quadrat of dimensions 25x25 m.

## Model inputs

Any forest water balance model needs information on **climate**,
**vegetation** and **soils** of the forest stand to be simulated.
Moreover, since the soil water balance in `medfate` differentiates
between species, **species-specific parameters** are also needed. Since
FontBlanche is one of the sites used for evaluating the model, and much
of the data can be found in Moreno et al. (2021). We can use a data list
`fb` with all the necessary inputs:

``` r

fb <- readRDS(paste0("~/OneDrive/mcaceres_work/model_development/medfate_evaluation/StandLevelEvaluation/data/input/FONBLA/FONBLA_", as.character(packageVersion("medfate")), ".rds"))
```

### Soil

We require information on the physical attributes of soil in
Font-Blanche, namely *soil depth*, *texture*, *bulk density* and *rock
fragment content*. Soil information needs to be entered as a
`data frame` with soil layers in rows and physical attributes in
columns. The model accepts one to five soil layers with arbitrary
widths. Because soil properties vary strongly at fine spatial scales,
ideally soil physical attributes should be measured on samples taken at
the forest stand to be simulated. For those users lacking such data,
soil properties modelled at larger scales are available via
soilgrids.org (see function `soilgridsParams()`). In our case soil
physical attributes are already defined in the data bundled for
FontBlanche:

``` r

spar <- fb$soilData
print(spar)
```

    ##   widths clay sand om   bd rfc
    ## 1    300   39   26  6 1.45  50
    ## 2    700   39   26  3 1.45  65
    ## 3   1000   39   26  1 1.45  90
    ## 4   2500   39   26  1 1.45  95

The soil input for function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md) is
actually an object of class `soil` that is created using a function with
the same name:

``` r

fb_soil <- soil(spar)
```

The [`print()`](https://rdrr.io/r/base/print.html) function for objects
`soil` provides a lot of information on soil physical properties and
water capacity:

``` r

print(fb_soil)
```

    ##   widths sand clay      usda om nitrogen ph   bd rfc   macro     Ksat VG_alpha
    ## 1    300   26   39 Clay loam  6       NA NA 1.45  50 0.07387 7232.425 44.14586
    ## 2    700   26   39 Clay loam  3       NA NA 1.45  65 0.07387 3481.917 61.34088
    ## 3   1000   26   39 Clay loam  1       NA NA 1.45  90 0.07387 1879.187 76.38182
    ## 4   2500   26   39 Clay loam  1       NA NA 1.45  95 0.07387 1879.187 76.38182
    ##       VG_n VG_theta_res VG_theta_sat W Temp
    ## 1 1.254346        0.041    0.4388377 1   NA
    ## 2 1.273896        0.041    0.4388377 1   NA
    ## 3 1.287757        0.041    0.4388377 1   NA
    ## 4 1.287757        0.041    0.4388377 1   NA

The soil object is also used to store the moisture degree of each soil
layer. In particular, `W` contains the state variable that represents
moisture content - the proportion of moisture **relative to field
capacity** - which is by default initialized to 1 for each layer:

``` r

fb_soil$W
```

    ## [1] 1 1 1 1

### Species parameters

Simulation models in `medfate` require a data frame with species
parameter values. The package provides a default data set of parameter
values for a number of Mediterranean species occurring in Spain (rows),
resulting from bibliographic search, fit to empirical data or
expert-based guesses:

``` r

data("SpParamsMED")
```

However, sometimes one may wish to override species defaults with custom
values. In the case of FontBlanche there is a table of preferred
parameters:

``` r

fb$customParams
```

    ##               Species GrowthForm Hmax Hmed  SLA Kmax_stemxylem LeafEPS LeafPI0
    ## 1 Phillyrea latifolia Tree/Shrub   NA   NA   NA             NA   12.38   -2.13
    ## 2    Pinus halepensis       Tree   NA   NA   NA           0.15    5.31   -1.50
    ## 3        Quercus ilex       Tree 1240  450 4.55           0.40   15.00   -2.50
    ##   LeafAF StemEPS StemPI0 StemAF    Al2As
    ## 1    0.5   12.38   -2.13    0.4       NA
    ## 2    0.6    5.00   -1.65    0.4  631.000
    ## 3    0.4   15.00   -2.50    0.4 1540.671

We can use function
[`modifySpParams()`](https://emf-creaf.github.io/medfate/reference/modifyParams.md)
to replace the values of parameters for the desired traits, leaving the
rest unaltered:

``` r

SpParamsFB <- modifySpParams(traits4models::SpParamsFR, fb$customParams)
SpParamsFB
```

    ##                     Name SpIndex        AcceptedName             Species
    ## 2574 Phillyrea latifolia    2573 Phillyrea latifolia Phillyrea latifolia
    ## 2631    Pinus halepensis    2630    Pinus halepensis    Pinus halepensis
    ## 2853        Quercus ilex    2853        Quercus ilex        Quercus ilex
    ##          Genus   Family    Order      Group GrowthForm     LifeForm LeafShape
    ## 2574 Phillyrea Oleaceae Lamiales Angiosperm Tree/Shrub Phanerophyte     Broad
    ## 2631     Pinus Pinaceae  Pinales Gymnosperm       Tree Phanerophyte    Needle
    ## 2853   Quercus Fagaceae  Fagales Angiosperm       Tree Phanerophyte     Broad
    ##      LeafSize      PhenologyType SenescencePeriod GrowthPeriod
    ## 2574   Medium oneflush-evergreen               NA           NA
    ## 2631    Small oneflush-evergreen               NA           NA
    ## 2853   Medium oneflush-evergreen               NA           NA
    ##      GrowthDeterminacy DispersalType Hmed      Hmax Dmax Z50       Z95 fHDmin
    ## 2574                NA    vertebrate 1500  850.7389   NA  NA  746.1321     NA
    ## 2631                NA          wind 2000 2200.0000   NA  NA 7500.0000     NA
    ## 2853                NA    vertebrate  450 1240.0000  200  NA 7000.0000     NA
    ##      fHDmax       a_ash    b_ash     a_bsh     b_bsh    a_btsh    b_btsh cr
    ## 2574     NA 0.008047873 2.865021 0.2502494 0.7300031 0.3105815 0.7509837 NA
    ## 2631     NA          NA       NA        NA        NA        NA        NA NA
    ## 2853     NA 1.857486249 1.885548 0.5238830 0.7337293 0.7327147 0.7375770 NA
    ##      BTsh      a_fbt    b_fbt       c_fbt a_cr b_1cr b_2cr b_3cr c_1cr c_2cr
    ## 2574   NA         NA       NA          NA   NA    NA    NA    NA    NA    NA
    ## 2631   NA 0.07607828 1.462411 -0.02280106   NA    NA    NA    NA    NA    NA
    ## 2853   NA 0.07848713 1.497670 -0.00309341   NA    NA    NA    NA    NA    NA
    ##           a_cw  b_cw      a_bt      b_bt LeafDuration BudFormationDays t0gdd
    ## 2574        NA    NA        NA        NA     2.500000               NA    NA
    ## 2631 0.6415296 0.731 0.5535741 1.1848613     2.083333               NA    NA
    ## 2853        NA    NA 0.5622245 0.9626839     2.000000               NA    NA
    ##      Sgdd Tbgdd Ugdd Ssen Phsen Tbsen xsen ysen Fsen     SLA LeafDensity
    ## 2574   NA    NA   NA   NA    NA    NA   NA   NA   NA 5.40000     0.56000
    ## 2631   NA    NA   NA   NA    NA    NA   NA   NA   NA 4.99002     0.28812
    ## 2853   NA    NA   NA   NA    NA    NA   NA   NA   NA 4.55000     0.66500
    ##      WoodDensity FineRootDensity conduit2sapwood     r635   pDead    Al2As
    ## 2574        0.70              NA              NA 1.917579      NA 1698.950
    ## 2631        0.60              NA              NA 1.964226 0.00050  631.000
    ## 2853        0.89              NA              NA 1.805872 0.00026 1540.671
    ##      Ar2Al LeafWidth      SRL RLD   maxFMC   minFMC  Ptlp LeafPI0 LeafEPS
    ## 2574    NA      0.40 902.1171  NA 106.2825 41.20000    NA   -2.13   12.38
    ## 2631    NA      0.10       NA  NA 122.1644 53.74550 -2.20   -1.50    5.31
    ## 2853    NA      1.63 735.7025  NA 109.2304 58.66087 -3.09   -2.50   15.00
    ##      LeafAF StemPI0 StemEPS StemAF  SAV HeatContent LeafLigninPercent
    ## 2574    0.5   -2.13   12.38    0.4   NA    17146.03                NA
    ## 2631    0.6   -1.65    5.00    0.4 6050    22150.00                NA
    ## 2853    0.4   -2.50   15.00    0.4 4050    19300.00                NA
    ##      WoodLigninPercent FineRootLigninPercent LeafAngle LeafAngleSD
    ## 2574                NA                    NA  35.14508    18.87379
    ## 2631                NA                    NA        NA          NA
    ## 2853                NA                    NA  36.00819    20.15953
    ##      ClumpingIndex gammaSWR alphaSWR kPAR  g Tmax_LAI Tmax_LAIsq Psi_Extract
    ## 2574            NA       NA       NA   NA NA       NA         NA          NA
    ## 2631            NA       NA       NA   NA NA       NA         NA          NA
    ## 2853            NA       NA       NA   NA NA       NA         NA          NA
    ##      Exp_Extract WUE WUE_par WUE_co2 WUE_vpd       Gswmin     Gswmax
    ## 2574          NA  NA      NA      NA      NA 0.0039671810 0.45220000
    ## 2631          NA  NA      NA      NA      NA 0.0006655459 0.04554533
    ## 2853          NA  NA      NA      NA      NA 0.0028095533 0.21575000
    ##      Gsw_Toptim_Jarvis Gsw_Tsens_Jarvis Gsw_AC_slope_Baldocchi
    ## 2574                NA               NA                     NA
    ## 2631                NA               NA                     NA
    ## 2853                NA               NA                     NA
    ##      Gsw_P50_Baldocchi Gsw_slope_Baldocchi VCleaf_kmax VCleaf_P12 VCleaf_P50
    ## 2574                NA                  NA          NA         NA         NA
    ## 2631                NA                  NA          NA         NA         NA
    ## 2853                NA                  NA    9.276068         NA       -3.5
    ##      VCleaf_P88 VCleaf_slope VCleafapo_P12 VCleafapo_P50 VCleafapo_P88
    ## 2574         NA           NA            NA            NA            NA
    ## 2631         NA           NA            NA            NA            NA
    ## 2853         NA           NA            NA            NA            NA
    ##      VCleafapo_slope Kmax_stemxylem VCstem_P12 VCstem_P50 VCstem_P88
    ## 2574              NA      0.4083769   -3.06267      -6.55  -12.17000
    ## 2631              NA      0.1500000   -4.50000      -5.14   -6.24103
    ## 2853              NA      0.4000000   -4.30000      -6.88   -7.90000
    ##      VCstem_slope Kmax_rootxylem VCroot_P12 VCroot_P50 VCroot_P88 VCroot_slope
    ## 2574     8.344927       1.889495   -3.12248      -5.30   -7.47752     17.45105
    ## 2631    43.652321       0.890000   -0.51103      -0.88   -1.24897    102.98940
    ## 2853    21.111111       0.470000   -0.43456      -2.39   -7.03000     11.52311
    ##       Vmax298  Jmax298   Nleaf Nsapwood Nfineroot  WoodC RERleaf RERsapwood
    ## 2574       NA       NA 15.7000     2.53        NA     NA      NA         NA
    ## 2631       NA       NA 11.7940       NA  12.90000 0.4981      NA         NA
    ## 2853 39.44837 93.30728 13.9021     5.63   2.45857 0.4751      NA         NA
    ##      RERfineroot CCleaf CCsapwood CCfineroot RGRleafmax RGRsapwoodmax
    ## 2574          NA   1.63        NA         NA         NA            NA
    ## 2631          NA     NA        NA         NA         NA            NA
    ## 2853          NA   1.43        NA         NA         NA            NA
    ##      RGRcambiummax RGRfinerootmax RGRbud SRsapwood SRfineroot RSSG
    ## 2574            NA             NA     NA        NA         NA   NA
    ## 2631            NA             NA     NA        NA         NA 0.25
    ## 2853            NA             NA     NA        NA         NA 0.25
    ##      MortalityBaselineRate SurvivalModelStep SurvivalB0 SurvivalB1
    ## 2574                    NA                NA         NA         NA
    ## 2631                    NA                NA         NA         NA
    ## 2853                    NA                NA         NA         NA
    ##      SeedProductionHeight SeedProductionDiameter  SeedMass SeedLongevity
    ## 2574                   NA                     NA   31.5244            NA
    ## 2631                   NA                     NA   19.6100            10
    ## 2853                   NA               10.64702 3225.8100            NA
    ##      DispersalDistance DispersalShape ProbRecr MinTempRecr MinMoistureRecr
    ## 2574                NA             NA       NA          NA              NA
    ## 2631                NA             NA       NA          NA              NA
    ## 2853                NA             NA       NA          NA              NA
    ##      MinFPARRecr RecrAge RecrTreeDBH RecrTreeHeight RecrShrubHeight
    ## 2574          NA      NA          NA             NA              NA
    ## 2631          NA      NA          NA             NA              NA
    ## 2853          NA      NA          NA             NA              NA
    ##      RecrTreeDensity RecrShrubCover RespFire RespDist RespClip
    ## 2574              NA             NA       NA       NA       NA
    ## 2631              NA             NA       NA       NA       NA
    ## 2853              NA             NA       NA       NA       NA
    ##      IngrowthTreeDensity IngrowthTreeDBH
    ## 2574                  NA              NA
    ## 2631                  NA              NA
    ## 2853                  NA              NA

Note that the function returns a subset of rows for the species
mentioned in `customParams`. Not all parameters are needed for the soil
water balance model. The user can find parameter definitions in the help
page of this data set. However, to fully understand the role of
parameters in the model, the user should read the details of model
design and formulation (<http://emf-creaf.github.io/medfate>).

### Vegetation

Models included in `medfate` were primarily designed to be ran on
**forest inventory plots**. In this kind of data, the vegetation of a
sampled area is described in terms of woody plants (trees and shrubs)
along with their size and species identity. Forest plots in `medfate`
are assumed to be in a format that follows closely the Spanish forest
inventory. Each forest plot is represented in an object of class
`forest`, a list that contains several elements. Among them, the most
important items are two data frames, `treeData` (for trees) and
`shrubData` for shrubs:

``` r

fb_forest <- emptyforest()
fb_forest
```

    ## $treeData
    ## [1] Species DBH     Height  N       Z50     Z95    
    ## <0 rows> (or 0-length row.names)
    ## 
    ## $shrubData
    ## [1] Species Height  Cover   Z50     Z95    
    ## <0 rows> (or 0-length row.names)
    ## 
    ## attr(,"class")
    ## [1] "forest" "list"

Trees are expected to be primarily described in terms of species,
diameter (DBH) and height, whereas shrubs are described in terms of
species, percent cover and mean height. In our case, we will for
simplicity avoid shrubs and concentrate on the main three tree species
in the Font-Blanche forest plot: *Phillyrea latifolia* (code 142),
*Pinus halepensis* (Alepo pine, code 148), and *Quercus ilex* (holm oak;
code 168). In order to run the model, one has to prepare a data table
like this one, already prepared for Font-Blanche:

``` r

fb$treeData
```

    ##               Species       DBH    Height    N Z50  Z95       LAI
    ## 1 Phillyrea latifolia  2.587859  323.0000 1248 200 1500 0.2238497
    ## 2    Pinus halepensis 26.759914 1195.7667  256 200  800 1.0767397
    ## 3        Quercus ilex  6.220031  495.5532 3104 300 2000 1.3994106

Trees have been grouped by species, so DBH and height values are means
(in cm), and `N` indicates the number of trees in each category. Package
`medfate` allows separating trees by size, but for simplicity we do not
distinguish here between tree sizes within each species. Columns `Z50`
and `Z95` indicate the depths (in mm) corresponding to cumulative 50%
and 95% of fine roots, respectively.

In order to use this data, we need to replace the part corresponding to
trees into the forest object that we created before:

``` r

fb_forest$treeData <- fb$treeData
fb_forest
```

    ## $treeData
    ##               Species       DBH    Height    N Z50  Z95       LAI
    ## 1 Phillyrea latifolia  2.587859  323.0000 1248 200 1500 0.2238497
    ## 2    Pinus halepensis 26.759914 1195.7667  256 200  800 1.0767397
    ## 3        Quercus ilex  6.220031  495.5532 3104 300 2000 1.3994106
    ## 
    ## $shrubData
    ## [1] Species Height  Cover   Z50     Z95    
    ## <0 rows> (or 0-length row.names)
    ## 
    ## attr(,"class")
    ## [1] "forest" "list"

Because the forest plot format is rather specific, `medfate` also allows
starting in an alternative way using two data frames, one with
**aboveground** information (i.e. the leave area and size of plants) and
the other with **belowground** information (i.e. root distribution). The
aboveground data frame does not distinguish between trees and shrubs. It
includes, for each plant cohort to be considered in rows, its *species
identity*, *height*, *leaf area index* (LAI) and *crown ratio*. While
users can build their input data themselves, we use internal function
[`forest2aboveground()`](https://emf-creaf.github.io/medfate/reference/forest2aboveground.md)
on the object `fb_forest` to show how should the data look like:

``` r

fb_above <- forest2aboveground(fb_forest, SpParamsFB)
fb_above
```

    ##           SP    N       DBH Cover         H        CR  LAI_live LAI_expanded
    ## T1_2573 2573 1248  2.587859    NA  323.0000 0.5510653 0.2238497    0.2238497
    ## T2_2630 2630  256 26.759914    NA 1195.7667 0.6001765 1.0767397    1.0767397
    ## T3_2853 2853 3104  6.220031    NA  495.5532 0.6059123 1.3994106    1.3994106
    ##         LAI_dead LAI_nocomp LAI_mistletoe Age ObsID
    ## T1_2573        0  0.2238497             0  NA  <NA>
    ## T2_2630        0  1.0767397             0  NA  <NA>
    ## T3_2853        0  1.3994106             0  NA  <NA>

Note that the call to
[`forest2aboveground()`](https://emf-creaf.github.io/medfate/reference/forest2aboveground.md)
included species parameters, because species-specific parameter values
are needed to calculate leaf area from tree diameters or shrub cover
using allometric relationships. Columns `N`, `DBH` and `Cover` are
required for simulating growth, but not for soil water balance, which
only requires columns `SP`, `H` (in cm), `CR` (i.e. the crown ratio),
`LAI_live`, `LAI_expanded` and `LAI_dead`. Here plant cohorts are given
unique codes that tell us whether they correspond to trees or shrubs. In
practice, the user only needs to worry to calculate the values for
`LAI_live`. `LAI_live` and `LAI_expanded` can contain the same LAI
values, and `LAI_dead` is normally zero.

We see that at Font-Blanche holm oaks (code 68) represent most of the
total leaf area. On the other hand, pines are taller than oaks.
`medfate` assumes that leaf distribution follows a truncated normal
curve between the crown base height and the total height. This can be
easily inspected using function
[`vprofile_leafAreaDensity()`](https://emf-creaf.github.io/medfate/reference/vprofile_leafAreaDensity.md):

``` r

vprofile_leafAreaDensity(fb_forest, SpParamsFB, byCohorts = T, bySpecies = T)
```

![](FontBlanche_files/figure-html/unnamed-chunk-12-1.png)

Regarding **belowground** information, the usuer should supply a matrix
describing for each plant cohort, the proportion of fine roots in each
soil layer. As before, we use internal function
[`forest2belowground()`](https://emf-creaf.github.io/medfate/reference/forest2aboveground.md)
on the object `fb_forest` to show how should the data look like:

``` r

fb_below <- forest2belowground(fb_forest, fb_soil, SpParamsFB)
fb_below
```

    ##                 1         2          3           4
    ## T1_2573 0.6505881 0.2719302 0.05418875 0.023292883
    ## T2_2630 0.7035917 0.2658370 0.02440189 0.006169504
    ## T3_2853 0.5075225 0.3714797 0.08507611 0.035921747

In our case, these proportions were implicitly specified in parameters
`Z50` and `Z95`. In fact, these values describe a continuous
distribution of fine roots along depth, which can be displayed using
function
[`vprofile_rootDistribution()`](https://emf-creaf.github.io/medfate/reference/vprofile_leafAreaDensity.md):

``` r

vprofile_rootDistribution(fb_forest, SpParamsFB, bySpecies = T)
```

![](FontBlanche_files/figure-html/unnamed-chunk-14-1.png)

Note that in Font-Blanche we set that the root system of Aleppo pines
(*Pinus halepensis*) would be more superficial than that of the other
two species. Moreover, holm oak trees are the ones who extend their
roots down to deepest soil layers.

### Meteorology

Water balance simulations of function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md)
require **daily weather** inputs. The weather variables that are
required depend on the complexity of the soil water balance model we are
using. In the simplest case, only **mean temperature**,
**precipitation** and **potential evapo-transpiration (PET)** is
required, but the more complex simulation model also requires radiation,
wind speed, min/max temparature and relative humitidy. Here we already
have a data frame with the daily meteorology measured at Font-Blanche
for year 2014:

``` r

fb_meteo <- fb$meteoData
head(fb_meteo)
```

    ##        dates MeanTemperature MinTemperature MaxTemperature MeanRelativeHumidity
    ## 1 2014-01-01        7.661856       5.988889       8.960000             87.78224
    ## 2 2014-01-02        9.525431       7.958333      11.550000             96.40669
    ## 3 2014-01-03        9.482417       8.176111      11.762220             93.05705
    ## 4 2014-01-04       10.016813       6.313000      11.010000             96.31667
    ## 5 2014-01-05        6.619919       4.766000       9.060555             57.77938
    ## 6 2014-01-06        8.923008       6.793889      12.329440             64.40477
    ##   MinRelativeHumidity MaxRelativeHumidity WindSpeed Precipitation Radiation
    ## 1            80.37265            98.48404  2.317495      0.000000 1.5050178
    ## 2            84.22588           100.00000  2.407691      0.000000 2.6173102
    ## 3            79.93501           100.00000  1.950114      0.000000 3.9089762
    ## 4            90.14023           100.00000  3.596797      2.590674 0.4753025
    ## 5            48.92043            65.71329  7.310334      0.000000 8.6224570
    ## 6            51.31975            74.46718  2.301697      0.000000 6.7835715

Simulation models in `medfate` have been designed to work along with
data generated from package `meteoland` (De Cáceres et al. 2018), which
specifies conventions for variable names and units. The user is strongly
recommended to resort to this package to obtain suitable weather input
for soil water balance simulations (see
<http://emf-creaf.github.io/meteoland>).

### Simulation control

Apart from data inputs, the behavior of simulation models can be
controlled using a set of **global parameters**. The default global
parameter values are obtained using function
[`defaultControl()`](https://emf-creaf.github.io/medfate/reference/defaultControl.md):

``` r

fb_control <- defaultControl(transpirationMode = "Sureau")
fb_control$subdailyResults <- TRUE
```

Where the following changes are set to control parameters:

1.  Transpiration is set `transpirationMode = "Sureau"`, which implies a
    greater complexity of plant hydraulics and energy balance
    calculations.
2.  Subdaily results generated by the model are kept.

### Water balance input object

A last step is needed before calling simulation functions. It consists
in the compilation of all aboveground and belowground parameters and the
specification of additional parameter values for each plant cohort, such
as their light extinction coefficient or their response to drought. If
one has a `forest` object, the `spwbInput` object can be generated in
directly from it, avoiding the need to explicitly build `fb_above` and
`fb_below` data frames:

``` r

fb_x <- spwbInput(fb_forest, fb_soil, SpParamsFB, fb_control)
```

Different species parameter variables will be drawn from `SpParamsMED`
depending on the value of `transpirationMode`. For the simple water
balance model, relatively few parameters are needed. All the input
information for forest data and species parameter values can be
inspected by printing the input object.

Finally, note that one can play with plant-specific parameters for soil
water balance (instead of using species-level values) by using function
[`modifyCohortParams()`](https://emf-creaf.github.io/medfate/reference/modifyParams.md).

## Running the model

Function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md)
requires two main objects as input:

- A `spwbInput` object with forest and soil parameters (`fb_x` in our
  case).
- A data frame with daily meteorology for the study period (`fb_meteo`
  in our case).

Now we are ready to call function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md):

``` r

fb_SWB <- spwb(fb_x, fb_meteo, elevation = 420, latitude = 43.24083)
```

    ## Initial plant water content (mm): 35.1976
    ## Initial soil water content (mm): 213.886
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2014]:............
    ##  [Year 2015]:............
    ##  [Year 2016]:............
    ##  [Year 2017]:............
    ##  [Year 2018]:............
    ## 
    ## Final plant water content (mm): 34.5749
    ## Final soil water content (mm): 250.173
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.622733
    ## Plant water balance result (mm): -7.24793
    ## Change in soil water content (mm): 36.2872
    ## Soil water balance result (mm): 36.2872
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): 0
    ## Water balance components:
    ##   Precipitation (mm) 3638 Rain (mm) 3638 Snow (mm) 1
    ##   Interception (mm) 568 Net rainfall (mm) 3070 Runon (mm) 0
    ##   Infiltration (mm) 2735 Infiltration excess (mm) 336 Saturation excess (mm) 426 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 130  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 1211  Mistletoe transpiration (mm) 0
    ##   Plant extraction from soil (mm) 1204  Plant water balance (mm) -7 Hydraulic redistribution (mm) 83
    ##   Runoff (mm) 762 Deep drainage (mm) 938

Console output provides the water balance totals for the period
considered, which may span several years. The output of function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md) is an
object of class with the same name, actually a list:

``` r

class(fb_SWB)
```

    ## [1] "spwb" "list"

If we inspect its elements, we realize that there are several
components:

``` r

names(fb_SWB)
```

    ##  [1] "latitude"      "topography"    "weather"       "spwbInput"    
    ##  [5] "spwbOutput"    "WaterBalance"  "EnergyBalance" "Temperature"  
    ##  [9] "Soil"          "Snow"          "Stand"         "Plants"       
    ## [13] "SunlitLeaves"  "ShadeLeaves"   "subdaily"

For example, `WaterBalance` contains water balance components in form of
a data frame with days in rows:

``` r

head(fb_SWB$WaterBalance)
```

    ##                  PET Precipitation     Rain Snow   NetRain Snowmelt Runon
    ## 2014-01-01 0.6209989      0.000000 0.000000    0 0.0000000        0     0
    ## 2014-01-02 0.5671238      0.000000 0.000000    0 0.0000000        0     0
    ## 2014-01-03 0.5418115      0.000000 0.000000    0 0.0000000        0     0
    ## 2014-01-04 0.6072565      2.590674 2.590674    0 0.6850664        0     0
    ## 2014-01-05 2.0447148      0.000000 0.000000    0 0.0000000        0     0
    ## 2014-01-06 0.9330456      0.000000 0.000000    0 0.0000000        0     0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2014-01-01    0.0000000                  0                0      0    0.0000000
    ## 2014-01-02    0.0000000                  0                0      0    0.0000000
    ## 2014-01-03    0.0000000                  0                0      0    0.0000000
    ## 2014-01-04    0.6850664                  0                0      0    0.1625982
    ## 2014-01-05    0.0000000                  0                0      0    0.0000000
    ## 2014-01-06    0.0000000                  0                0      0    0.0000000
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2014-01-01               0          0.2397383     0.000000       0.2151227
    ## 2014-01-02               0          0.2109763     0.000000       0.1964596
    ## 2014-01-03               0          0.2246855     0.000000       0.1876911
    ## 2014-01-04               0          2.0794742     1.905607       0.1738434
    ## 2014-01-05               0          0.6088971     0.000000       0.2720512
    ## 2014-01-06               0          0.4067717     0.000000       0.1506880
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2014-01-01                 0    0.0234958585  2.461560e-02
    ## 2014-01-02                 0    0.0144841266  1.451665e-02
    ## 2014-01-03                 0    0.0357446793  3.699440e-02
    ## 2014-01-04                 0   -0.0002881711  2.363583e-05
    ## 2014-01-05                 0    0.3205646215  3.368459e-01
    ## 2014-01-06                 0    0.2588755337  2.560838e-01
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2014-01-01                      0            0.0000000000
    ## 2014-01-02                      0            0.0007454724
    ## 2014-01-03                      0            0.0019186949
    ## 2014-01-04                      0            0.0056910988
    ## 2014-01-05                      0            0.0000000000
    ## 2014-01-06                      0            0.0000000000

## Comparing results with observations

Before examining the results of the model, it is important to compare
its predictions against observed data, if available. The following
observations are available from the experimental forest plot for year
2014:

- Stand total evapotranspiration estimated using an Eddy-covariance flux
  tower.
- Soil moisture content of the first 0-30 cm layer.
- Cohort transpiration estimates derived from sapflow measurements
  for Q. ilex and P. halepensis.
- Pre-dawn and midday leaf water potentials for Q. ilex and P.
  halepensis.

We first load the measured data into the workspace and filter for the
dates used in the simulation:

``` r

fb_observed <- fb$measuredData
fb_observed <- fb_observed[fb_observed$dates %in% fb_meteo$dates,]
row.names(fb_observed) <- fb_observed$dates
head(fb_observed)
```

    ##                 dates       SWC E_T2_2630 E_T2_2630_err E_T3_2853 E_T3_2853_err
    ## 2014-01-01 2014-01-01 0.5813407        NA            NA        NA            NA
    ## 2014-01-02 2014-01-02 0.6507478        NA            NA        NA            NA
    ## 2014-01-03 2014-01-03 0.6224243        NA            NA        NA            NA
    ## 2014-01-04 2014-01-04        NA        NA            NA        NA            NA
    ## 2014-01-05 2014-01-05 0.6285134        NA            NA        NA            NA
    ## 2014-01-06 2014-01-06 0.6035415        NA            NA        NA            NA
    ##            PD_T2_2630 PD_T2_2630_err PD_T3_2853 PD_T3_2853_err MD_T2_2630
    ## 2014-01-01         NA             NA         NA             NA         NA
    ## 2014-01-02         NA             NA         NA             NA         NA
    ## 2014-01-03         NA             NA         NA             NA         NA
    ## 2014-01-04         NA             NA         NA             NA         NA
    ## 2014-01-05         NA             NA         NA             NA         NA
    ## 2014-01-06         NA             NA         NA             NA         NA
    ##            MD_T2_2630_err MD_T3_2853 MD_T3_2853_err LE  H
    ## 2014-01-01             NA         NA             NA NA NA
    ## 2014-01-02             NA         NA             NA NA NA
    ## 2014-01-03             NA         NA             NA NA NA
    ## 2014-01-04             NA         NA             NA NA NA
    ## 2014-01-05             NA         NA             NA NA NA
    ## 2014-01-06             NA         NA             NA NA NA

### Stand evapotranspiration

Package `medfate` contains several functions to assist the evaluation of
model results. We can first compare the observed vs modelled latent
heat. We can plot the two time series:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "LE", plotType="dynamics")+
  theme(legend.position = c(0.8,0.85))
```

![](FontBlanche_files/figure-html/unnamed-chunk-23-1.png)

The relationship can be shown in a scatter plot:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "LE", plotType="scatter")
```

![](FontBlanche_files/figure-html/unnamed-chunk-24-1.png) Where we see
that the model tends to underestimate total evapotranspiration during
seasons with low evaporative demand. Function
[`evaluation_stats()`](https://emf-creaf.github.io/medfate/reference/evaluation.md)
allows us to generate evaluation statistics:

``` r

evaluation_stats(fb_SWB, fb_observed, type = "LE")
```

    ##            n         Bias     Bias.rel          MAE      MAE.rel            r 
    ## 1026.0000000   -0.5730345  -19.5440356    1.7302957   59.0138280    0.3384201 
    ##          NSE      NSE.abs 
    ##   -0.3454001   -0.1439333

### Soil moisture

We can compare observed vs modelled soil moisture content in a similar
way as we did for total evapotranspiration:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "SWC", plotType="dynamics")
```

![](FontBlanche_files/figure-html/unnamed-chunk-26-1.png)

As before, we can generate a scatter plot:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "SWC", plotType="scatter")
```

![](FontBlanche_files/figure-html/unnamed-chunk-27-1.png)

or examine evaluation statistics:

``` r

evaluation_stats(fb_SWB, fb_observed, type = "SWC")
```

    ##             n          Bias      Bias.rel           MAE       MAE.rel 
    ## 1760.00000000   -0.13894266  -31.17594665    0.13937653   31.27329772 
    ##             r           NSE       NSE.abs 
    ##    0.88091861   -0.25512156   -0.09045666

### Plant transpiration

The following plots display the observed and predicted transpiration
dynamics for *Pinus halepensis* and *Quercus ilex*:

``` r

g1<-evaluation_plot(fb_SWB, fb_observed, 
                            cohort = "T2_2630",
                            type="E", plotType = "dynamics")+
  theme(legend.position = c(0.85,0.83))
g2<-evaluation_plot(fb_SWB, fb_observed, 
                            cohort = "T3_2853",
                            type="E", plotType = "dynamics")+
  theme(legend.position = c(0.85,0.83))
plot_grid(g1, g2, ncol=1)
```

![](FontBlanche_files/figure-html/unnamed-chunk-29-1.png)

In general, the agreement is quite good, but the model seems to
overestimate the transpiration of *P. halepensis* in early summer and
after the first drought period. The transpiration of *Q. ilex* seems
also overestimated in spring and autumn. We can also inspect the
evaluation statistics for both species using:

``` r

evaluation_stats(fb_SWB, fb_observed, cohort = "T2_2630", type="E")
```

    ##             n          Bias      Bias.rel           MAE       MAE.rel 
    ## 300.000000000   0.001567047   0.761916820   0.142187085  69.133051472 
    ##             r           NSE       NSE.abs 
    ##   0.199853778  -1.510387442  -0.436873670

``` r

evaluation_stats(fb_SWB, fb_observed, cohort = "T3_2853", type="E")
```

    ##            n         Bias     Bias.rel          MAE      MAE.rel            r 
    ## 309.00000000   0.04343765  15.00651604   0.10913508  37.70317115   0.84898882 
    ##          NSE      NSE.abs 
    ##   0.16604351   0.25648566

### Leaf water potentials

Finally, we can compare observed with predicted water potentials. In
this case measurements are available for three dates, but they include
the standard deviation of several measurements.

``` r

g1<-evaluation_plot(fb_SWB, fb_observed, 
                            cohort = "T2_2630",
                            type="WP", plotType = "dynamics")+
  theme(legend.position = c(0.85,0.23))
g2<-evaluation_plot(fb_SWB, fb_observed, 
                            cohort = "T3_2853",
                            type="WP", plotType = "dynamics")+
  theme(legend.position = c(0.85,0.23))
plot_grid(g1, g2, ncol=1)
```

![](FontBlanche_files/figure-html/unnamed-chunk-31-1.png)

The model seems to overestimate water potentials (i.e. it predicts less
negative values than those observed) for P. halepensis during the
drought season.

## Drawing plots

Package `medfate` provides a simple `plot` function for objects of class
`spwb`. Here we will use this function to display the seasonal variation
predicted by the model, as well as the variation at higher temporal
resolution within four different selected 3-day periods that we define
here:

``` r

d1 = seq(as.Date("2014-03-01"), as.Date("2014-03-03"), by="day")
d2 = seq(as.Date("2014-06-01"), as.Date("2014-06-03"), by="day")
d3 = seq(as.Date("2014-08-01"), as.Date("2014-08-03"), by="day")
d4 = seq(as.Date("2014-10-01"), as.Date("2014-10-03"), by="day")
```

### Meteorological input and input/output water flows

Function [`plot()`](https://rdrr.io/r/graphics/plot.default.html) can be
used to show the meteorological input:

``` r

plot(fb_SWB, type = "PET_Precipitation")
```

![](FontBlanche_files/figure-html/unnamed-chunk-33-1.png) It is apparent
the climatic drought period between april and august 2014. This should
have an impact on soil moisture and plant stress.

If we are interested in forest hydrology, we can plot the amount of
water that the model predicts to leave the forest via surface runoff or
drainage to lower water compartments.

``` r

plot(fb_SWB, type = "Export")
```

![](FontBlanche_files/figure-html/unnamed-chunk-34-1.png) As expected,
water exported from the forest plot was only relevant for the autumn and
winter periods. Note also that the model predicts some runoff during
convective storms during autumn, whereas winter events occur when the
soil is already full, so that most exported water is assumed to be lost
via deep drainage. One can also display the evapotranspiration flows,
which we do in the following plot that also combines the two previous:

``` r

g1<-plot(fb_SWB)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = "none")
g2<-plot(fb_SWB, "Evapotranspiration")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.13,0.73))
g3<-plot(fb_SWB, "Export")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.35,0.60))
plot_grid(g1,g2, g3, ncol=1, rel_heights = c(0.4,1,0.6))
```

![](FontBlanche_files/figure-html/unnamed-chunk-35-1.png)

### Soil moisture dynamics and hydraulic redistribution

It is also useful to plot the dynamics of soil state variables by layer,
such as the percentage of moisture in relation to field capacity:

``` r

plot(fb_SWB, type="SoilTheta")
```

![](FontBlanche_files/figure-html/unnamed-chunk-36-1.png) Note that the
model predicts soil drought to occur earlier in the season for the first
three layers (0-200 cm) whereas the bottom layer dries out much more
slowly. At this point is important to mention that the water balance
model incorporates. We can also display the dynamics of the
corresponding soil layer water potentials:

``` r

plot(fb_SWB, type="SoilPsi")
```

![](FontBlanche_files/figure-html/unnamed-chunk-37-1.png) or draw a
composite plot including absolute soil water volume:

``` r

g1<-plot(fb_SWB)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = "none")
g2<-plot(fb_SWB, "SoilVol")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.08,0.65))
g3<-plot(fb_SWB, "SoilPsi")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.08,0.5))
plot_grid(g1, g2,  g3, rel_heights = c(0.4,0.8,0.8), ncol=1)
```

![](FontBlanche_files/figure-html/unnamed-chunk-38-1.png)

### Root water uptake and hydraulic redistribution

The following composite plot shows the daily root water uptake (or
release) from different soil layers, and the daily amount of water
entering soil layers due to hydraulic redistribution:

``` r

g1<-plot(fb_SWB, "SoilPsi")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = "none")+ylab("Soil wp (MPa)")
g2<-plot(fb_SWB, "PlantExtraction")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.08,0.68))
g3<-plot(fb_SWB, "HydraulicRedistribution")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.08,0.5))
plot_grid(g1, g2,  g3, rel_heights = c(0.4,0.8,0.8), ncol=1)
```

![](FontBlanche_files/figure-html/unnamed-chunk-39-1.png)

If we create a composite plot including subdaily water uptake/release
patterns, we can further understand the redistribution flows created by
the model during different periods:

``` r

g0<-plot(fb_SWB, "PlantExtraction")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.08,0.68))
g1<-plot(fb_SWB, "PlantExtraction", subdaily = T, dates = d1)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylim(c(-0.05,0.13))
g2<-plot(fb_SWB, "PlantExtraction", subdaily = T, dates = d2)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(-0.05,0.13))
g3<-plot(fb_SWB, "PlantExtraction", subdaily = T, dates = d3)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(-0.05,0.13))
g4<-plot(fb_SWB, "PlantExtraction", subdaily = T, dates = d4)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(-0.05,0.13))
plot_grid(g0,plot_grid(g1, g2, g3, g4, ncol=4),ncol=1)
```

![](FontBlanche_files/figure-html/unnamed-chunk-40-1.png)

### Plant transpiration

We can use function
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to display the
seasonal dynamics of cohort-level variables, such as plant transpiration
per leaf area:

``` r

par(mar=c(5,5,1,1))
plot(fb_SWB, type="TranspirationPerLeaf", bySpecies = T)
```

![](FontBlanche_files/figure-html/unnamed-chunk-41-1.png) Where we can
observe that some species transpire more than others due to their
vertical position within the canopy.

``` r

g1<-plot(fb_SWB)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = "none")
g2<-plot(fb_SWB, "TranspirationPerLeaf", bySpecies = T)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.1,0.75))
g21<-plot(fb_SWB, "LeafTranspiration", subdaily = T, dates = d1)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylim(c(0,0.32))
g22<-plot(fb_SWB, "LeafTranspiration", subdaily = T, dates = d2)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(0,0.32))
g23<-plot(fb_SWB, "LeafTranspiration", subdaily = T, dates = d3)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(0,0.32))
g24<-plot(fb_SWB, "LeafTranspiration", subdaily = T, dates = d4)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(0,0.32))
plot_grid(g1, g2,  
          plot_grid(g21,g22,g23,g24, ncol=4), 
          ncol=1, rel_heights = c(0.4,0.8,0.8))
```

![](FontBlanche_files/figure-html/unnamed-chunk-42-1.png)

### Plant stress

In the model, reduction of (whole-plant) plant transpiration is what
used to define drought stress, which depends on the species identity:

``` r

plot(fb_SWB, type="PlantStress", bySpecies = T)
```

![](FontBlanche_files/figure-html/unnamed-chunk-43-1.png)

To examine the impact of drought on plants, one can inspect the
whole-plant conductance (from which the stress index is derived) or the
stem percent loss of conductance derived from embolism, as we do in the
following composite plot:

``` r

g1<-plot(fb_SWB)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = "none")
g2<-plot(fb_SWB, "SoilPlantConductance", bySpecies = T)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+
  ylab(expression(paste("Soil-plant conductance ",(mmol%.%m^{-2}%.%s^{-1}))))+
  theme(legend.position = "none")
g3<-plot(fb_SWB, "StemPLC", bySpecies = T)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.2,0.75))
plot_grid(g1, g2,g3,                          
          ncol=1, rel_heights = c(0.4,0.8,0.8))
```

![](FontBlanche_files/figure-html/unnamed-chunk-44-1.png)

### Leaf water potentials

``` r

g1<-plot(fb_SWB)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = "none")
g2<-plot(fb_SWB, "LeafPsiRange", bySpecies = T)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.1,0.25)) + ylab("Leaf water potential (MPa)")
g21<-plot(fb_SWB, "LeafPsi", subdaily = T, dates = d1)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylim(c(-7,0))
g22<-plot(fb_SWB, "LeafPsi", subdaily = T, dates = d2)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(-7,0))
g23<-plot(fb_SWB, "LeafPsi", subdaily = T, dates = d3)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(-7,0))
g24<-plot(fb_SWB, "LeafPsi", subdaily = T, dates = d4)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(-7,0))
plot_grid(g1, g2,                          
          plot_grid(g21,g22,g23,g24, ncol=4), 
          ncol=1, rel_heights = c(0.4,0.8,0.8))
```

![](FontBlanche_files/figure-html/unnamed-chunk-45-1.png)

### Stomatal conductance

``` r

g1<-plot(fb_SWB)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = "none")
g2<-plot(fb_SWB, "GSWMax_SL", bySpecies = T)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.5,0.74))+ylab("Sunlit leaf stomatal conductance")+ylim(c(0,0.3))
g21<-plot(fb_SWB, "LeafStomatalConductance", subdaily = T, dates = d1)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylim(c(0,0.2))
g22<-plot(fb_SWB, "LeafStomatalConductance", subdaily = T, dates = d2)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(0,0.2))
g23<-plot(fb_SWB, "LeafStomatalConductance", subdaily = T, dates = d3)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(0,0.2))
g24<-plot(fb_SWB, "LeafStomatalConductance", subdaily = T, dates = d4)+scale_x_datetime(date_breaks = "1 day",  date_labels = "%m/%d")+theme(legend.position = "none")+ylab("")+ylim(c(0,0.2))
plot_grid(g1, g2,
          plot_grid(g21,g22,g23,g24, ncol=4),
          ncol=1, rel_heights = c(0.4,0.8,0.8))
```

![](FontBlanche_files/figure-html/unnamed-chunk-46-1.png)

## Generating output summaries

While the water balance model operates at daily and sub-daily steps,
users will normally be interested in outputs at larger time scales. The
package provides a `summary` for objects of class `spwb`. This function
can be used to summarize the model’s output at different temporal steps
(i.e. weekly, monthly or annual). For example, to obtain the average
soil moisture and water potentials by months one can use:

``` r

summary(fb_SWB, freq="months",FUN=sum, output="WaterBalance")
```

    ##                  PET Precipitation      Rain Snow      NetRain Snowmelt Runon
    ## 2014-01-01  27.03414     205.04814 205.04814  0.0 182.21891719      0.0     0
    ## 2014-02-01  37.11592     181.09641 181.09641  0.0 154.86636936      0.0     0
    ## 2014-03-01  80.49737      44.61248  44.61248  0.0  39.82086698      0.0     0
    ## 2014-04-01 109.24874      15.00000  15.00000  0.0   7.20513386      0.0     0
    ## 2014-05-01 147.99639      21.60000  21.60000  0.0  16.29630898      0.0     0
    ## 2014-06-01 167.27898      33.60000  33.60000  0.0  25.81618662      0.0     0
    ## 2014-07-01 183.99299       0.60000   0.60000  0.0   0.14352019      0.0     0
    ## 2014-08-01 159.66330      60.40000  60.40000  0.0  52.81961926      0.0     0
    ## 2014-09-01 103.42793     137.60000 137.60000  0.0 125.87946030      0.0     0
    ## 2014-10-01  63.53896      50.60000  50.60000  0.0  41.91508112      0.0     0
    ## 2014-11-01  30.12083     222.60000 222.60000  0.0 198.01806495      0.0     0
    ## 2014-12-01  26.01617      93.00000  93.00000  0.0  78.70784400      0.0     0
    ## 2015-01-01  33.33412      97.40000  97.40000  0.0  86.41058179      0.0     0
    ## 2015-02-01  38.17936      67.60000  67.60000  0.0  47.44367914      0.0     0
    ## 2015-03-01  78.27482      53.00000  53.00000  0.0  36.73646326      0.0     0
    ## 2015-04-01 109.20859      41.20000  41.20000  0.0  35.27149099      0.0     0
    ## 2015-05-01 157.73664       0.80000   0.80000  0.0   0.19716630      0.0     0
    ## 2015-06-01 176.08186      71.80000  71.80000  0.0  58.24747382      0.0     0
    ## 2015-07-01 199.29900       0.20000   0.20000  0.0   0.04929158      0.0     0
    ## 2015-08-01 158.13941       6.40000   6.40000  0.0   4.16694789      0.0     0
    ## 2015-09-01 107.67144      50.40000  50.40000  0.0  41.03981869      0.0     0
    ## 2015-10-01  54.97221      94.00000  94.00000  0.0  79.01195147      0.0     0
    ## 2015-11-01  40.87647      32.20000  32.20000  0.0  26.16340247      0.0     0
    ## 2015-12-01  18.36149      23.80000  23.80000  0.0  12.52936595      0.0     0
    ## 2016-01-01  25.98382      21.40000  21.40000  0.0  13.12334576      0.0     0
    ## 2016-02-01  42.43210      87.80000  87.80000  0.0  74.31832849      0.0     0
    ## 2016-03-01  72.70349      19.40000  19.40000  0.0  11.66800343      0.0     0
    ## 2016-04-01 110.36442      13.00000  13.00000  0.0   5.24821513      0.0     0
    ## 2016-05-01 151.45145      50.20000  50.20000  0.0  37.50894960      0.0     0
    ## 2016-06-01 168.69569       8.20000   8.20000  0.0   5.52895478      0.0     0
    ## 2016-07-01 194.66192      41.40000  41.40000  0.0  36.33128037      0.0     0
    ## 2016-08-01 167.03376       6.00000   6.00000  0.0   2.07730727      0.0     0
    ## 2016-09-01 113.75042      87.40000  87.40000  0.0  79.12256195      0.0     0
    ## 2016-10-01  59.38548     109.80000 109.80000  0.0  95.79241647      0.0     0
    ## 2016-11-01  31.62855     154.40000 154.40000  0.0 140.17216856      0.0     0
    ## 2016-12-01  25.74714      16.60000  16.60000  0.0  13.13569173      0.0     0
    ## 2017-01-01  30.97017      51.20000  51.20000  0.0  40.79283860      0.0     0
    ## 2017-02-01  39.00523      27.20000  27.20000  0.0  13.28325570      0.0     0
    ## 2017-03-01  75.39796      59.00000  59.00000  0.0  48.83006602      0.0     0
    ## 2017-04-01 115.11745      71.00000  71.00000  0.0  59.70351533      0.0     0
    ## 2017-05-01 148.68280      21.20000  21.20000  0.0  13.30109554      0.0     0
    ## 2017-06-01 178.51316      14.20000  14.20000  0.0  11.11308991      0.0     0
    ## 2017-07-01 200.88791       0.60000   0.60000  0.0   0.15707887      0.0     0
    ## 2017-08-01 171.08660       1.40000   1.40000  0.0   0.39112601      0.0     0
    ## 2017-09-01 109.87609      19.80000  19.80000  0.0  14.97238705      0.0     0
    ## 2017-10-01  75.77870       0.40000   0.40000  0.0   0.13767760      0.0     0
    ## 2017-11-01  40.14167      48.00000  48.00000  0.0  45.87929621      0.0     0
    ## 2017-12-01  26.10486      57.40000  57.40000  0.0  45.82343306      0.0     0
    ## 2018-01-01  34.17212      68.60000  68.60000  0.0  59.21116058      0.0     0
    ## 2018-02-01  32.99313      27.60000  27.00000  0.6  16.20940385      0.0     0
    ## 2018-03-01  65.00212     105.60000 105.60000  0.0  87.01457364      0.6     0
    ## 2018-04-01 106.80480      89.40000  89.40000  0.0  80.65872821      0.0     0
    ## 2018-05-01 110.14567     101.20000 101.20000  0.0  82.74485073      0.0     0
    ## 2018-06-01 149.84895      51.60000  51.60000  0.0  38.49169593      0.0     0
    ## 2018-07-01 188.24784      12.00000  12.00000  0.0   9.41068515      0.0     0
    ## 2018-08-01 161.32739      89.00000  89.00000  0.0  82.57776081      0.0     0
    ## 2018-09-01 113.71008       1.00000   1.00000  0.0   0.34419461      0.0     0
    ## 2018-10-01  61.59182     299.40000 299.40000  0.0 279.91630762      0.0     0
    ## 2018-11-01  28.64739     153.00000 153.00000  0.0 132.08348013      0.0     0
    ## 2018-12-01  25.83913      48.20000  48.20000  0.0  41.85688078      0.0     0
    ##            Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2014-01-01 182.21891719           0.000000        90.353928  90.353928
    ## 2014-02-01 154.86636936           0.000000       122.592068 122.592068
    ## 2014-03-01  39.82086698           0.000000         0.000000   0.000000
    ## 2014-04-01   7.20513386           0.000000         0.000000   0.000000
    ## 2014-05-01  16.29630898           0.000000         0.000000   0.000000
    ## 2014-06-01  25.81618662           0.000000         0.000000   0.000000
    ## 2014-07-01   0.14352019           0.000000         0.000000   0.000000
    ## 2014-08-01  43.75153702           9.068082         0.000000   9.068082
    ## 2014-09-01  98.01236828          27.867092         0.000000  27.867092
    ## 2014-10-01  31.30131212          10.613769         0.000000  10.613769
    ## 2014-11-01 144.50610950          53.511955        21.263979  74.775934
    ## 2014-12-01  78.70784400           0.000000        55.160132  55.160132
    ## 2015-01-01  86.41058179           0.000000         6.017072   6.017072
    ## 2015-02-01  47.44367914           0.000000         3.760259   3.760259
    ## 2015-03-01  36.73646326           0.000000         0.000000   0.000000
    ## 2015-04-01  35.27149099           0.000000         0.000000   0.000000
    ## 2015-05-01   0.19716630           0.000000         0.000000   0.000000
    ## 2015-06-01  58.24747382           0.000000         0.000000   0.000000
    ## 2015-07-01   0.04929158           0.000000         0.000000   0.000000
    ## 2015-08-01   4.16694789           0.000000         0.000000   0.000000
    ## 2015-09-01  37.94683452           3.092984         0.000000   3.092984
    ## 2015-10-01  73.15756356           5.854388         0.000000   5.854388
    ## 2015-11-01  24.84460899           1.318793         0.000000   1.318793
    ## 2015-12-01  12.52936595           0.000000         0.000000   0.000000
    ## 2016-01-01  13.12334576           0.000000         0.000000   0.000000
    ## 2016-02-01  74.31832849           0.000000         0.000000   0.000000
    ## 2016-03-01  11.66800343           0.000000         0.000000   0.000000
    ## 2016-04-01   5.24821513           0.000000         0.000000   0.000000
    ## 2016-05-01  37.50894960           0.000000         0.000000   0.000000
    ## 2016-06-01   5.52895478           0.000000         0.000000   0.000000
    ## 2016-07-01  34.03901702           2.292263         0.000000   2.292263
    ## 2016-08-01   2.07730727           0.000000         0.000000   0.000000
    ## 2016-09-01  72.44848565           6.674076         0.000000   6.674076
    ## 2016-10-01  71.67022158          24.122195         0.000000  24.122195
    ## 2016-11-01  97.95997067          42.212198         9.359735  51.571933
    ## 2016-12-01  13.13569173           0.000000         0.000000   0.000000
    ## 2017-01-01  40.79283860           0.000000         0.000000   0.000000
    ## 2017-02-01  13.28325570           0.000000         0.000000   0.000000
    ## 2017-03-01  48.83006602           0.000000         0.000000   0.000000
    ## 2017-04-01  59.70351533           0.000000         0.000000   0.000000
    ## 2017-05-01  13.30109554           0.000000         0.000000   0.000000
    ## 2017-06-01  11.11308991           0.000000         0.000000   0.000000
    ## 2017-07-01   0.15707887           0.000000         0.000000   0.000000
    ## 2017-08-01   0.39112601           0.000000         0.000000   0.000000
    ## 2017-09-01  14.97238705           0.000000         0.000000   0.000000
    ## 2017-10-01   0.13767760           0.000000         0.000000   0.000000
    ## 2017-11-01  39.69193096           6.187365         0.000000   6.187365
    ## 2017-12-01  45.82343306           0.000000         0.000000   0.000000
    ## 2018-01-01  59.21116058           0.000000         0.000000   0.000000
    ## 2018-02-01  16.20940385           0.000000         0.000000   0.000000
    ## 2018-03-01  87.61457364           0.000000         0.000000   0.000000
    ## 2018-04-01  80.65872821           0.000000        17.042554  17.042554
    ## 2018-05-01  82.74485073           0.000000         7.907588   7.907588
    ## 2018-06-01  38.49169593           0.000000         0.000000   0.000000
    ## 2018-07-01   9.41068515           0.000000         0.000000   0.000000
    ## 2018-08-01  63.98773626          18.590025         0.000000  18.590025
    ## 2018-09-01   0.34419461           0.000000         0.000000   0.000000
    ## 2018-10-01 187.34522936          92.571078        30.899218 123.470296
    ## 2018-11-01 100.15513457          31.928346        52.634563  84.562908
    ## 2018-12-01  41.85688078           0.000000         9.305614   9.305614
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2014-01-01   27.4057795               0          30.545864   22.8292233
    ## 2014-02-01   40.1478462               0          34.172384   26.2300390
    ## 2014-03-01   44.4494011               0          25.505837    4.7916094
    ## 2014-04-01   19.4817664               0          31.051741    7.7948661
    ## 2014-05-01    0.0000000               0          40.321502    5.3036910
    ## 2014-06-01    0.0000000               0          57.183630    7.7838134
    ## 2014-07-01    0.0000000               0          45.164377    0.4564798
    ## 2014-08-01    0.0000000               0          42.344913    7.5803807
    ## 2014-09-01    0.5718093               0          35.779103   11.7205387
    ## 2014-10-01   13.5379689               0          24.740167    8.6849189
    ## 2014-11-01   38.7139945               0          32.780150   24.5819351
    ## 2014-12-01   44.4494011               0          21.997162   14.2921560
    ## 2015-01-01   44.4403098               0          19.886230   10.9894182
    ## 2015-02-01   40.1478462               0          28.175949   20.1563209
    ## 2015-03-01   44.4494011               0          35.751365   16.2635367
    ## 2015-04-01   29.6271190               0          30.643410    5.9285090
    ## 2015-05-01    6.0416296               0          41.655182    0.6028337
    ## 2015-06-01    3.5757122               0          64.327109   13.5525252
    ## 2015-07-01    0.0000000               0          64.007093    0.1507084
    ## 2015-08-01    0.0000000               0          26.710825    2.2330521
    ## 2015-09-01    0.0000000               0          27.544534    9.3601813
    ## 2015-10-01    0.0000000               0          27.084666   14.9880475
    ## 2015-11-01    0.0000000               0          17.814915    6.0365975
    ## 2015-12-01    0.8493024               0          17.343636   11.2706340
    ## 2016-01-01    3.8750047               0          14.387254    8.2766542
    ## 2016-02-01   13.2169687               0          24.201027   13.4816715
    ## 2016-03-01   35.5085117               0          22.491557    7.7319966
    ## 2016-04-01    0.0000000               0          29.746999    7.7517849
    ## 2016-05-01    0.9277800               0          46.679022   12.6910504
    ## 2016-06-01    0.0000000               0          48.482087    2.6710452
    ## 2016-07-01    0.0000000               0          60.338844    5.0687196
    ## 2016-08-01    0.0000000               0          41.108497    3.9226927
    ## 2016-09-01    0.0000000               0          28.355737    8.2774381
    ## 2016-10-01    6.8741223               0          30.498641   14.0075835
    ## 2016-11-01   14.3503349               0          22.696004   14.2278314
    ## 2016-12-01   44.4494011               0          13.798456    3.4643083
    ## 2017-01-01   19.4269649               0          18.007837   10.4071614
    ## 2017-02-01   27.0994195               0          23.845350   13.9167443
    ## 2017-03-01   23.2834479               0          28.808920   10.1699340
    ## 2017-04-01   34.6040262               0          37.357026   11.2964847
    ## 2017-05-01   11.1844141               0          46.505872    7.8989045
    ## 2017-06-01    0.0000000               0          53.473929    3.0869101
    ## 2017-07-01    0.0000000               0          47.285416    0.4429211
    ## 2017-08-01    0.0000000               0          13.370785    1.0088740
    ## 2017-09-01    0.0000000               0          13.451341    4.8276129
    ## 2017-10-01    0.0000000               0           6.044858    0.2623224
    ## 2017-11-01    0.0000000               0           8.748361    2.1207038
    ## 2017-12-01    0.0000000               0          16.256684   11.5765669
    ## 2018-01-01    2.8590464               0          17.866397    9.3888394
    ## 2018-02-01    3.4265797               0          18.883082   10.7905961
    ## 2018-03-01   43.0076908               0          32.919527   18.5854264
    ## 2018-04-01   43.0155495               0          30.274571    8.7412718
    ## 2018-05-01   44.4494011               0          39.813724   18.4551493
    ## 2018-06-01   43.0155495               0          46.087218   13.1083041
    ## 2018-07-01    8.4547547               0          50.030904    2.5893148
    ## 2018-08-01    0.8268817               0          49.910022    6.4222392
    ## 2018-09-01    0.0000000               0          27.776693    0.6558054
    ## 2018-10-01   29.2860778               0          36.021493   19.4836924
    ## 2018-11-01   43.0155495               0          25.425609   20.9165199
    ## 2018-12-01   44.4494011               0          14.938733    6.3431192
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2014-01-01     4.395375328                 0        3.303881      3.321265
    ## 2014-02-01     3.185813579                 0        4.719379      4.756531
    ## 2014-03-01     3.385384523                 0       17.236817     17.328843
    ## 2014-04-01     0.616410070                 0       22.530781     22.640465
    ## 2014-05-01     0.423976489                 0       34.480101     34.593834
    ## 2014-06-01     0.324957517                 0       48.816930     49.074859
    ## 2014-07-01     0.150206608                 0       43.247614     44.557691
    ## 2014-08-01     0.460383426                 0       34.547733     34.304149
    ## 2014-09-01     0.791594515                 0       23.614032     23.266970
    ## 2014-10-01     3.406499349                 0       12.679192     12.648748
    ## 2014-11-01     3.817777277                 0        4.399523      4.380437
    ## 2014-12-01     3.324426990                 0        4.365928      4.380579
    ## 2015-01-01     2.472811521                 0        6.424608      6.424000
    ## 2015-02-01     2.732162958                 0        5.279627      5.287465
    ## 2015-03-01     5.903042259                 0       13.503769     13.584786
    ## 2015-04-01     1.277873201                 0       23.346379     23.437028
    ## 2015-05-01     0.566715087                 0       40.302349     40.485633
    ## 2015-06-01     2.046971060                 0       48.638701     48.727613
    ## 2015-07-01     0.159123100                 0       62.423496     63.697261
    ## 2015-08-01     0.056742474                 0       23.275901     24.421030
    ## 2015-09-01     0.150328594                 0       19.082950     18.034024
    ## 2015-10-01     2.111280587                 0       10.298164      9.985337
    ## 2015-11-01     3.136961106                 0        8.640883      8.641356
    ## 2015-12-01     3.489278420                 0        2.596816      2.583724
    ## 2016-01-01     2.839387829                 0        3.257899      3.271212
    ## 2016-02-01     4.283692556                 0        6.438564      6.435663
    ## 2016-03-01     2.608324553                 0       12.128338     12.151236
    ## 2016-04-01     0.713061067                 0       21.245979     21.282153
    ## 2016-05-01     1.340183807                 0       32.601959     32.647788
    ## 2016-06-01     0.360519149                 0       45.166491     45.450522
    ## 2016-07-01     0.213987666                 0       54.706472     55.056137
    ## 2016-08-01     0.099789121                 0       35.504008     37.086015
    ## 2016-09-01     1.131238319                 0       19.818338     18.947061
    ## 2016-10-01     4.053347196                 0       12.465053     12.437710
    ## 2016-11-01     3.054173755                 0        5.426885      5.413999
    ## 2016-12-01     3.518071328                 0        6.796048      6.816076
    ## 2017-01-01     2.464223389                 0        5.143252      5.136452
    ## 2017-02-01     3.942807197                 0        5.977291      5.985799
    ## 2017-03-01     3.846811182                 0       14.761615     14.792175
    ## 2017-04-01     2.888016298                 0       23.162441     23.172525
    ## 2017-05-01     2.982928641                 0       35.526808     35.624039
    ## 2017-06-01     0.193531014                 0       49.941331     50.193487
    ## 2017-07-01     0.099218927                 0       45.011136     46.743276
    ## 2017-08-01     0.006720004                 0       10.743977     12.355191
    ## 2017-09-01     0.042844705                 0        8.908071      8.580884
    ## 2017-10-01     0.030398689                 0        5.596645      5.752137
    ## 2017-11-01     0.273196763                 0        7.057287      6.354460
    ## 2017-12-01     1.515509031                 0        3.256899      3.164608
    ## 2018-01-01     3.890886211                 0        4.601138      4.586671
    ## 2018-02-01     4.582227861                 0        3.506711      3.510258
    ## 2018-03-01     7.133102433                 0        7.188199      7.200998
    ## 2018-04-01     2.124557063                 0       19.387099     19.408742
    ## 2018-05-01     4.726506070                 0       16.621853     16.632069
    ## 2018-06-01     2.952209009                 0       29.963061     30.026705
    ## 2018-07-01     0.359105315                 0       46.970397     47.082484
    ## 2018-08-01     1.692129232                 0       41.809653     41.795653
    ## 2018-09-01     0.258883952                 0       26.762034     26.862004
    ## 2018-10-01     5.363311026                 0       11.263003     11.174490
    ## 2018-11-01     1.556495256                 0        2.948527      2.952594
    ## 2018-12-01     4.016478337                 0        4.558093      4.579136
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2014-01-01                      0              0.31402257
    ## 2014-02-01                      0              0.32050834
    ## 2014-03-01                      0              0.16335767
    ## 2014-04-01                      0              0.28005688
    ## 2014-05-01                      0              0.32848982
    ## 2014-06-01                      0              1.01806443
    ## 2014-07-01                      0              1.39514948
    ## 2014-08-01                      0              4.57535546
    ## 2014-09-01                      0              6.49682236
    ## 2014-10-01                      0              1.91431287
    ## 2014-11-01                      0              0.50609181
    ## 2014-12-01                      0              0.39487745
    ## 2015-01-01                      0              0.20845705
    ## 2015-02-01                      0              0.41043723
    ## 2015-03-01                      0              0.24153349
    ## 2015-04-01                      0              0.17136761
    ## 2015-05-01                      0              0.37457496
    ## 2015-06-01                      0              1.00359047
    ## 2015-07-01                      0              2.09016241
    ## 2015-08-01                      0              1.21845756
    ## 2015-09-01                      0              6.19178416
    ## 2015-10-01                      0              6.63344896
    ## 2015-11-01                      0              1.08117711
    ## 2015-12-01                      0              0.88369060
    ## 2016-01-01                      0              2.32060891
    ## 2016-02-01                      0              0.86890341
    ## 2016-03-01                      0              0.06545428
    ## 2016-04-01                      0              0.05286983
    ## 2016-05-01                      0              0.17213671
    ## 2016-06-01                      0              0.31888001
    ## 2016-07-01                      0              4.56115113
    ## 2016-08-01                      0              0.40043597
    ## 2016-09-01                      0              7.56260562
    ## 2016-10-01                      0              3.41729526
    ## 2016-11-01                      0              1.19132688
    ## 2016-12-01                      0              0.17622040
    ## 2017-01-01                      0              0.08181758
    ## 2017-02-01                      0              0.07564421
    ## 2017-03-01                      0              0.05630399
    ## 2017-04-01                      0              0.12509753
    ## 2017-05-01                      0              0.01612544
    ## 2017-06-01                      0              2.55186824
    ## 2017-07-01                      0              1.34651181
    ## 2017-08-01                      0              0.06179886
    ## 2017-09-01                      0              1.46556860
    ## 2017-10-01                      0              0.21950607
    ## 2017-11-01                      0              8.34896075
    ## 2017-12-01                      0              4.87431288
    ## 2018-01-01                      0              1.67207698
    ## 2018-02-01                      0              0.02877751
    ## 2018-03-01                      0              0.05393424
    ## 2018-04-01                      0              0.09771792
    ## 2018-05-01                      0              0.05414739
    ## 2018-06-01                      0              0.03520693
    ## 2018-07-01                      0              0.21304440
    ## 2018-08-01                      0              0.62683303
    ## 2018-09-01                      0              0.64306422
    ## 2018-10-01                      0              0.64412374
    ## 2018-11-01                      0              0.17953758
    ## 2018-12-01                      0              0.25010845

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of the three tree species by months:

``` r

summary(fb_SWB, freq="months",FUN=mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Phillyrea latifolia Pinus halepensis Quercus ilex
    ## 2014-01-01          0.08017912     0.0001081269  0.002361513
    ## 2014-02-01          0.08057793     0.0001081143  0.002370308
    ## 2014-03-01          0.08434290     0.0001287368  0.002590731
    ## 2014-04-01          0.08775242     0.0001399277  0.002791190
    ## 2014-05-01          0.09283954     0.0001652265  0.003212696
    ## 2014-06-01          0.10352502     0.0002077553  0.004283342
    ## 2014-07-01          0.13735414     0.0002800530  0.020944152
    ## 2014-08-01          0.19586620     0.0004633708  0.091462232
    ## 2014-09-01          0.19274398     0.0002173203  0.096859916
    ## 2014-10-01          0.17598415     0.0001231343  0.082479753
    ## 2014-11-01          0.15542776     0.0001099366  0.064647401
    ## 2014-12-01          0.13454457     0.0001089919  0.046548041
    ## 2015-01-01          0.11397309     0.0001145816  0.028591165
    ## 2015-02-01          0.09446188     0.0001095851  0.011564567
    ## 2015-03-01          0.08424985     0.0001223920  0.002635809
    ## 2015-04-01          0.08780937     0.0001405034  0.002796292
    ## 2015-05-01          0.09240752     0.0001723882  0.003246606
    ## 2015-06-01          0.10034424     0.0001998523  0.003776099
    ## 2015-07-01          0.12635667     0.0003840159  0.013170574
    ## 2015-08-01          0.20524421     0.0021233072  0.097714639
    ## 2015-09-01          0.23055652     0.0044442256  0.136794337
    ## 2015-10-01          0.22214022     0.0006657703  0.129879094
    ## 2015-11-01          0.20224650     0.0001208266  0.112366833
    ## 2015-12-01          0.18160607     0.0001097796  0.094241989
    ## 2016-01-01          0.16040905     0.0001107674  0.075889545
    ## 2016-02-01          0.14030209     0.0001164499  0.058524504
    ## 2016-03-01          0.12057524     0.0001230603  0.041181432
    ## 2016-04-01          0.10113507     0.0001398563  0.023932154
    ## 2016-05-01          0.09295211     0.0001561593  0.007952405
    ## 2016-06-01          0.09783889     0.0001935287  0.003715783
    ## 2016-07-01          0.16384628     0.0003757079  0.025488635
    ## 2016-08-01          0.19858970     0.0005953991  0.068786958
    ## 2016-09-01          0.25420538     0.0028655151  0.159886058
    ## 2016-10-01          0.24813428     0.0002056098  0.157771519
    ## 2016-11-01          0.22705493     0.0001131331  0.139378350
    ## 2016-12-01          0.20561383     0.0001160679  0.120664451
    ## 2017-01-01          0.18435827     0.0001150265  0.102007529
    ## 2017-02-01          0.16425773     0.0001156427  0.084413060
    ## 2017-03-01          0.14466582     0.0001273647  0.067130510
    ## 2017-04-01          0.12503489     0.0001402265  0.049616899
    ## 2017-05-01          0.10614800     0.0001562991  0.032470587
    ## 2017-06-01          0.11163171     0.0002797394  0.017472821
    ## 2017-07-01          0.17656627     0.0011841423  0.035165269
    ## 2017-08-01          0.27250108     0.0090528147  0.218412337
    ## 2017-09-01          0.35796610     0.0128921921  0.388595297
    ## 2017-10-01          0.36069009     0.0128467520  0.394831438
    ## 2017-11-01          0.35553154     0.0108946521  0.390598450
    ## 2017-12-01          0.34023369     0.0049882687  0.376635747
    ## 2018-01-01          0.31938202     0.0002829264  0.357219167
    ## 2018-02-01          0.29790934     0.0001125145  0.336996841
    ## 2018-03-01          0.27644625     0.0001156398  0.316758497
    ## 2018-04-01          0.25505554     0.0001343714  0.296329373
    ## 2018-05-01          0.23415147     0.0001236757  0.276273224
    ## 2018-06-01          0.21354820     0.0001446187  0.256478857
    ## 2018-07-01          0.19542508     0.0002038250  0.238018329
    ## 2018-08-01          0.18287143     0.0002071400  0.220658005
    ## 2018-09-01          0.16788373     0.0001918222  0.202746511
    ## 2018-10-01          0.15875245     0.0001433243  0.185179894
    ## 2018-11-01          0.13866157     0.0001051863  0.165886213
    ## 2018-12-01          0.11802898     0.0001110850  0.146607616

In this case, the `summary` function aggregates the output by species
using LAI values as weights.

## Bibliography

- De Caceres M, Martin-StPaul N, Turco M, et al (2018) Estimating daily
  meteorological data and downscaling climate models over landscapes.
  Environ Model Softw 108:186–196.
  <https://doi.org/10.1016/j.envsoft.2018.08.003>

- De Caceres M, Martinez-Vilalta J, Coll L, et al (2015) Coupling a
  water balance model with forest inventory data to predict drought
  stress: the role of forest structural changes vs. climate changes.
  Agric For Meteorol 213:77–90.
  <https://doi.org/10.1016/j.agrformet.2015.06.012>

- Simioni G, Durand-gillmann M, Huc R, et al (2013) Asymmetric
  competition increases leaf inclination effect on light absorption in
  mixed canopies. Ann For Sci 70:123–131.
  <https://doi.org/10.1007/s13595-012-0246-8>

- Moreno, M., Simioni, G., Cailleret, M., Ruffault, J., Badel, E.,
  Carrière, S., Davi, H., Gavinet, J., Huc, R., Limousin, J.-M.,
  Marloie, O., Martin, L., Rodríguez-Calcerrada, J., Vennetier, M.,
  Martin-StPaul, N., 2021. Consistently lower sap velocity and growth
  over nine years of rainfall exclusion in a Mediterranean mixed
  pine-oak forest. Agric. For. Meteorol. 308–309, 108472.
  <https://doi.org/10.1016/j.agrformet.2021.108472>
