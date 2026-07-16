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

fb <- medfatereports::load_list("FONBLA")
```

``` r

fb <- readRDS("fb_data.rds")
names(fb)
```

    ##  [1] "siteData"       "treeData"       "shrubData"      "customParams"  
    ##  [5] "measuredData"   "meteoData"      "miscData"       "soilData"      
    ##  [9] "terrainData"    "remarks"        "sp_params"      "forest_object1"

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

    ##               Species VCleaf_P12 VCleaf_P50 VCleaf_P88 VCleaf_slope VCstem_P12
    ## 1 Phillyrea latifolia         NA         NA         NA           NA  -1.971750
    ## 2    Pinus halepensis         NA         NA         NA           NA  -3.707158
    ## 3        Quercus ilex         NA         NA         NA           NA  -4.739642
    ##   VCstem_P50 VCstem_P88 VCstem_slope VCroot_P12 VCroot_P50 VCroot_P88
    ## 1      -6.50 -11.028250           11         NA         NA         NA
    ## 2      -4.79  -5.872842           46         -1  -1.741565  -2.301482
    ## 3      -6.40  -8.060358           30         NA         NA         NA
    ##   VCroot_slope VCleaf_kmax LeafEPS LeafPI0 LeafAF StemEPS StemPI0 StemAF Gswmin
    ## 1           NA        3.00   12.38   -2.13    0.5   12.38   -2.13    0.4  0.002
    ## 2           NA        4.00    5.31   -1.50    0.6    5.00   -1.65    0.4  0.001
    ## 3           NA        2.63   15.00   -2.50    0.4   15.00   -2.50    0.4  0.002
    ##   Gswmax    Gs_P50 Gs_slope    Al2As
    ## 1 0.2200 -2.207094 89.41176       NA
    ## 2 0.2175 -1.871216 97.43590  631.000
    ## 3 0.2200 -2.114188 44.70588 1540.671

We can use function
[`modifySpParams()`](https://emf-creaf.github.io/medfate/reference/modifyParams.md)
to replace the values of parameters for the desired traits, leaving the
rest unaltered:

``` r

SpParamsFB <- modifySpParams(SpParamsMED, fb$customParams)
SpParamsFB
```

    ##                    Name SpIndex IFNcodes        AcceptedName
    ## 143 Phillyrea latifolia     142        8 Phillyrea latifolia
    ## 149    Pinus halepensis     148       24    Pinus halepensis
    ## 169        Quercus ilex     168   45/245        Quercus ilex
    ##                 Species     Genus   Family    Order      Group GrowthForm
    ## 143 Phillyrea latifolia Phillyrea Oleaceae Lamiales Angiosperm       Tree
    ## 149    Pinus halepensis     Pinus Pinaceae  Pinales Gymnosperm       Tree
    ## 169        Quercus ilex   Quercus Fagaceae  Fagales Angiosperm       Tree
    ##         LifeForm LeafShape LeafSize      PhenologyType DispersalType Hmed Hmax
    ## 143 Phanerophyte     Broad   Medium oneflush-evergreen    vertebrate   80  872
    ## 149 Phanerophyte    Needle    Small oneflush-evergreen          wind  850 1870
    ## 169 Phanerophyte     Broad   Medium oneflush-evergreen    vertebrate   48  415
    ##     Dmax Z50     Z95 fHDmin fHDmax       a_ash    b_ash     a_bsh     b_bsh
    ## 143   NA  NA 1396.48     35     90 0.008047873 2.865021 0.2502494 0.7300031
    ## 149   NA  NA 7500.00     80    160          NA       NA        NA        NA
    ## 169  200  NA 7000.00     NA     NA 1.857486249 1.885548 0.5238830 0.7337293
    ##        a_btsh    b_btsh cr BTsh      a_fbt    b_fbt        c_fbt a_cr b_1cr
    ## 143 0.3105815 0.7509837 NA   NA         NA       NA           NA   NA    NA
    ## 149        NA        NA NA   NA 0.05684887 1.521218 -0.024946517   NA    NA
    ## 169 0.7327147 0.7375770 NA   NA 0.06310630 1.545032 -0.005288476   NA    NA
    ##     b_2cr b_3cr c_1cr c_2cr      a_cw  b_cw      a_bt      b_bt LeafDuration
    ## 143    NA    NA    NA    NA        NA    NA        NA        NA     2.500000
    ## 149    NA    NA    NA    NA 0.6415296 0.731 0.5535741 1.1848613     2.083333
    ## 169    NA    NA    NA    NA        NA    NA 0.5622245 0.9626839     2.000000
    ##     t0gdd Sgdd Tbgdd Ssen Phsen Tbsen xsen ysen      SLA LeafDensity
    ## 143    NA   NA    NA   NA    NA    NA   NA   NA 5.400000     0.56000
    ## 149    NA   NA    NA   NA    NA    NA   NA   NA 5.140523     0.28812
    ## 169    NA   NA    NA   NA    NA    NA   NA   NA 6.340000     0.66500
    ##     WoodDensity FineRootDensity conduit2sapwood     r635   pDead    Al2As Ar2Al
    ## 143        0.70              NA              NA 1.917579      NA 1700.000    NA
    ## 149        0.60              NA              NA 1.964226 0.00050  631.000    NA
    ## 169        0.89              NA              NA 1.805872 0.00026 1540.671    NA
    ##     LeafWidth      SRL RLD   maxFMC   minFMC  Ptlp LeafPI0 LeafEPS LeafAF
    ## 143 0.4000000 902.1171  NA 106.2825 41.20000    NA   -2.13   12.38    0.5
    ## 149 0.1384772       NA  NA 122.1644 53.74550 -2.20   -1.50    5.31    0.6
    ## 169 1.7674359 735.7025  NA 109.2304 58.66087 -3.09   -2.50   15.00    0.4
    ##     StemPI0 StemEPS StemAF  SAV HeatContent LeafLigninPercent WoodLigninPercent
    ## 143   -2.13   12.38    0.4   NA    17146.03                NA                NA
    ## 149   -1.65    5.00    0.4 6050    22150.00                NA                NA
    ## 169   -2.50   15.00    0.4 4050    19300.00                NA                NA
    ##     FineRootLigninPercent LeafAngle LeafAngleSD ClumpingIndex gammaSWR alphaSWR
    ## 143                    NA  35.14508    18.87379            NA       NA       NA
    ## 149                    NA        NA          NA            NA       NA       NA
    ## 169                    NA  36.00819    20.15953            NA       NA       NA
    ##     kPAR  g  Tmax_LAI   Tmax_LAIsq Psi_Extract Exp_Extract      WUE   WUE_par
    ## 143   NA NA        NA           NA          NA          NA       NA        NA
    ## 149   NA NA 0.1869849 -0.008372458  -0.9218219    1.504542 8.525550 0.5239136
    ## 169   NA NA 0.1251027 -0.005601615  -1.9726871    1.149052 8.968208 0.1412266
    ##         WUE_co2    WUE_vpd Gswmin Gswmax Gsw_Toptim_Jarvis Gsw_Tsens_Jarvis
    ## 143          NA         NA  0.002 0.2200                NA               NA
    ## 149 0.002586327 -0.2647169  0.001 0.2175                NA               NA
    ## 169 0.002413091 -0.5664879  0.002 0.2200                NA               NA
    ##     Gsw_AC_slope_Baldocchi Gsw_P50_Baldocchi Gsw_slope_Baldocchi VCleaf_kmax
    ## 143                     NA                NA                  NA        3.00
    ## 149                     NA                NA                  NA        4.00
    ## 169                     NA                NA                  NA        2.63
    ##     VCleaf_P12 VCleaf_P50 VCleaf_P88 VCleaf_slope Kmax_stemxylem VCstem_P12
    ## 143         NA         NA         NA           NA      0.4083769  -1.971750
    ## 149 -1.9793246  -2.303772  -2.547056    133.86620      0.1500000  -3.707158
    ## 169 -0.5559123  -1.964085  -4.525766     19.14428      0.4000000  -4.739642
    ##     VCstem_P50 VCstem_P88 VCstem_slope Kmax_rootxylem VCroot_P12 VCroot_P50
    ## 143      -6.50 -11.028250           11             NA -3.1224800  -5.300000
    ## 149      -4.79  -5.872842           46             NA -1.0000000  -1.741565
    ## 169      -6.40  -8.060358           30             NA -0.4766469  -1.684034
    ##     VCroot_P88 VCroot_slope  Vmax298  Jmax298   Nleaf Nsapwood Nfineroot  WoodC
    ## 143  -7.477520     17.45105       NA       NA 15.7000     2.53        NA     NA
    ## 149  -2.301482    103.96607 72.19617 124.1687 11.7940       NA  12.90000 0.4981
    ## 169  -3.880455     22.32794 68.51600 118.7863 13.9021     5.63   2.45857 0.4751
    ##     RERleaf RERsapwood RERfineroot CCleaf CCsapwood CCfineroot RGRleafmax
    ## 143      NA         NA          NA   1.63        NA         NA         NA
    ## 149      NA         NA          NA     NA        NA         NA         NA
    ## 169      NA         NA          NA   1.43        NA         NA         NA
    ##     RGRsapwoodmax RGRcambiummax RGRfinerootmax RGRbud SRsapwood SRfineroot RSSG
    ## 143            NA  0.0006653797             NA     NA        NA         NA   NA
    ## 149            NA  0.0026280949             NA     NA        NA         NA 1.35
    ## 169            NA            NA             NA     NA        NA         NA 3.02
    ##     MortalityBaselineRate SurvivalModelStep SurvivalB0 SurvivalB1
    ## 143           0.001622378                NA         NA         NA
    ## 149           0.005000000                10   7.311515 -0.6532989
    ## 169           0.001000000                10   7.484348 -0.5420550
    ##     SeedProductionHeight SeedProductionDiameter  SeedMass SeedLongevity
    ## 143                   NA                     NA   31.5244            NA
    ## 149                   NA                     NA   19.6100            10
    ## 169                   NA               10.64702 3225.8100            NA
    ##     DispersalDistance DispersalShape   ProbRecr MinTempRecr MinMoistureRecr
    ## 143                NA             NA 0.04459023   -2.570181      0.05070956
    ## 149                NA             NA 0.02473379    1.083300      0.10154153
    ## 169                NA             NA 0.03005748   -3.744526      0.09657161
    ##     MinFPARRecr RecrAge RecrTreeDBH RecrTreeHeight RecrShrubHeight
    ## 143   0.4943654      NA          NA       52.54367              NA
    ## 149   4.5625766      NA          NA       56.93647              NA
    ## 169   0.1307250      NA          NA       47.23629              NA
    ##     RecrTreeDensity RecrShrubCover RespFire RespDist RespClip
    ## 143              NA             NA      0.9     0.95     0.96
    ## 149              NA             NA       NA       NA       NA
    ## 169              NA             NA      0.9     0.95     0.96
    ##     IngrowthTreeDensity IngrowthTreeDBH
    ## 143            235.1347              NA
    ## 149            246.2793              NA
    ## 169            352.2668              NA

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
    ## 1 Phillyrea latifolia  2.587859  323.0000 1248 390 1470 0.2581029
    ## 2    Pinus halepensis 26.759914 1195.7667  256 300 1200 1.0035486
    ## 3        Quercus ilex  6.220031  495.5532 3104 500 2287 1.4383485

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
    ## 1 Phillyrea latifolia  2.587859  323.0000 1248 390 1470 0.2581029
    ## 2    Pinus halepensis 26.759914 1195.7667  256 300 1200 1.0035486
    ## 3        Quercus ilex  6.220031  495.5532 3104 500 2287 1.4383485
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

    ##         SP    N       DBH Cover         H        CR  LAI_live LAI_expanded
    ## T1_142 142 1248  2.587859    NA  323.0000 0.5510653 0.2581029    0.2581029
    ## T2_148 148  256 26.759914    NA 1195.7667 0.6001765 1.0035486    1.0035486
    ## T3_168 168 3104  6.220031    NA  495.5532 0.6059123 1.4383485    1.4383485
    ##        LAI_dead LAI_nocomp LAI_mistletoe Age ObsID
    ## T1_142        0  0.2581029             0  NA  <NA>
    ## T2_148        0  1.0035486             0  NA  <NA>
    ## T3_168        0  1.4383485             0  NA  <NA>

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

![](FontBlanche_files/figure-html/unnamed-chunk-14-1.png)

Regarding **belowground** information, the usuer should supply a matrix
describing for each plant cohort, the proportion of fine roots in each
soil layer. As before, we use internal function
[`forest2belowground()`](https://emf-creaf.github.io/medfate/reference/forest2aboveground.md)
on the object `fb_forest` to show how should the data look like:

``` r

fb_below <- forest2belowground(fb_forest, fb_soil, SpParamsFB)
fb_below
```

    ##                1         2          3          4
    ## T1_142 0.3602157 0.5332967 0.08477533 0.02171222
    ## T2_148 0.5016024 0.4291685 0.05479894 0.01443019
    ## T3_168 0.2752236 0.5286425 0.14537757 0.05075634

In our case, these proportions were implicitly specified in parameters
`Z50` and `Z95`. In fact, these values describe a continuous
distribution of fine roots along depth, which can be displayed using
function
[`vprofile_rootDistribution()`](https://emf-creaf.github.io/medfate/reference/vprofile_leafAreaDensity.md):

``` r

vprofile_rootDistribution(fb_forest, SpParamsFB, bySpecies = T)
```

![](FontBlanche_files/figure-html/unnamed-chunk-16-1.png)

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

fb_control <- defaultControl()
fb_control$transpirationMode <- "Sperry"
fb_control$subdailyResults <- TRUE
fb_control$stemCavitationRecovery <- "rate"
fb_control$leafCavitationRecovery <- "total"
fb_control$fracRootResistance <- 0.4
```

Where the following changes are set to control parameters:

1.  Transpiration is set `transpirationMode = "Sperry"`, which implies a
    greater complexity of plant hydraulics and energy balance
    calculations.
2.  Soil water retention curves are calculated using Van Genuchten’s
    equations.
3.  Subdaily results generated by the model are kept.
4.  Coarse root resistance is assumed to be 40% of total plant
    resistance

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

    ## Initial plant water content (mm): 31.3483
    ## Initial soil water content (mm): 213.886
    ## Initial snowpack content (mm): 0
    ## Performing daily simulations
    ## 
    ##  [Year 2014]:............
    ## 
    ## Final plant water content (mm): 31.3067
    ## Final soil water content (mm): 234.755
    ## Final snowpack content (mm): 0
    ## Change in plant water content (mm): -0.0416382
    ## Plant water balance result (mm): 2.68215e-16
    ## Change in soil water content (mm): 20.8696
    ## Soil water balance result (mm): 20.8696
    ## Change in snowpack water content (mm): 0
    ## Snowpack water balance result (mm): 0
    ## Water balance components:
    ##   Precipitation (mm) 1066 Rain (mm) 1066 Snow (mm) 0
    ##   Interception (mm) 141 Net rainfall (mm) 925
    ##   Infiltration (mm) 832 Infiltration excess (mm) 93 Saturation excess (mm) 268 Capillarity rise (mm) 0
    ##   Soil evaporation (mm) 23  Herbaceous transpiration (mm) 0  Woody plant transpiration (mm) 327  Mistletoe transpiration (mm) 0
    ##   Plant extraction from soil (mm) 327  Plant water balance (mm) 0 Hydraulic redistribution (mm) 36
    ##   Runoff (mm) 361 Deep drainage (mm) 193

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

    ##                  PET Precipitation     Rain Snow   NetRain Snowmelt
    ## 2014-01-01 0.6209989      0.000000 0.000000    0 0.0000000        0
    ## 2014-01-02 0.5671238      0.000000 0.000000    0 0.0000000        0
    ## 2014-01-03 0.5418115      0.000000 0.000000    0 0.0000000        0
    ## 2014-01-04 0.6072565      2.590674 2.590674    0 0.7213133        0
    ## 2014-01-05 2.0447148      0.000000 0.000000    0 0.0000000        0
    ## 2014-01-06 0.9330456      0.000000 0.000000    0 0.0000000        0
    ##            Infiltration InfiltrationExcess SaturationExcess Runoff DeepDrainage
    ## 2014-01-01    0.0000000                  0                0      0    0.0000000
    ## 2014-01-02    0.0000000                  0                0      0    0.0000000
    ## 2014-01-03    0.0000000                  0                0      0    0.0000000
    ## 2014-01-04    0.7213133                  0                0      0    0.1890744
    ## 2014-01-05    0.0000000                  0                0      0    0.0000000
    ## 2014-01-06    0.0000000                  0                0      0    0.0000000
    ##            CapillarityRise Evapotranspiration Interception SoilEvaporation
    ## 2014-01-01               0          0.2361114      0.00000       0.2145403
    ## 2014-01-02               0          0.1959278      0.00000       0.1959278
    ## 2014-01-03               0          0.1966164      0.00000       0.1871830
    ## 2014-01-04               0          2.0496595      1.86936       0.1802993
    ## 2014-01-05               0          0.9738787      0.00000       0.2947672
    ## 2014-01-06               0          0.6695831      0.00000       0.1403825
    ##            HerbTranspiration PlantExtraction Transpiration
    ## 2014-01-01                 0    2.157104e-02   0.021571043
    ## 2014-01-02                 0   -3.523657e-19   0.000000000
    ## 2014-01-03                 0    9.433443e-03   0.009433443
    ## 2014-01-04                 0   -2.195509e-18   0.000000000
    ## 2014-01-05                 0    6.791115e-01   0.679111496
    ## 2014-01-06                 0    5.292006e-01   0.529200608
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2014-01-01                      0             0.000000000
    ## 2014-01-02                      0             0.002564064
    ## 2014-01-03                      0             0.003646811
    ## 2014-01-04                      0             0.004565272
    ## 2014-01-05                      0             0.000000000
    ## 2014-01-06                      0             0.001544657

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

    ##                 dates       SWC SWC.err       ETR E_T2_148 E_T2_148_err
    ## 2014-01-01 2014-01-01 0.5813407      NA 0.2259528       NA           NA
    ## 2014-01-02 2014-01-02 0.6507478      NA 0.2337668       NA           NA
    ## 2014-01-03 2014-01-03 0.6224243      NA 0.5229000       NA           NA
    ## 2014-01-04 2014-01-04        NA      NA 0.1117191       NA           NA
    ## 2014-01-05 2014-01-05 0.6285134      NA 0.8132403       NA           NA
    ## 2014-01-06 2014-01-06 0.6035415      NA 0.6012234       NA           NA
    ##            E_T3_168 E_T3_168_err PD_T2_148 PD_T2_148_err PD_T3_168
    ## 2014-01-01       NA           NA        NA            NA        NA
    ## 2014-01-02       NA           NA        NA            NA        NA
    ## 2014-01-03       NA           NA        NA            NA        NA
    ## 2014-01-04       NA           NA        NA            NA        NA
    ## 2014-01-05       NA           NA        NA            NA        NA
    ## 2014-01-06       NA           NA        NA            NA        NA
    ##            PD_T3_168_err MD_T2_148 MD_T2_148_err MD_T3_168 MD_T3_168_err
    ## 2014-01-01            NA        NA            NA        NA            NA
    ## 2014-01-02            NA        NA            NA        NA            NA
    ## 2014-01-03            NA        NA            NA        NA            NA
    ## 2014-01-04            NA        NA            NA        NA            NA
    ## 2014-01-05            NA        NA            NA        NA            NA
    ## 2014-01-06            NA        NA            NA        NA            NA

### Stand evapotranspiration

Package `medfate` contains several functions to assist the evaluation of
model results. We can first compare the observed vs modelled total
evapotranspiration. We can plot the two time series:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "ETR", plotType="dynamics")+
  theme(legend.position = c(0.8,0.85))
```

![](FontBlanche_files/figure-html/unnamed-chunk-25-1.png)

It is easy to see that in rainy days the predicted evapotranspiration is
much higher than that of the observed data. We repeat the comparison but
excluding the intercepted water from modeled results:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "SE+TR", plotType="dynamics")+
  theme(legend.position = c(0.8,0.85))
```

![](FontBlanche_files/figure-html/unnamed-chunk-26-1.png)

The relationship can be shown in a scatter plot:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "SE+TR", plotType="scatter")
```

![](FontBlanche_files/figure-html/unnamed-chunk-27-1.png) Where we see a
reasonably good relationship, but the model tends to underestimate total
evapotranspiration during seasons with low evaporative demand. Function
[`evaluation_stats()`](https://emf-creaf.github.io/medfate/reference/evaluation.md)
allows us to generate evaluation statistics:

``` r

evaluation_stats(fb_SWB, fb_observed, type = "SE+TR")
```

    ##            n         Bias     Bias.rel          MAE      MAE.rel            r 
    ## 365.00000000  -0.37458163 -28.11376353   0.47417785  35.58883461   0.68150269 
    ##          NSE      NSE.abs 
    ##  -0.00628706   0.04767776

### Soil moisture

We can compare observed vs modelled soil moisture content in a similar
way as we did for total evapotranspiration:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "SWC", plotType="dynamics")
```

![](FontBlanche_files/figure-html/unnamed-chunk-29-1.png)

As before, we can generate a scatter plot:

``` r

evaluation_plot(fb_SWB, fb_observed, type = "SWC", plotType="scatter")
```

![](FontBlanche_files/figure-html/unnamed-chunk-30-1.png)

or examine evaluation statistics:

``` r

evaluation_stats(fb_SWB, fb_observed, type = "SWC")
```

    ##            n         Bias     Bias.rel          MAE      MAE.rel            r 
    ## 364.00000000  -0.13329240 -31.07487385   0.13329240  31.07487385   0.95519253 
    ##          NSE      NSE.abs 
    ##  -0.10516414  -0.04428937

### Plant transpiration

The following plots display the observed and predicted transpiration
dynamics for *Pinus halepensis* and *Quercus ilex*:

``` r

g1<-evaluation_plot(fb_SWB, fb_observed, 
                            cohort = "T2_148",
                            type="E", plotType = "dynamics")+
  theme(legend.position = c(0.85,0.83))
g2<-evaluation_plot(fb_SWB, fb_observed, 
                            cohort = "T3_168",
                            type="E", plotType = "dynamics")+
  theme(legend.position = c(0.85,0.83))
plot_grid(g1, g2, ncol=1)
```

![](FontBlanche_files/figure-html/unnamed-chunk-32-1.png)

In general, the agreement is quite good, but the model seems to
overestimate the transpiration of *P. halepensis* in early summer and
after the first drought period. The transpiration of *Q. ilex* seems
also overestimated in spring and autumn. We can also inspect the
evaluation statistics for both species using:

``` r

evaluation_stats(fb_SWB, fb_observed, cohort = "T2_148", type="E")
```

    ##           n        Bias    Bias.rel         MAE     MAE.rel           r 
    ## 300.0000000   0.2419407 117.6344214   0.2693672 130.9695180   0.7432842 
    ##         NSE     NSE.abs 
    ##  -8.4643333  -1.7220938

``` r

evaluation_stats(fb_SWB, fb_observed, cohort = "T3_168", type="E")
```

    ##            n         Bias     Bias.rel          MAE      MAE.rel            r 
    ## 309.00000000   0.03209984  11.08961119   0.07788813  26.90821018   0.90667188 
    ##          NSE      NSE.abs 
    ##   0.62729587   0.46936452

### Leaf water potentials

Finally, we can compare observed with predicted water potentials. In
this case measurements are available for three dates, but they include
the standard deviation of several measurements.

``` r

g1<-evaluation_plot(fb_SWB, fb_observed, 
                            cohort = "T2_148",
                            type="WP", plotType = "dynamics")+
  theme(legend.position = c(0.85,0.23))
g2<-evaluation_plot(fb_SWB, fb_observed, 
                            cohort = "T3_168",
                            type="WP", plotType = "dynamics")+
  theme(legend.position = c(0.85,0.23))
plot_grid(g1, g2, ncol=1)
```

![](FontBlanche_files/figure-html/unnamed-chunk-34-1.png)

The model seems to underestimate water potentials (i.e. it predicts more
negative values than those observed) during the drought season.

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

![](FontBlanche_files/figure-html/unnamed-chunk-36-1.png) It is apparent
the climatic drought period between april and august 2014. This should
have an impact on soil moisture and plant stress.

If we are interested in forest hydrology, we can plot the amount of
water that the model predicts to leave the forest via surface runoff or
drainage to lower water compartments.

``` r

plot(fb_SWB, type = "Export")
```

![](FontBlanche_files/figure-html/unnamed-chunk-37-1.png) As expected,
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

![](FontBlanche_files/figure-html/unnamed-chunk-38-1.png)

### Soil moisture dynamics and hydraulic redistribution

It is also useful to plot the dynamics of soil state variables by layer,
such as the percentage of moisture in relation to field capacity:

``` r

plot(fb_SWB, type="SoilTheta")
```

![](FontBlanche_files/figure-html/unnamed-chunk-39-1.png) Note that the
model predicts soil drought to occur earlier in the season for the first
three layers (0-200 cm) whereas the bottom layer dries out much more
slowly. At this point is important to mention that the water balance
model incorporates. We can also display the dynamics of the
corresponding soil layer water potentials:

``` r

plot(fb_SWB, type="SoilPsi")
```

![](FontBlanche_files/figure-html/unnamed-chunk-40-1.png) or draw a
composite plot including absolute soil water volume:

``` r

g1<-plot(fb_SWB)+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = "none")
g2<-plot(fb_SWB, "SoilVol")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.08,0.65))
g3<-plot(fb_SWB, "SoilPsi")+scale_x_date(date_breaks = "1 month", date_labels = "%m")+theme(legend.position = c(0.08,0.5))
plot_grid(g1, g2,  g3, rel_heights = c(0.4,0.8,0.8), ncol=1)
```

![](FontBlanche_files/figure-html/unnamed-chunk-41-1.png)

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

![](FontBlanche_files/figure-html/unnamed-chunk-42-1.png)

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

![](FontBlanche_files/figure-html/unnamed-chunk-43-1.png)

### Plant transpiration

We can use function
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to display the
seasonal dynamics of cohort-level variables, such as plant transpiration
per leaf area:

``` r

par(mar=c(5,5,1,1))
plot(fb_SWB, type="TranspirationPerLeaf", bySpecies = T)
```

![](FontBlanche_files/figure-html/unnamed-chunk-44-1.png) Where we can
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

![](FontBlanche_files/figure-html/unnamed-chunk-45-1.png)

### Plant stress

In the model, reduction of (whole-plant) plant transpiration is what
used to define drought stress, which depends on the species identity:

``` r

plot(fb_SWB, type="PlantStress", bySpecies = T)
```

![](FontBlanche_files/figure-html/unnamed-chunk-46-1.png)

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

![](FontBlanche_files/figure-html/unnamed-chunk-47-1.png)

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

![](FontBlanche_files/figure-html/unnamed-chunk-48-1.png)

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

![](FontBlanche_files/figure-html/unnamed-chunk-49-1.png)

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

    ##                  PET Precipitation      Rain Snow     NetRain Snowmelt
    ## 2014-01-01  27.03414     205.04814 205.04814    0 182.5767907        0
    ## 2014-02-01  37.11592     181.09641 181.09641    0 155.2573002        0
    ## 2014-03-01  80.49737      44.61248  44.61248    0  39.8917051        0
    ## 2014-04-01 109.24874      15.00000  15.00000    0   7.2713589        0
    ## 2014-05-01 147.99639      21.60000  21.60000    0  16.3281633        0
    ## 2014-06-01 167.27898      33.60000  33.60000    0  25.8839490        0
    ## 2014-07-01 183.99299       0.60000   0.60000    0   0.1428946        0
    ## 2014-08-01 159.66330      60.40000  60.40000    0  52.8025568        0
    ## 2014-09-01 103.42793     137.60000 137.60000    0 125.8957242        0
    ## 2014-10-01  63.53896      50.60000  50.60000    0  41.9066889        0
    ## 2014-11-01  30.12083     222.60000 222.60000    0 198.0975096        0
    ## 2014-12-01  26.01617      93.00000  93.00000    0  78.7534201        0
    ##            Infiltration InfiltrationExcess SaturationExcess     Runoff
    ## 2014-01-01  182.5767907           0.000000        86.696614  86.696614
    ## 2014-02-01  155.2573002           0.000000       122.076919 122.076919
    ## 2014-03-01   39.8917051           0.000000         0.000000   0.000000
    ## 2014-04-01    7.2713589           0.000000         0.000000   0.000000
    ## 2014-05-01   16.3281633           0.000000         0.000000   0.000000
    ## 2014-06-01   25.8839490           0.000000         0.000000   0.000000
    ## 2014-07-01    0.1428946           0.000000         0.000000   0.000000
    ## 2014-08-01   43.8641244           8.938432         0.000000   8.938432
    ## 2014-09-01  101.0233304          24.872394         0.000000  24.872394
    ## 2014-10-01   32.3395984           9.567091         0.000000   9.567091
    ## 2014-11-01  148.2069996          49.890510         9.060547  58.951057
    ## 2014-12-01   78.7534201           0.000000        50.210110  50.210110
    ##            DeepDrainage CapillarityRise Evapotranspiration Interception
    ## 2014-01-01     27.43226               0           34.15932   22.4713498
    ## 2014-02-01     40.14785               0           40.34980   25.8391081
    ## 2014-03-01     39.78432               0           47.04269    4.7207713
    ## 2014-04-01      7.45117               0           52.09614    7.7286411
    ## 2014-05-01      0.00000               0           47.72911    5.2718367
    ## 2014-06-01      0.00000               0           51.95191    7.7160510
    ## 2014-07-01      0.00000               0           24.24042    0.4571054
    ## 2014-08-01      0.00000               0           42.58144    7.5974432
    ## 2014-09-01      0.00000               0           39.68645   11.7042748
    ## 2014-10-01      0.00000               0           43.00754    8.6933111
    ## 2014-11-01     33.76513               0           39.75339   24.5024904
    ## 2014-12-01     44.44940               0           27.94647   14.2465799
    ##            SoilEvaporation HerbTranspiration PlantExtraction Transpiration
    ## 2014-01-01       4.2151687                 0         7.47280       7.47280
    ## 2014-02-01       2.9061658                 0        11.60453      11.60453
    ## 2014-03-01       3.0583127                 0        39.26361      39.26361
    ## 2014-04-01       0.4436067                 0        43.92389      43.92389
    ## 2014-05-01       0.2200408                 0        42.23723      42.23723
    ## 2014-06-01       0.1508599                 0        44.08500      44.08500
    ## 2014-07-01       0.1172265                 0        23.66608      23.66608
    ## 2014-08-01       0.2549165                 0        34.72908      34.72908
    ## 2014-09-01       0.8386994                 0        27.14347      27.14347
    ## 2014-10-01       2.8421025                 0        31.47213      31.47213
    ## 2014-11-01       4.2636858                 0        10.98722      10.98722
    ## 2014-12-01       3.4948016                 0        10.20509      10.20509
    ##            MistletoeTranspiration HydraulicRedistribution
    ## 2014-01-01                      0               0.3900485
    ## 2014-02-01                      0               0.4678481
    ## 2014-03-01                      0               0.8975168
    ## 2014-04-01                      0               1.4812743
    ## 2014-05-01                      0               3.9017283
    ## 2014-06-01                      0               4.1778300
    ## 2014-07-01                      0               2.2033719
    ## 2014-08-01                      0               7.5679965
    ## 2014-09-01                      0               9.7632514
    ## 2014-10-01                      0               3.3999496
    ## 2014-11-01                      0               0.6338105
    ## 2014-12-01                      0               0.7962439

Parameter `output` is used to indicate the element of the `spwb` object
for which we desire summaries. Similarly, it is possible to calculate
the average stress of the three tree species by months:

``` r

summary(fb_SWB, freq="months",FUN=mean, output="PlantStress", bySpecies = TRUE)
```

    ##            Phillyrea latifolia Pinus halepensis Quercus ilex
    ## 2014-01-01        0.0002005823     9.213874e-05  0.007289001
    ## 2014-02-01        0.0004241553     1.801736e-04  0.012831286
    ## 2014-03-01        0.0021615694     1.558130e-03  0.037302946
    ## 2014-04-01        0.0055229219     1.916886e-02  0.059477525
    ## 2014-05-01        0.0584843596     2.475189e-01  0.098947305
    ## 2014-06-01        0.2596275575     4.129932e-01  0.206159134
    ## 2014-07-01        0.5508210211     5.917795e-01  0.535269874
    ## 2014-08-01        0.3855427230     3.661730e-01  0.487188044
    ## 2014-09-01        0.2530723414     2.897205e-01  0.328572652
    ## 2014-10-01        0.0135198518     1.955654e-03  0.034691130
    ## 2014-11-01        0.0083424116     2.206970e-04  0.011217083
    ## 2014-12-01        0.0045925906     9.153584e-05  0.007592380

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
