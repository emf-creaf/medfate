# Preparing model inputs

## About this article

A companion article [*Understanding model
inputs*](https://emf-creaf.github.io/medfate/articles/intro/UnderstandingInputs.html)
explained the vegetation, soil and weather structures needed to run the
simulation models included in **medfate**. Preparing inputs for
simulations with **medfate** is not straightforward, because it requires
obtaining and reshaping data for vegetation, soil and weather.
Therefore, this article illustrates some common issues that arise in the
process of preparing inputs, so that the user is aware of them when
processing his/her own data.

We begin by loading packages **medfate** and **meteoland**:

``` r
library(medfate)
#> Package 'medfate' [ver. 4.8.5]
library(meteoland)
#> Package 'meteoland' [ver. 2.2.4]
```

## Building/manipulating forest objects

In this section we show how to build and manipulate objects of class
`forest`, for their use in package medfate, starting from a table
containing forest inventory data.

### Poblet tree data set

Package **medfate** includes a small dataset of tree data, corresponding
to a dense holm oak forest in Poblet (Catalonia, Spain). As a result of
the abandonment of previous coppicing exploitation, there is a high
density of stems per individual.

We begin by loading the tree data from Poblet:

``` r
data("poblet_trees")
```

and we inspect its content, for example using:

``` r
summary(poblet_trees)
#>   Plot.Code            Indv.Ref       Species           Diameter.cm   
#>  Length:717         Min.   :  1.0   Length:717         Min.   : 7.50  
#>  Class :character   1st Qu.: 45.0   Class :character   1st Qu.: 9.10  
#>  Mode  :character   Median : 97.0   Mode  :character   Median :11.10  
#>                     Mean   :103.4                      Mean   :11.62  
#>                     3rd Qu.:156.0                      3rd Qu.:13.40  
#>                     Max.   :261.0                      Max.   :26.00
```

The data frame includes tree data corresponding to three forest
inventories:

``` r
table(poblet_trees$Plot.Code)
#> 
#>     POBL_CTL POBL_THI_AFT POBL_THI_BEF 
#>          267          189          261
```

`POBL_CTL` corresponds to an oak forest where no treatment was done
(control), whereas `POBL_THI_BEF` and `POBL_THI_AFT` are two forest
inventories conducted on the same forest plot, before and after a
thinning intervention to reduce the number of stems.

### Mapping trees from the control forest

We initialize an empty forest object using function
[`emptyforest()`](https://emf-creaf.github.io/medfate/reference/emptyforest.md)
from package **medfate**:

``` r
pobl_ctl <- emptyforest()
pobl_ctl
#> $treeData
#> [1] Species DBH     Height  N       Z50     Z95    
#> <0 rows> (or 0-length row.names)
#> 
#> $shrubData
#> [1] Species Height  Cover   Z50     Z95    
#> <0 rows> (or 0-length row.names)
#> 
#> $herbCover
#> [1] NA
#> 
#> $herbHeight
#> [1] NA
#> 
#> $seedBank
#> [1] Species Percent
#> <0 rows> (or 0-length row.names)
#> 
#> attr(,"class")
#> [1] "forest" "list"
```

Now we will fill in data for element `treeData` in the `forest` object.
For that, we need to define a mapping from column names in
`poblet_trees` to variables in `treeData`. The mapping can be defined
using a **named string vector**, i.e. a vector where element names are
variable names in `treeData` and vector elements are strings of the
variable names in `poblet_trees`:

``` r
mapping <- c("Species.name" = "Species", "DBH" = "Diameter.cm")
```

We can now replace the empty `treeData` in `pobl_ctl` using functions
[`subset()`](https://rdrr.io/r/base/subset.html) and
[`forest_mapTreeTable()`](https://emf-creaf.github.io/medfate/reference/forest_mapWoodyTables.md):

``` r
pobl_ctl$treeData <- forest_mapTreeTable(subset(poblet_trees, Plot.Code=="POBL_CTL"), 
                                         mapping_x = mapping, SpParams = SpParamsMED)
```

We can inspect the result using:

``` r
summary(pobl_ctl$treeData)
#>    Species                N      Height             DBH          Z50         
#>  Length:267         Min.   :1   Mode:logical   Min.   : 7.50   Mode:logical  
#>  Class :character   1st Qu.:1   NA's:267       1st Qu.: 9.00   NA's:267      
#>  Mode  :character   Median :1                  Median :10.70                 
#>                     Mean   :1                  Mean   :11.53                 
#>                     3rd Qu.:1                  3rd Qu.:13.30                 
#>                     Max.   :1                  Max.   :26.00                 
#>    Z95         
#>  Mode:logical  
#>  NA's:267      
#>                
#>                
#>                
#> 
```

Some data are missing, but we will not worry about it now. One way to
evaluate if the tree data is correctly specified is to display a summary
of the `forest` object using the `summary` function defined in
**medfate** for this object class:

``` r
summary(pobl_ctl, SpParamsMED)
#> Tree BA (m2/ha): 3.0179815  adult trees: 3.0179815  saplings: 0 
#> Density (ind/ha) adult trees: 267  saplings: 0  shrubs (estimated): 0 
#> Cover (%) adult trees: 42.1205627  saplings: 0  shrubs: 0  herbs: 0 
#> LAI (m2/m2) total: 0.544959  adult trees: 0.544959  saplings: 0  shrubs: 0  herbs: 0 
#> Fuel loading (kg/m2) total: 0.1421746  adult trees: 0.1421746  saplings: 0  shrubs: 0  herbs: 0 
#> PAR ground (%): NA  SWR ground (%): NA
```

The values of stand density and stand basal area are too low for such a
dense forest, which indicates that something needs to be corrected. At
this point, it is important to remember that `forest` objects need the
density of trees to be specified as *stems per hectare*. We conducted
our tree data mapping without indicating the area of the sampled plot.
We are told that forest stand sampling was done using a circular plot
whose radius was 15 m. We can calculate the sampled area using:

``` r
sampled_area <- pi*15^2
```

and use this information to map the tree data again, where we specify
parameter `plot_size_x`:

``` r
pobl_ctl$treeData <- forest_mapTreeTable(subset(poblet_trees, Plot.Code=="POBL_CTL"),
                                         mapping_x = mapping, SpParams = SpParamsMED, 
                                         plot_size_x = sampled_area)
```

We run again the summary:

``` r
summary(pobl_ctl, SpParamsMED)
#> Tree BA (m2/ha): 42.6957047  adult trees: 42.6957047  saplings: 0 
#> Density (ind/ha) adult trees: 3777.277316  saplings: 0  shrubs (estimated): 0 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 0  herbs: 0 
#> LAI (m2/m2) total: 5.6770407  adult trees: 5.6770407  saplings: 0  shrubs: 0  herbs: 0 
#> Fuel loading (kg/m2) total: 1.493419  adult trees: 1.493419  saplings: 0  shrubs: 0  herbs: 0 
#> PAR ground (%): NA  SWR ground (%): NA
```

which results in a much higher basal area and density, as should be
expected for a dense oak forest resulting from an abandoned old coppice.

Another issue that we see is the percentage of PAR and SWR that reaches
the ground, which have missing values. This indicates that medfate
cannot calculate the light extinction profile, in our case because tree
heights are missing. Thus, we should somehow estimate tree heights, for
example using an allometric relationship:

``` r
poblet_trees$Height.cm <- 100 * 1.806*poblet_trees$Diameter.cm^0.518
summary(poblet_trees$Height.cm)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   512.9   566.9   628.3   638.0   692.7   976.5
```

So trees are between 5 and 10 m height. Once tree heights are defined,
we can include them in our mapping:

``` r
mapping = c("Species.name" = "Species", "DBH" = "Diameter.cm", "Height" = "Height.cm")
```

and rerun the tree data mapping:

``` r
pobl_ctl$treeData <- forest_mapTreeTable(subset(poblet_trees, Plot.Code=="POBL_CTL"),
                                         mapping_x = mapping, SpParams = SpParamsMED, 
                                         plot_size_x = sampled_area)
```

Now the summary of the control forest stand looks like:

``` r
summary(pobl_ctl, SpParamsMED)
#> Tree BA (m2/ha): 42.6957047  adult trees: 42.6957047  saplings: 0 
#> Density (ind/ha) adult trees: 3777.277316  saplings: 0  shrubs (estimated): 0 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 0  herbs: 0 
#> LAI (m2/m2) total: 5.6770407  adult trees: 5.6770407  saplings: 0  shrubs: 0  herbs: 0 
#> Fuel loading (kg/m2) total: 1.493419  adult trees: 1.493419  saplings: 0  shrubs: 0  herbs: 0 
#> PAR ground (%): 4.4052535  SWR ground (%): 9.8976935
```

The fraction of PAR/SWR reaching the ground is low, as would be expected
for a dense forest.

### Mapping trees from the managed forest

Here we can repeat our mapping for the managed forest plot, which has
two codes corresponding to before and after the thinning intervention.
Let us first address the pre-thinning state:

``` r
pobl_thi_bef  <- emptyforest()
pobl_thi_bef$treeData <- forest_mapTreeTable(subset(poblet_trees, Plot.Code=="POBL_THI_BEF"),
                                             mapping_x = mapping, SpParams = SpParamsMED, 
                                             plot_size_x = sampled_area)
#> Warning in forest_mapTreeTable(subset(poblet_trees, Plot.Code ==
#> "POBL_THI_BEF"), : Taxon names that were not matched: Quercus humilis.
```

A warning is raised that not all species names could be parsed. In this
case, the reason is that the name used for the downy oak (*Quercus
humilis*) is a synonym and needs to be replaced by its accepted name
(*Quercus pubescens*), which we can do:

``` r
poblet_trees$Species[poblet_trees$Species=="Quercus humilis"] <- "Quercus pubescens"
```

Now we repeat our mapping:

``` r
pobl_thi_bef$treeData <- forest_mapTreeTable(subset(poblet_trees, Plot.Code=="POBL_THI_BEF"),
                                             mapping_x = mapping, SpParams = SpParamsMED, 
                                             plot_size_x = sampled_area)

summary(pobl_thi_bef, SpParamsMED)
#> Tree BA (m2/ha): 40.9224267  adult trees: 40.9224267  saplings: 0 
#> Density (ind/ha) adult trees: 3692.3946797  saplings: 0  shrubs (estimated): 0 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 0  herbs: 0 
#> LAI (m2/m2) total: 5.5833511  adult trees: 5.5833511  saplings: 0  shrubs: 0  herbs: 0 
#> Fuel loading (kg/m2) total: 1.4629714  adult trees: 1.4629714  saplings: 0  shrubs: 0  herbs: 0 
#> PAR ground (%): 4.6382035  SWR ground (%): 10.2827898
```

Like the control plot, these statistics indicate a dense oak forest. We
can repeat the same operations with the forest plot after the thinning
intervention:

``` r
pobl_thi_aft = emptyforest()
pobl_thi_aft$treeData <- forest_mapTreeTable(subset(poblet_trees, Plot.Code=="POBL_THI_AFT"),
                                             mapping_x = mapping, SpParams = SpParamsMED, 
                                             plot_size_x = sampled_area)
summary(pobl_thi_aft, SpParamsMED)
#> Tree BA (m2/ha): 31.6162035  adult trees: 31.6162035  saplings: 0 
#> Density (ind/ha) adult trees: 2673.8030439  saplings: 0  shrubs (estimated): 0 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 0  herbs: 0 
#> LAI (m2/m2) total: 4.5328748  adult trees: 4.5328748  saplings: 0  shrubs: 0  herbs: 0 
#> Fuel loading (kg/m2) total: 1.1915321  adult trees: 1.1915321  saplings: 0  shrubs: 0  herbs: 0 
#> PAR ground (%): 8.2654902  SWR ground (%): 15.7752682
```

Note the decrease in tree density and basal area, and the increase in
light reaching the ground, despite the estimated leaf area index is
still high.

### Reducing the number of woody cohorts

So far we have considered that each tree record should correspond to a
woody cohort. We can check the number of tree cohorts in each `forest`
structure using:

``` r
nrow(pobl_ctl$treeData)
#> [1] 267
nrow(pobl_thi_bef$treeData)
#> [1] 261
nrow(pobl_thi_aft$treeData)
#> [1] 189
```

This large amount of cohorts can slow done simulations considerably.
Hence, it is advisable to lump them into coarser woody cohorts. One way
of doing this is via function
[`forest_mergeTrees()`](https://emf-creaf.github.io/medfate/reference/forest_simplification.md)
from package **medfate**:

``` r
pobl_ctl <- forest_mergeTrees(pobl_ctl)
pobl_thi_bef <- forest_mergeTrees(pobl_thi_bef)
pobl_thi_aft <- forest_mergeTrees(pobl_thi_aft)
```

By default, the function will pool tree cohorts of the same species and
diameter class (defined every 5 cm). We can check the new number of tree
cohorts using again:

``` r
nrow(pobl_ctl$treeData)
#> [1] 9
nrow(pobl_thi_bef$treeData)
#> [1] 11
nrow(pobl_thi_aft$treeData)
#> [1] 8
```

We can check whether stand properties were altered using the
[`summary()`](https://rdrr.io/r/base/summary.html) function:

``` r
summary(pobl_thi_aft, SpParamsMED)
#> Tree BA (m2/ha): 31.6162035  adult trees: 31.6162035  saplings: 0 
#> Density (ind/ha) adult trees: 2673.8030439  saplings: 0  shrubs (estimated): 0 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 0  herbs: 0 
#> LAI (m2/m2) total: 4.0969956  adult trees: 4.0969956  saplings: 0  shrubs: 0  herbs: 0 
#> Fuel loading (kg/m2) total: 1.0724731  adult trees: 1.0724731  saplings: 0  shrubs: 0  herbs: 0 
#> PAR ground (%): 10.5046983  SWR ground (%): 18.8407832
```

Function
[`forest_mergeTrees()`](https://emf-creaf.github.io/medfate/reference/forest_simplification.md)
will preserve the stand density and basal area that the stand
description had before merging cohorts. Other properties like leaf area
index may be slightly modified.

In general, it is advisable to reduce the number of woody cohorts before
running simulation models in **medfate**.

## Retrieving SoilGrids data

Because soil properties vary strongly at fine spatial scales, ideally
soil physical attributes should be measured on samples taken at the
forest stand to be simulated. For those users lacking such data, soil
properties modelled at larger scales are available via SoilGrids.org.

Retrieval of soil properties from SoilGrids can be done using function
`add_soilgrids()` from package **medfateland**. Assuming we know the
plot coordinates, we first create an object `sf` (see package **sf**):

``` r
sf_pt <- sf::st_sfc(sf::st_point(c(1.0219, 41.3443)), crs = 4326)
```

With this function we will obtain, for each location, a data frame of
soil properties:

``` r
pobl_soil_props
#>   widths     clay     sand       om       bd  rfc
#> 1    300 26.43333 31.06667 4.133333 1.166667 18.0
#> 2    700 30.40000 29.75000 0.900000 1.440000 19.2
#> 3   1000 31.60000 29.60000 0.610000 1.500000 20.9
```

This data frame is a physical description of the soil. Initialization of
additional parameters and state variables is done using function
[`soil()`](https://emf-creaf.github.io/medfate/reference/soil.md):

``` r
pobl_soil <- soil(pobl_soil_props)
```

We can inspect the soil definition using:

``` r
summary(pobl_soil)
#> Soil depth (mm): 2000 
#> 
#> Layer  1  [ 0  to  300 mm ] 
#>     clay (%): 26 silt (%): 38 sand (%): 31 organic matter (%): 4 [ Loam ]
#>     Rock fragment content (%): 18 Macroporosity (%): 22 
#>     Theta WP (%): 18 Theta FC (%): 34 Theta SAT (%): 51 Theta current (%) 34 
#>     Vol. WP (mm): 44 Vol. FC (mm): 83 Vol. SAT (mm): 125 Vol. current (mm): 83 
#>     Temperature (Celsius): NA 
#> 
#> Layer  2  [ 300  to  1000 mm ] 
#>     clay (%): 30 silt (%): 39 sand (%): 30 organic matter (%): 1 [ Clay loam ]
#>     Rock fragment content (%): 19 Macroporosity (%): 9 
#>     Theta WP (%): 19 Theta FC (%): 33 Theta SAT (%): 44 Theta current (%) 33 
#>     Vol. WP (mm): 106 Vol. FC (mm): 186 Vol. SAT (mm): 250 Vol. current (mm): 186 
#>     Temperature (Celsius): NA 
#> 
#> Layer  3  [ 1000  to  2000 mm ] 
#>     clay (%): 32 silt (%): 38 sand (%): 30 organic matter (%): 1 [ Clay loam ]
#>     Rock fragment content (%): 21 Macroporosity (%): 6 
#>     Theta WP (%): 19 Theta FC (%): 33 Theta SAT (%): 44 Theta current (%) 33 
#>     Vol. WP (mm): 152 Vol. FC (mm): 263 Vol. SAT (mm): 347 Vol. current (mm): 263 
#>     Temperature (Celsius): NA 
#> 
#> Total soil saturated capacity (mm): 723 
#> Total soil water holding capacity (mm): 532 
#> Total soil extractable water (mm): 280 
#> Total soil current Volume (mm): 532 
#> Saturated water depth (mm): NA
```

It is important to remember that SoilGrids may underestimate the amount
of rocks in the soil. This is because soil samples (which were used to
generate the global database) do not normally contain large stones or
blocks. Hence, realistic simulations should reduce the soil water
holding capacity by increasing the column `rfc`. For example, here we
will assume that the third layer contains 80% of rocks:

``` r
pobl_soil_props$rfc[3] <- 80
```

If we rebuild the soil object and inspect its properties we will see the
effect on the soil water holding capacity and soil extractable water:

``` r
pobl_soil <- soil(pobl_soil_props)
summary(pobl_soil)
#> Soil depth (mm): 2000 
#> 
#> Layer  1  [ 0  to  300 mm ] 
#>     clay (%): 26 silt (%): 38 sand (%): 31 organic matter (%): 4 [ Loam ]
#>     Rock fragment content (%): 18 Macroporosity (%): 22 
#>     Theta WP (%): 18 Theta FC (%): 34 Theta SAT (%): 51 Theta current (%) 34 
#>     Vol. WP (mm): 44 Vol. FC (mm): 83 Vol. SAT (mm): 125 Vol. current (mm): 83 
#>     Temperature (Celsius): NA 
#> 
#> Layer  2  [ 300  to  1000 mm ] 
#>     clay (%): 30 silt (%): 39 sand (%): 30 organic matter (%): 1 [ Clay loam ]
#>     Rock fragment content (%): 19 Macroporosity (%): 9 
#>     Theta WP (%): 19 Theta FC (%): 33 Theta SAT (%): 44 Theta current (%) 33 
#>     Vol. WP (mm): 106 Vol. FC (mm): 186 Vol. SAT (mm): 250 Vol. current (mm): 186 
#>     Temperature (Celsius): NA 
#> 
#> Layer  3  [ 1000  to  2000 mm ] 
#>     clay (%): 32 silt (%): 38 sand (%): 30 organic matter (%): 1 [ Clay loam ]
#>     Rock fragment content (%): 80 Macroporosity (%): 6 
#>     Theta WP (%): 19 Theta FC (%): 33 Theta SAT (%): 44 Theta current (%) 33 
#>     Vol. WP (mm): 38 Vol. FC (mm): 67 Vol. SAT (mm): 88 Vol. current (mm): 67 
#>     Temperature (Celsius): NA 
#> 
#> Total soil saturated capacity (mm): 463 
#> Total soil water holding capacity (mm): 336 
#> Total soil extractable water (mm): 179 
#> Total soil current Volume (mm): 336 
#> Saturated water depth (mm): NA
```

## Interpolating weather

While soil information is often scarce and uncertain, obtaining daily
weather data suitable for simulations is not straightforward either.
Here we illustrate one way of obtaining such data using package
**meteoland**. We begin by adding topographic variables into a `sf`
object:

``` r
pobl_spt <- sf::st_sf(sf_pt) |>
            dplyr::mutate(elevation = 850,
                          slope = 15.1, 
                          aspect = 15)

pobl_spt
#> Simple feature collection with 1 feature and 3 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 1.0219 ymin: 41.3443 xmax: 1.0219 ymax: 41.3443
#> Geodetic CRS:  WGS 84
#>                    sf_pt elevation slope aspect
#> 1 POINT (1.0219 41.3443)       850  15.1     15
```

The more difficult part of using package **meteoland** is to assemble
reference weather data from surface weather stations into the socalled
**interpolator** object (of class `stars`). Please see the [meteoland
package documentation](https://emf-creaf.github.io/meteoland/index.html)
to learn how to create interpolator objects. Here we will assume that
such an object is already available, by using the example object
provided in the **meteoland** package.

``` r
data("meteoland_interpolator_example")
```

Once we have this interpolator, obtaining interpolated weather for a set
of target points is rather straightforward using function
[`interpolate_data()`](https://emf-creaf.github.io/meteoland/reference/interpolate_data.html)
from **meteoland**:

``` r
meteo <- interpolate_data(pobl_spt, meteoland_interpolator_example)
#> ℹ Starting interpolation...
#> ℹ Temperature interpolation is needed also...
#> • Interpolating temperature...
#> ℹ Precipitation interpolation is needed also...
#> • Interpolating precipitation...
#> ℹ Relative humidity interpolation is needed also...
#> • Interpolating relative humidity...
#> ℹ Radiation calculation is needed also...
#> • Calculating radiation...
#> ℹ Wind interpolation is needed also...
#> • Interpolating wind...
#> • Calculating PET...
#> ✔ Interpolation done...
```

The output of function
[`interpolate_data()`](https://emf-creaf.github.io/meteoland/reference/interpolate_data.html)
is an object of class `sf`:

``` r
meteo
#> Simple feature collection with 1 feature and 4 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 1.0219 ymin: 41.3443 xmax: 1.0219 ymax: 41.3443
#> Geodetic CRS:  WGS 84
#> # A tibble: 1 × 5
#>              sf_pt elevation slope aspect interpolated_data 
#>        <POINT [°]>     <dbl> <dbl>  <dbl> <list>            
#> 1 (1.0219 41.3443)       850  15.1     15 <tibble [30 × 13]>
```

We can access the weather data frame by subsetting the appropriate
element of `interpolated_data`:

``` r
pobl_weather <- meteo$interpolated_data[[1]]
head(pobl_weather)
#> # A tibble: 6 × 13
#>   dates                 DOY MeanTemperature MinTemperature MaxTemperature
#>   <dttm>              <dbl>           <dbl>          <dbl>          <dbl>
#> 1 2022-04-01 00:00:00    91            3.37         -2.21            6.99
#> 2 2022-04-02 00:00:00    92            3.60         -4.01            8.54
#> 3 2022-04-03 00:00:00    93            2.33         -7.67            8.83
#> 4 2022-04-04 00:00:00    94            4.16         -4.46            9.76
#> 5 2022-04-05 00:00:00    95            5.86         -5.21           13.1 
#> 6 2022-04-06 00:00:00    96            9.16          0.472          14.8 
#> # ℹ 8 more variables: Precipitation <dbl>, MeanRelativeHumidity <dbl>,
#> #   MinRelativeHumidity <dbl>, MaxRelativeHumidity <dbl>, Radiation <dbl>,
#> #   WindSpeed <dbl>, WindDirection <dbl>, PET <dbl>
```
