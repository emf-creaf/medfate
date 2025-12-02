# Creation of an empty forest

Creates an empty
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object.

## Usage

``` r
emptyforest(ntree = 0, nshrub = 0, nseedling = 0, nseed = 0, addcolumns = NULL)
```

## Arguments

- ntree, nshrub:

  Number of tree and shrub cohorts, respectively.

- nseedling:

  Number of species in the seedling bank.

- nseed:

  Number of species in the seed bank.

- addcolumns:

  A character vector with additional columns. Currently allowed are:

  - `Z100`: A numeric vector with maximum rooting depth in mm.

  - `LAI`: Leaf area index (m2/m2).

  - `FoliarBiomass`: Standing dry biomass of leaves (kg/m2).

  - `FuelLoading`: Fine fuel loading (kg/m2).

  - `CrownRatio`: The ratio between crown length and total height
    (between 0 and 1)

  - `Age`: A numeric vector indicating age of cohorts in years.

  - `ObsID` A character vector to label specific cohorts in simulations
    of forest dynamics.

## Value

An empty
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object.

## See also

[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`tree2forest`](https://emf-creaf.github.io/medfate/reference/tree2forest.md),
[`summary.forest`](https://emf-creaf.github.io/medfate/reference/summary.forest.md),
[`forest_mapWoodyTables`](https://emf-creaf.github.io/medfate/reference/forest_mapWoodyTables.md),
[`forest_mergeTrees`](https://emf-creaf.github.io/medfate/reference/forest_simplification.md),
[`plot.forest`](https://emf-creaf.github.io/medfate/reference/plot.forest.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
# Initializes forest with 2 tree cohorts and 1 shrub cohort
emptyforest(ntree = 2, nshrub = 1)
#> $treeData
#>   Species DBH Height  N Z50 Z95
#> 1    <NA>  NA     NA NA  NA  NA
#> 2    <NA>  NA     NA NA  NA  NA
#> 
#> $shrubData
#>   Species Height Cover Z50 Z95
#> 1    <NA>     NA    NA  NA  NA
#> 
#> $herbCover
#> [1] NA
#> 
#> $herbHeight
#> [1] NA
#> 
#> $seedlingBank
#> [1] Species Percent Age     Z50     Z95    
#> <0 rows> (or 0-length row.names)
#> 
#> $seedBank
#> [1] Species Percent
#> <0 rows> (or 0-length row.names)
#> 
#> attr(,"class")
#> [1] "forest" "list"  

# Initializes with extra column for leaf area index (LAI)
emptyforest(ntree = 2, nshrub = 1, addcolumns = "LAI")
#> $treeData
#>   Species DBH Height  N Z50 Z95 LAI
#> 1    <NA>  NA     NA NA  NA  NA  NA
#> 2    <NA>  NA     NA NA  NA  NA  NA
#> 
#> $shrubData
#>   Species Height Cover Z50 Z95 LAI
#> 1    <NA>     NA    NA  NA  NA  NA
#> 
#> $herbCover
#> [1] NA
#> 
#> $herbHeight
#> [1] NA
#> 
#> $seedlingBank
#> [1] Species Percent Age     Z50     Z95    
#> <0 rows> (or 0-length row.names)
#> 
#> $seedBank
#> [1] Species Percent
#> <0 rows> (or 0-length row.names)
#> 
#> attr(,"class")
#> [1] "forest" "list"  
```
