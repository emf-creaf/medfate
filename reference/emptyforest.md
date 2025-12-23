# Creation of an empty forest

Creates an empty
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object.

## Usage

``` r
emptyforest(
  ntree = 0,
  nshrub = 0,
  nherb = 0,
  nseedling = 0,
  nseed = 0,
  nsnag = 0,
  nlitter = 0,
  addcolumns = NULL,
  SOC = FALSE
)
```

## Arguments

- ntree, nshrub, nherb:

  Number of tree, shrub and herb cohorts, respectively.

- nseedling:

  Number of species in the seedling bank.

- nseed:

  Number of species in the seed bank.

- nsnag:

  Number of snag (dead standing) cohorts.

- nlitter:

  Number of items in the litter compartment.

- addcolumns:

  A character vector with additional columns. Currently allowed are (see
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)):

  - `Z100`: A numeric vector with maximum rooting depth in mm.

  - `LAI`: Leaf area index (m2/m2).

  - `FoliarBiomass`: Standing dry biomass of leaves (kg/m2).

  - `FuelLoading`: Fine fuel loading (kg/m2).

  - `CrownRatio`: The ratio between crown length and total height
    (between 0 and 1)

  - `Age`: A numeric vector indicating age of cohorts in years.

  - `ObsID` A character vector to label specific cohorts in simulations
    of forest dynamics.

- SOC:

  A boolean flag to initialize SOC data (see
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)).

## Value

An empty
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object.

## Details

List elements `treeData` and `shrubData` are always created, regardless
of the number of cohorts. In contrast, list elements `herbData`,
`seedBank` and `seedlingBank` are only created if `nherb`, `nseed` and
`nseedling` are non-zero, respectively.

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
# Initializes forest with 2 tree cohorts
emptyforest(ntree = 2)
#> $treeData
#>   Species DBH Height  N Z50 Z95
#> 1    <NA>  NA     NA NA  NA  NA
#> 2    <NA>  NA     NA NA  NA  NA
#> 
#> $shrubData
#> [1] Species Height  Cover   Z50     Z95    
#> <0 rows> (or 0-length row.names)
#> 
#> attr(,"class")
#> [1] "forest" "list"  

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
#> attr(,"class")
#> [1] "forest" "list"  

# Initializes forest with 2 tree cohorts, 1 shrub cohort and 1 herbaceous cohort
emptyforest(ntree = 2, nshrub = 1, nherb = 1)
#> $treeData
#>   Species DBH Height  N Z50 Z95
#> 1    <NA>  NA     NA NA  NA  NA
#> 2    <NA>  NA     NA NA  NA  NA
#> 
#> $shrubData
#>   Species Height Cover Z50 Z95
#> 1    <NA>     NA    NA  NA  NA
#> 
#> $herbData
#>   Species Height Cover Z50 Z95
#> 1    <NA>     NA    NA  NA  NA
#> 
#> attr(,"class")
#> [1] "forest" "list"  
```
