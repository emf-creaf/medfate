# Single-cohort forests

Creates a
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object with a single plant cohort

## Usage

``` r
tree2forest(
  Species,
  Height,
  LAI = NA,
  N = NA,
  DBH = NA,
  Z50 = NA,
  Z95 = NA,
  Z100 = NA,
  CrownRatio = NA,
  FoliarBiomass = NA,
  FuelLoading = NA
)

shrub2forest(
  Species,
  Height,
  LAI = NA,
  Cover = NA,
  Z50 = NA,
  Z95 = NA,
  Z100 = NA,
  CrownRatio = NA,
  FoliarBiomass = NA,
  FuelLoading = NA
)
```

## Arguments

- Species:

  String with species (taxon) name or a non-negative integer for species
  identity (i.e., 0,1,2,...) matching SpParams.

- Height:

  Plant height (cm).

- LAI:

  Leaf area index (m2/m2)

- N:

  Tree density (ind/ha)

- DBH:

  Tree DBH (cm).

- Z50:

  Depth (in mm) corresponding to 50% of fine roots.

- Z95:

  Depth (in mm) corresponding to 95% of fine roots.

- Z100:

  Depth (in mm) corresponding to 100% of fine roots.

- CrownRatio:

  Crown ratio (fraction of total height)

- FoliarBiomass:

  Standing dry biomass of leaves (kg/m2)

- FuelLoading:

  Fine fuel loading (kg/m2)

- Cover:

  Percent cover

## Value

An object of class
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)

## See also

[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`emptyforest`](https://emf-creaf.github.io/medfate/reference/emptyforest.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
oak_forest <-tree2forest("Quercus ilex", Height= 200, LAI = 2)
oak_forest
#> $treeData
#>        Species DBH Height  N Z50 Z95 LAI
#> 1 Quercus ilex  NA    200 NA  NA  NA   2
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
#> $seedlingBank
#> [1] Species Percent Age     Z50     Z95     Z100   
#> <0 rows> (or 0-length row.names)
#> 
#> $seedBank
#> [1] Species Percent
#> <0 rows> (or 0-length row.names)
#> 
#> attr(,"class")
#> [1] "forest" "list"  
```
