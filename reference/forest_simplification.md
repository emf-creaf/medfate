# Forest complexity reduction

Functions `forest_mergeTrees` and `forest_mergeShrubs` merge cohorts of
a [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object. Function `forest_reduceToDominant` performs a strongest
simplification of plant cohorts (see details).

## Usage

``` r
forest_mergeTrees(x, byDBHclass = TRUE, keepCohortsWithObsID = FALSE)

forest_mergeShrubs(x, byHeightclass = TRUE, keepCohortsWithObsID = FALSE)

forest_reduceToDominant(x, SpParams)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- byDBHclass:

  Logical flag to indicate that 5-cm tree DBH classes should be kept
  separated.

- keepCohortsWithObsID:

  Logical flag to indicate that cohorts with non-missin ObsID should be
  spared from merging.

- byHeightclass:

  Boolean flag to indicate that 10-cm shrub height classes should be
  kept separated.

- SpParams:

  A data frame with species parameters (see
  [`SpParamsDefinition`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  and
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

## Value

Another
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object with simplified structure/composition, depending on the function.

## Details

Tree DBH classes are defined in 5-cm intervals, whereas shrub height
classes are defined in 10-cm intervals. Tree DBH and shrub height
classes are defined up to a specific size (i.e. larger plants are not
merged) corresponding to 52.5 cm and 90 cm, respectively.

Function `forest_reduceToDominant` simplifies the input forest to the
tree cohort of highest LAI, among those of the tree species with highest
LAI. The leaf area index of the whole tree layer will be attributed to
the chosen cohort. The same is performed for the shrub layer.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`forest_mapWoodyTables`](https://emf-creaf.github.io/medfate/reference/forest_mapWoodyTables.md),
[`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md),
[`summary.forest`](https://emf-creaf.github.io/medfate/reference/summary.forest.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
# Example forest data
data("exampleforest")

# Reduce to dominant tree and dominant shrub
reduced <- forest_reduceToDominant(exampleforest, SpParamsMED)

# Check that overall LAI does not change
stand_LAI(exampleforest, SpParamsMED)
#> [1] 1.584948
stand_LAI(reduced, SpParamsMED)
#> [1] 1.584948
```
