# Summary of forest structure

Displays a summary of forest structure

## Usage

``` r
# S3 method for class 'forest'
summary(object, SpParams, ...)

# S3 method for class 'summary.forest'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- object:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- ...:

  Additional parameters for functions
  [`summary`](https://rdrr.io/r/base/summary.html) and
  [`print`](https://rdrr.io/r/base/print.html).

- x:

  The object returned by `summary.forest`.

- digits:

  Minimal number of significant digits.

## Value

Function `summary.forest` returns a list with several structural
attributes, such as the basal area and LAI of the forest.

## Details

Function `summary.forest` can be used to summarize a `forest` object in
the console.

## See also

[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`forest_mapWoodyTables`](https://emf-creaf.github.io/medfate/reference/forest_mapWoodyTables.md),
[`forest_mergeTrees`](https://emf-creaf.github.io/medfate/reference/forest_simplification.md),
[`plot.forest`](https://emf-creaf.github.io/medfate/reference/plot.forest.md),
[`tree2forest`](https://emf-creaf.github.io/medfate/reference/tree2forest.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
# Summary of example forests
summary(exampleforest, SpParamsMED)
#> Tree BA (m2/ha): 25.0333016  adult trees: 25.0333016  saplings: 0 
#> Density (ind/ha) adult trees: 552  saplings: 0  shrubs (estimated): 749.4923076 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 3.75  herbs: 0 
#> LAI (m2/m2) total: 1.4516629  adult trees: 1.3983416  saplings: 0  shrubs: 0.0533213  herbs: 0  mistletoe: 0 
#> Fuel loading (kg/m2) total: 0.4889222  adult trees: 0.471855  saplings: 0  shrubs: 0.0170672  herbs: 0 
#> PAR ground (%): 46.7337671  SWR ground (%): 56.9222511 
summary(exampleforest2, SpParamsMED)
#> Tree BA (m2/ha): 0  adult trees: 0  saplings: 0 
#> Density (ind/ha) adult trees: 0  saplings: 0  shrubs (estimated): 0 
#> Cover (%) adult trees: 0  saplings: 0  shrubs: 0  herbs: 0 
#> LAI (m2/m2) total: 1.58  adult trees: 0  saplings: 0  shrubs: 0.03  herbs: 0.25  mistletoe: 0 
#> Fuel loading (kg/m2) total: 0.4856741  adult trees: 0  saplings: 0  shrubs: 0.0096025  herbs: 0.0277778 
#> PAR ground (%): 50.0824269  SWR ground (%): 59.9162626 
```
