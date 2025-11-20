# Plot forest attributes

Convenient wrappers for vertical forest profiles (see
[`vprofile_leafAreaDensity`](https://emf-creaf.github.io/medfate/reference/vprofile_leafAreaDensity.md)).

## Usage

``` r
# S3 method for class 'forest'
plot(
  x,
  SpParams,
  type = "LeafAreaDensity",
  byCohorts = FALSE,
  bySpecies = FALSE,
  includeHerbs = FALSE,
  ...
)

# S3 method for class 'forest'
shinyplot(x, SpParams, ...)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- type:

  A string of the plot type: "LeafAreaDensity", "RootDistribution",
  "FuelBulkDensity", "PARExtinction", "SWRExtinction" or
  "WindExtinction".

- byCohorts:

  A logical flag to separate profiles for each cohort.

- bySpecies:

  A logical flag to aggregate results by species.

- includeHerbs:

  A logical flag to include herbaceous layer in the profile.

- ...:

  Additional parameters to vertical profiles

## Value

A ggplot or a shiny application, depending on the function.

## See also

[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`summary.forest`](https://emf-creaf.github.io/medfate/reference/summary.forest.md),
[`vprofile_leafAreaDensity`](https://emf-creaf.github.io/medfate/reference/vprofile_leafAreaDensity.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
data(exampleforest)
data(SpParamsMED)
plot(exampleforest, SpParamsMED)

```
