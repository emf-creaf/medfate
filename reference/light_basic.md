# Radiation extinction functions used in basic transpiration sub-model

Radiation extinction functions used in basic transpiration sub-model

## Usage

``` r
light_PARcohort(x, SpParams, gdd = NA_real_)

light_PARground(x, SpParams, gdd = NA_real_)

light_SWRground(x, SpParams, gdd = NA_real_)

light_cohortAbsorbedSWRFraction(z, x, SpParams, gdd = NA_real_)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- gdd:

  Growth degree days.

- z:

  A numeric vector with height values.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`light_advanced`](https://emf-creaf.github.io/medfate/reference/light_advanced.md)

## Author

Miquel De Cáceres Ainsa, CREAF
