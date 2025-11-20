# Water use efficiency

Calculates plant water use efficiency (WUE), at different temporal
scales, from simulation results.

## Usage

``` r
waterUseEfficiency(
  x,
  type = "Plant Ag/E",
  leaves = "average",
  freq = "days",
  draw = TRUE,
  ylim = NULL
)
```

## Arguments

- x:

  An object of class
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md) or
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

- type:

  A string to indicate the scale of WUE calculation. Either:

  - `"Leaf iWUE"`: Leaf intrinsic WUE, i.e. instantaneous ratio between
    photosynthesis and stomatal conductance (only for simulations with
    `transpirationMode = "Sperry"` or `transpirationMode = "Sureau"` and
    `subdailyResults = TRUE`).

  - `"Leaf Ci"`: Leaf intercellular CO2 concentration (only for
    simulations with `transpirationMode = "Sperry"` or
    `transpirationMode = "Sureau"` and `subdailyResults = TRUE`).

  - `"Plant An/E"`: Plant (cohort) net photosynthesis over plant
    transpiration (only for simulations with
    `transpirationMode = "Sperry"` or `transpirationMode = "Sureau"`)

  - `"Stand An/E"`: Stand net photosynthesis over stand transpiration
    (only for simulations with `transpirationMode = "Sperry"` or
    `transpirationMode = "Sureau"`)

  - `"Plant Ag/E"`: Plant (cohort) gross photosynthesis over plant
    transpiration

  - `"Stand Ag/E"`: Stand gross photosynthesis over stand transpiration

- leaves:

  Either `"sunlit"`, `"shade"` or `"average"`. Refers to the WUE of
  different leaf types or the average (with weights according to the LAI
  of sunlit and shade leaves). Only relevant for `type = "iWUE"`.

- freq:

  Frequency of summary statistics (see
  [`cut.Date`](https://rdrr.io/r/base/cut.POSIXt.html)).

- draw:

  A boolean flag to indicate that a plot should be returned.

- ylim:

  Range of values for y.

## Value

If `draw=TRUE` a plot is returned. Otherwise, the function returns a
matrix with WUE values, where rows are dates (at the desired temporal
scale), and columns are plant cohorts. In the case of
`type = "Plant Ag/E"`, `type = "Stand Ag/E"`, `type = "Plant An/E"` and
`type = "Stand An/E"` values are in gC/L. In the case of
`type = "Leaf iWUE"` values are in micromol of carbon per mmol of water.

## Details

Temporal aggregation of WUE values is done differently depending on the
value of `type`. For `type = "Plant Ag/E"`, `type = "Stand Ag/E"`,
`type = "Plant An/E"` and `type = "Stand An/E"` sums or daily
photosynthesis and transpiration are first calculated at the desired
temporal scale and the ratio is calculated afterwards. For
`type = "Leaf iWUE"` intrinsic WUE values are first calculated at the
daily scale (as averages of instantaneous An/gs ratios weighted by An)
and then they are aggregated to the desired scale by calculating
weighted averages, where weights are given by daily photosynthesis.

## See also

[`droughtStress`](https://emf-creaf.github.io/medfate/reference/droughtStress.md)

## Author

Miquel De Cáceres Ainsa, CREAF
