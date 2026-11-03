# Soil-plant resistances

Calculates and draws rhizosphere, root, stem and leaf resistances for
simulation time steps

## Usage

``` r
resistances(
  x,
  cohort,
  relative = FALSE,
  draw = FALSE,
  cumulative = FALSE,
  xlab = NULL,
  ylab = NULL
)
```

## Arguments

- x:

  An object of class
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/pwb.md),
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md) or
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).
  The function only works with the result of simulations with
  `transpirationMode = "Sperry"` or `transpirationMode = "Sureau"`.

- cohort:

  An string indicating the cohort for which resistances are desired.

- relative:

  A boolean flag to indicate that relative percentages are desired as
  output.

- draw:

  A boolean flag to indicate that a plot should be drawn (only pathway
  resistances, without discriminating between soil layers).

- cumulative:

  A flag to indicate that drawn series should be cumulative.

- xlab:

  x-axis label.

- ylab:

  y-axis label.

## Value

If `draw = FALSE`, the function returns list with three items:

- `pathway`: A data frame with dates in rows and resistance segments in
  columns (Rhizosphere, Root, Stem and Leaf).

- `root`: A data frame with dates in rows and root resistances for soil
  layers in columns.

- `rhizosphere`: A data frame with dates in rows and rhizosphere
  resistances for soil layers in columns.

Values depend on whether `relative = TRUE` (percentages) or
`relative = FALSE` (absolute resistance values).

If `draw = TRUE`, a plot object is returned showing the time series of
pathway resistances.

## Details

The function makes internal calls to
[`hydraulics_soilPlantResistancesWeibull`](https://emf-creaf.github.io/medfate/reference/hydraulics_scalingconductance.md)
or
[`hydraulics_soilPlantResistancesSigmoid`](https://emf-creaf.github.io/medfate/reference/hydraulics_scalingconductance.md)
depending on the value of `transpirationMode` in `x`.

## See also

[`waterUseEfficiency`](https://emf-creaf.github.io/medfate/reference/waterUseEfficiency.md),
[`droughtStress`](https://emf-creaf.github.io/medfate/reference/droughtStress.md)

## Author

Miquel De Cáceres Ainsa, CREAF

Léa Veuillen, INRAE-URFM
