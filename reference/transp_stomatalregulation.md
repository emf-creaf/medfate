# Stomatal regulation

Set of high-level functions used in the calculation of stomatal
conductance and transpiration. Function `transp_profitMaximization`
calculates gain and cost functions, as well as profit maximization from
supply and photosynthesis input functions. Function
`transp_stomatalRegulationPlot` produces a plot with the cohort supply
functions against water potential and a plot with the cohort
photosynthesis functions against water potential, both with the maximum
profit values indicated.

## Usage

``` r
transp_profitMaximization(
  supplyFunction,
  photosynthesisFunction,
  Gswmin,
  Gswmax
)
```

## Arguments

- supplyFunction:

  Water supply function (see
  [`hydraulics_supplyFunctionNetwork`](https://emf-creaf.github.io/medfate/reference/hydraulics_supplyfunctions.md)).

- photosynthesisFunction:

  Function returned by `photo_photosynthesisFunction()`.

- Gswmin, Gswmax:

  Minimum and maximum stomatal conductance to water vapour
  (mol·m-2·s-1).

## Value

Function `transp_profitMaximization` returns a list with the following
elements:

- `Cost`: Cost function \[0-1\].

- `Gain`: Gain function \[0-1\].

- `Profit`: Profit function \[0-1\].

- `iMaxProfit`: Index corresponding to maximum profit (starting from 0).

## References

Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S.
Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to
the environment from the optimization of photosynthetic gain and
hydraulic cost. Plant Cell and Environment 40, 816-830 (doi:
10.1111/pce.12852).

## See also

[`transp_transpirationSperry`](https://emf-creaf.github.io/medfate/reference/transp_modes.md),
[`hydraulics_supplyFunctionNetwork`](https://emf-creaf.github.io/medfate/reference/hydraulics_supplyfunctions.md),
[`biophysics_leafTemperature`](https://emf-creaf.github.io/medfate/reference/biophysics.md),
[`photo_photosynthesis`](https://emf-creaf.github.io/medfate/reference/photo.md),
[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
[`plot.spwb_day`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md)

## Author

Miquel De Cáceres Ainsa, CREAF
