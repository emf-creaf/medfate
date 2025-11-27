# Drought stress indicators

Calculates plant drought stress indices, at different temporal scales,
from simulation results.

## Usage

``` r
droughtStress(x, index = "NDD", freq = "years", bySpecies = FALSE, draw = TRUE)
```

## Arguments

- x:

  An object of class
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/pwb.md),
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md) or
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

- index:

  A string with the index to be calculated, either `"DI"`, `"NDD"`,
  `"ADS"`, `"MDS"` or `"WSI"` (see details).

- freq:

  Frequency of stress statistics (see
  [`cut.Date`](https://rdrr.io/r/base/cut.POSIXt.html)). Normally,
  either `"years"` or `"months"` for yearly-based or monthly-based
  indices.

- bySpecies:

  Allows aggregating output by species.

- draw:

  A boolean flag to indicate that a plot should be returned.

## Value

A data frame with periods (e.g., years or months) in rows and plant
cohorts (or species) in columns. Values are the calculated stress index.
If `draw=TRUE` a ggplot is returned instead.

## Details

The currently available drought stress indices are:

- `"ADS"`: Average of daily drought stress values for the period
  considered.

- `"MDS"`: Maximum daily drought stress during the period considered.

- `"DI"`: Drought intensity, as defined in De Cáceres et al. (2015).

- `"NDD"`: Number of drought days, as defined in De Cáceres et al.
  (2015).

- `"WSI"`: Water stress integral, as defined in Myers (1988).

## References

De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
forest inventory data to predict drought stress: the role of forest
structural changes vs. climate changes. Agricultural and Forest
Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).

Myers BJ (1988) Water stress integral - a link between short-term stress
and long-term growth. Tree Physiology 4: 315–323 (doi:
10.1093/treephys/4.4.315)

## See also

[`summary.spwb`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md),
[`waterUseEfficiency`](https://emf-creaf.github.io/medfate/reference/waterUseEfficiency.md)

## Author

Miquel De Cáceres Ainsa, CREAF
