# Hydraulic-related defoliation

Functions to calculate the proportion of crown defoliation due to
hydraulic disconnection.

## Usage

``` r
hydraulics_proportionDefoliationSigmoid(
  psiLeaf,
  P50,
  slope,
  PLC_crit = 0.88,
  P50_cv = 10
)

hydraulics_proportionDefoliationWeibull(
  psiLeaf,
  c,
  d,
  PLC_crit = 0.88,
  P50_cv = 10
)
```

## Arguments

- psiLeaf:

  Leaf water potential (in MPa).

- P50, slope:

  Parameters of the Sigmoid function.

- PLC_crit:

  Critical leaf PLC corresponding to defoliation

- P50_cv:

  Coefficient of variation (in percent) of leaf P50, to describe the
  variability in hydraulic vulnerability across crown leaves.

- c, d:

  Parameters of the Weibull function.

## Value

The proportion of crown defoliation.

## Details

The functions assume that crowns are made of a population of leaves
whose hydraulic vulnerability (i.e. the water potential corresponding to
50% loss of conductance) follows a Gaussian distribution centered on the
input P50 and with a known coefficient of variation (`P50_cv`). The
slope parameter (or the c exponent in the case of a Weibull function) is
considered constant. Leaves are hydraulically disconnected, and shedded,
when their embolism rate exceeds a critical value (`PLC_crit`).

## See also

[`hydraulics_conductancefunctions`](https://emf-creaf.github.io/medfate/reference/hydraulics_conductancefunctions.md)

## Author

Hervé Cochard, INRAE

Miquel De Cáceres Ainsa, CREAF
