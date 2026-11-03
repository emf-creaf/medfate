# Mortality

A simple sigmoid function to determine a daily mortality likelihood
according to the value of a stress variable.

## Usage

``` r
mortality_dailyProbability(stressValue, stressThreshold)
```

## Arguments

- stressValue:

  Current value of the stress variable (0 to 1, with higher values
  indicate stronger stress).

- stressThreshold:

  Threshold to indicate 50% annual mortality probability.

## Value

Returns a probability (between 0 and 1)

## See also

[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
[`regeneration`](https://emf-creaf.github.io/medfate/reference/regeneration.md)

## Author

Miquel De Cáceres Ainsa, CREAF
