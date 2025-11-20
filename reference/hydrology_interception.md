# Rainfall interception

Function `hydrology_rainInterception` calculates the amount of rainfall
intercepted daily by the canopy, given a rainfall and canopy
characteristics. Two canopy interception models are currently available:
the sparse Gash (1995) model and the Liu (2001) model. In both cases the
current implementation assumes no trunk interception.

## Usage

``` r
hydrology_rainfallIntensity(month, prec, rainfallIntensityPerMonth)

hydrology_rainInterception(Rainfall, Cm, p, ER = 0.05, model = "Gash1995")

hydrology_interceptionPlot(
  x,
  SpParams,
  ER = 0.05,
  gdd = NA,
  throughfall = FALSE,
  model = "Gash1995"
)
```

## Arguments

- month:

  Month of the year (from 1 to 12).

- prec:

  Precipitation for a given day (mm).

- rainfallIntensityPerMonth:

  A vector with twelve positions with average intensity of rainfall (in
  mm/h) for each month.

- Rainfall:

  A numeric vector of (daily) rainfall.

- Cm:

  Canopy water storage capacity.

- p:

  Proportion of throughfall (normally 1 - c, where c is the canopy
  cover).

- ER:

  The ratio of evaporation rate to rainfall rate.

- model:

  Rainfall interception model (either `"Gash1995"` or `"Liu2001"`).

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  and
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- gdd:

  Growth degree days (in Celsius).

- throughfall:

  Boolean flag to plot relative throughfall instead of percentage of
  intercepted rainfall.

## Value

Function `hydrology_rainInterception` returns a vector of the same
length as `Rainfall` containing intercepted rain values.

Function `hydrology_rainfallIntensity` returns a scalar with the
rainfall intensity.

## Details

Function `hydrology_rainInterception` can accept either vectors or
scalars as parameters `Cm`, `p` and `ER`. If they are supplied as
vectors they should be of the same length as `Rainfall`.

Function `hydrology_rainfallIntensity` estimates the rainfall intensity
(mm/h) for input values of rainfall and seasonal variation in rainfall
intensity (mm/h).

## References

Liu (2001). Evaluation of the Liu model for predicting rainfall
interception in forests world-wide. - Hydrol. Process. 15: 2341-2360.

Gash (1979). An analytical model of rainfall interception by forests. -
Quarterly Journal of the Royal Meteorological Society.

Gash et al. (1995). Estimating sparse forest rainfall interception with
an analytical model. - Journal of Hydrology.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Draw rainfall interception for two values of the E/R ratio
hydrology_interceptionPlot(exampleforest, SpParamsMED, ER = c(0.05, 0.2))

```
