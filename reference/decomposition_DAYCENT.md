# DAYCENT decomposition

Functions implementing the DAYCENT carbon decomposition model, following
the description in Bonan (2019). Function `decompositionDAYCENTlitter`
conducts litter decomposition only, whereas function
`decomposition_DAYCENT` performs the whole model for carbon
decomposition.

## Usage

``` r
decomposition_DAYCENTlitter(
  structuralLitter,
  paramsLitterDecomposition,
  baseAnnualRates,
  sand,
  clay,
  soilTemperature,
  soilMoisture,
  soilPH,
  soilO2 = 1,
  cultfac = 1,
  tstep = 1
)

decomposition_DAYCENT(
  structuralLitter,
  SOC,
  paramsLitterDecomposition,
  baseAnnualRates,
  annualTurnoverRate,
  sand,
  clay,
  soilTemperature,
  soilMoisture,
  soilPH,
  soilO2 = 1,
  cultfac = 1,
  tstep = 1
)
```

## Arguments

- structuralLitter:

  A data frame with structural carbon pools corresponding to plant
  cohorts, in g C/m2 (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- paramsLitterDecomposition:

  A data frame of species-specific decomposition parameters (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- baseAnnualRates:

  A named vector of annual decomposition rates, in yr-1 (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- sand, clay:

  Soil texture (sand and sand) in percent volume (%).

- soilTemperature:

  Soil temperature (in Celsius).

- soilMoisture:

  Soil moisture content, relative to saturation (0-1).

- soilPH:

  Soil pH (0-14).

- soilO2:

  Soil oxygen factor (0-1).

- cultfac:

  Cultivation factor (0-1).

- tstep:

  Time step in days. By default, one day. For annual time steps, use
  `tstep = 365.25`.

- SOC:

  A named numeric vector with metabolic, active, slow and passive carbon
  pools for surface and soil, in g C/m2 (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- annualTurnoverRate:

  Annual turnover rate, in yr-1 (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

## Value

Function `decomposition_DAYCENTlitter` returns a vector containing
transfer carbon flows to SOC pools and heterotrophic respiration from
litter decomposition. Function `decomposition_DAYCENT` returns scalar
value with heterotrophic respiration (litter + soil), in g C/m2.

## Details

Each call to function `decomposition_DAYCENTlitter` conducts one time
step of the litter dynamics in DAYCENT. Each call to function
`decomposition_DAYCENT` conducts one time step of the whole DAYCENT
model and returns the heterotrophic respiration for that day. Both
functions modify input data `structuralLitter` (and in case case of
`decomposition_DAYCENT` also `SOC`) according to decomposition rates and
carbon transfer rates. When used as part of
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)
simulations, soil physical and chemical characteristics correspond to
the uppermost soil layer.

## References

Bonan, G. (2019). Climate change and terrestrial ecosystem modeling.
Cambridge University Press, Cambridge, UK.

## See also

[`decomposition_temperatureEffect`](https://emf-creaf.github.io/medfate/reference/decomposition_annualLitterDecompositionRate.md),
[`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)

## Author

Miquel De Cáceres Ainsa, CREAF
