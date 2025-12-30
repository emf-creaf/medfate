# DAYCENT decomposition

Functions implementing a modification of the DAYCENT carbon
decomposition model (Parton et al. 1988, 1993, 1998), inspired by the
description given in Chapter 18 of Bonan (2019). Functions
`decompositionDAYCENTsnags` and `decompositionDAYCENTlitter` conduct
snag and litter decomposition, respectively, whereas function
`decomposition_DAYCENT` performs the whole model for carbon
decomposition.

## Usage

``` r
decomposition_DAYCENTsnags(
  snags,
  baseAnnualRates,
  airTemperature,
  airRelativeHumidity,
  tstep = 1
)

decomposition_DAYCENTlitter(
  litter,
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
  snags,
  litter,
  SOC,
  paramsLitterDecomposition,
  baseAnnualRates,
  annualTurnoverRate,
  airTemperature,
  airRelativeHumidity,
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

- snags:

  A data frame with dead standing (snag) cohort information (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- baseAnnualRates:

  A named vector of annual decomposition rates, in yr-1 (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- airTemperature:

  Mean daily air temperature (in Celsius).

- airRelativeHumidity:

  Mean daily relative humidity (%).

- tstep:

  Time step in days. By default, one day. For annual time steps, use
  `tstep = 365.25`.

- litter:

  A data frame with aboveground and belowground structural carbon pools
  corresponding to plant cohorts, in g C/m2 (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- paramsLitterDecomposition:

  A data frame of species-specific litter decomposition parameters (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

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

- SOC:

  A named numeric vector with metabolic, active, slow and passive carbon
  pools for surface and soil, in g C/m2 (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- annualTurnoverRate:

  Annual turnover rate, in yr-1 (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

## Value

Function `decomposition_DAYCENTsnags` returns a vector containing
transfer carbon flows to SOC pools and heterotrophic respiration from
snag decomposition. Function `decomposition_DAYCENTlitter` returns a
vector containing transfer carbon flows to SOC pools and heterotrophic
respiration from litter decomposition. Function `decomposition_DAYCENT`
returns scalar value with heterotrophic respiration (snags + litter +
soil), in g C/m2.

## Details

Each call to functions `decomposition_DAYCENTlitter` or
`decomposition_DAYCENTsnags` conducts one time step of the snag or
litter dynamics, respectively. Each call to function
`decomposition_DAYCENT` conducts one time step of the whole DAYCENT
model and returns the heterotrophic respiration for that day.

*IMPORTANT NOTE*: Decomposition functions modify the input data (i.e.
`snags`, `litter` and/or `SOC`) according to decomposition rates and
carbon transfer rates. When used as part of
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)
simulations, soil physical and chemical parameters correspond to the
uppermost soil layer.

## References

Bonan, G. (2019). Climate change and terrestrial ecosystem modeling.
Cambridge University Press, Cambridge, UK.

Parton WJ, Steward JWB, Cole CV (1988). Dynamics of C, N, P and S in
grassland soils: a model. Biogeochemistry 5: 109-131.

Parton WJ, Scurlock JMO, Ojima DS, Gilmanov TG, Scoles RJ et al. (1993).
Observations and modeling of biomass and soil organic matter dynamics
for the grassland biome worldwide. Global Biogeochemical Cycles 7:
785-809.

Parton WJ, Hartman M, Ojima DS, Schimel D (1998). DAYCENT and its land
surface submodel: Description and testing. Global and Planetary Change,
19: 35-48.

## See also

[`decomposition_temperatureEffect`](https://emf-creaf.github.io/medfate/reference/decomposition_annualLitterDecompositionRate.md),
[`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)

## Author

Miquel De Cáceres Ainsa, CREAF
