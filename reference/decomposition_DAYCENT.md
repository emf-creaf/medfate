# DAYCENT decomposition

Function implementing a modification of the DAYCENT carbon decomposition
model (Parton et al. 1988, 1993, 1998), inspired by the description
given in Chapter 18 of Bonan (2019).

## Usage

``` r
decomposition_DAYCENT(
  snags,
  litter,
  SOC,
  paramsLitterDecomposition,
  baseAnnualRates,
  annualTurnoverRate,
  environmentalConditions,
  litterProduction,
  sand,
  clay,
  soilPH,
  soilO2 = 1,
  cultfac = 1,
  tstep = 1
)
```

## Arguments

- snags:

  A data frame with initial dead standing (snag) cohort information (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- litter:

  A data frame with initial aboveground and belowground structural
  carbon pools corresponding to plant cohorts, in g C/m2 (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- SOC:

  A named numeric vector with initial metabolic, active, slow and
  passive carbon pools for surface and soil, in g C/m2 (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- paramsLitterDecomposition:

  A data frame of species-specific litter decomposition parameters (see
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)).

- baseAnnualRates:

  A named vector of annual decomposition rates, in yr-1 (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- annualTurnoverRate:

  Annual turnover rate, in yr-1 (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- environmentalConditions:

  A data frame containing environmental conditions for each time step to
  simulate:

  - `AirTempperature`: Mean air temperature (in Celsius).

  - `AirRelativeHumidity`: Mean relative humidity (percent).

  - `SoilTemperature`: Mean soil temperature (in Celsius).

  - `SoilMoisture`: Mean soil moisture, relative to saturation (0-1).

- litterProduction:

  A data frame containing litter inputs corresponding to time steps:

  - `Step`: Integer indicating the time step where litter is generated,
    corresponding to a row in `environmentalConditions`.

  - `Species`: Mean relative humidity (percent).

  - `Leaves`: Leaf litter production (g C·m-2).

  - `Twigs`: Twig litter production (g C·m-2).

  - `SmallBranches`: Small branch litter production (g C·m-2).

  - `LargeWood`: Large wood litter production (g C·m-2).

  - `CoarseRoots`: Coarse root litter production (g C·m-2).

  - `FineRoots`: Fine root litter production (g C·m-2).

- sand, clay:

  Soil texture (sand and sand) in percent volume (percent).

- soilPH:

  Soil pH (0-14).

- soilO2:

  Soil oxygen factor (0-1).

- cultfac:

  Cultivation factor (0-1).

- tstep:

  Time step in days. By default, one day. For annual time steps, use
  `tstep = 365.25`.

## Value

A list with two elements:

- `"HeterotrophicRespiration"`: A numeric vector with the heterotrophic
  respiration (g C·m-2) corresponding to the decomposition in each time
  step.

- `"DecompositionPools"`: A numeric matrix with the mass of different
  decomposition carbon pools, all in g C · m-2.

## Details

Species names must match between inputs `paramsLitterDecomposition`,
`litter` and `litterProduction`.

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

Roberto Molowny-Horas, CREAF

Inés Delsman Valderrama, CREAF
