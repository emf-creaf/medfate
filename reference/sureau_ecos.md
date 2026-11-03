# Sureau-ECOS inner functions for testing only

Function `initSureauNetworks` initializes hydraulic networks for all
plant cohorts in x Function `semi_implicit_integration` updates water
potentials and cavitation across the hydraulic network

## Usage

``` r
initSureauNetworks(x)

semi_implicit_integration(
  network,
  dt,
  opt,
  stemCavitationRecovery = "annual",
  leafCavitationRecovery = "total"
)
```

## Arguments

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  or
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  created using `transpirationMode = "Sureau"`.

- network:

  A hydraulic network element of the list returned by
  `initSureauNetworks`

- dt:

  Smallest time step (seconds)

- opt:

  Option flag vector

- stemCavitationRecovery, leafCavitationRecovery:

  A string indicating how refilling of embolized conduits is done:

  - "none" - no refilling.

  - "annual" - every first day of the year.

  - "rate" - following a rate of new sapwood formation.

  - "total" - instantaneous complete refilling.

## Value

Function `initSureauNetworks` returns a vector of length equal to the
number of cohorts. Each element is a list with Sureau-ECOS parameters.
Function `semi_implicit_integration` does not return anything, but
modifies input parameter `network`.

## References

Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022)
SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations
of plant water status and drought-induced mortality at the ecosystem
level. Geoscientific Model Development 15, 5593-5626
(doi:10.5194/gmd-15-5593-2022).

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)

## Author

- Miquel De Cáceres Ainsa, CREAF

- Nicolas Martin-StPaul, URFM-INRAE
