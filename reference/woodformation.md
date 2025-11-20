# Wood formation

Functions to initialize and expand a ring of tracheids to simulate
secondary growth.

## Usage

``` r
woodformation_initRing()

woodformation_temperatureEffect(
  Tc,
  Y_T = 5,
  DHa = 87500,
  DSd = 1090,
  DHd = 333000
)

woodformation_relativeExpansionRate(psi, Tc, pi, phi, Y_P, Y_T)

woodformation_growRing(
  ring,
  psi,
  Tc,
  Nc = 8.85,
  phi0 = 0.13,
  pi0 = -0.8,
  CRD0 = 8.3,
  Y_P = 0.05,
  Y_T = 5,
  h = 0.043 * 1.8,
  s = 1.8
)

woodformation_relativeGrowthRate(dbh1, dbh2, yeardiff, lower = -2, upper = 8)
```

## Arguments

- Tc:

  Temperature in Celsius.

- Y_T:

  Temperature yield threshold (in Celsius)

- DHa, DSd, DHd:

  Enthalpy of activation, enthalpy difference and entropy difference
  between the catalytically active and inactive states of the enzymatic
  system (Parent et al. 2010).

- psi:

  Water potential (in MPa).

- pi:

  Osmotic potential (in MPa)

- phi:

  Cell extensibility (in MPa-1 day-1)

- Y_P:

  Turgor pressure yield threshold (in MPa)

- ring:

  An object of class [`ring`](https://rdrr.io/r/grDevices/plotmath.html)
  returned by function `woodformation_initRing`.

- Nc:

  Number of active cells in the cambium.

- phi0:

  Initial value of cell extensibility (in MPa-1 day-1)

- pi0:

  Initial value of cell osmotic potential (in MPa)

- CRD0:

  Initial value of cell radial diameter

- h:

  Cell wall hardening coefficient (in day-1)

- s:

  Cell wall softening coefficient (unitless)

- dbh1, dbh2:

  Initial and final diameter at breast height.

- yeardiff:

  Interval between dbh measurements, in years.

- lower, upper:

  Lower and upper bounds for root finding.

## Value

Function `woodformation_initRing()` returns a list of class 'ring', that
is a list containing a data frame `cells` and two vectors: `P` and `SA`.
Dataframe `cells` contains the columns "formation_date", "phi", "pi" and
"CRD" and as many rows as dates processed. Vectors `P` and `SA` contain,
respectively, the number of cells produced and the sapwood area
corresponding to the ring of cells (assuming a tangencial radius of 20
micrometers).

Function `woodformation_growRing()` modifies the input 'ring' object
according to the environmental conditions given as input.

Function `woodformation_relativeExpansionRate()` returns a numeric
scalar with the relative expansion rate.

Function `woodformation_temperatureEffect()` returns a scalar between 0
and 1 reflecting the temperature effect on tissue formation rate.

Function `woodformation_relativeGrowthRate` returns the annual growth
rate, relative to cambium perimeter, estimated from initial and final
diameter values.

## Note

Code modified from package xylomod by Antoine Cabon, available at GitHub

## References

Cabon A, Fernández-de-Uña L, Gea-Izquierdo G, Meinzer FC, Woodruff DR,
Martínez-Vilalta J, De Cáceres M. 2020a. Water potential control of
turgor-driven tracheid enlargement in Scots pine at its xeric
distribution edge. New Phytologist 225: 209–221.

Cabon A, Peters RL, Fonti P, Martínez-Vilalta J, De Cáceres M. 2020b.
Temperature and water potential co-limit stem cambial activity along a
steep elevational gradient. New Phytologist: nph.16456.

Parent, B., O. Turc, Y. Gibon, M. Stitt, and F. Tardieu. 2010. Modelling
temperature-compensated physiological rates, based on the co-ordination
of responses to temperature of developmental processes. Journal of
Experimental Botany 61:2057–2069.

## See also

[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)

## Author

Antoine Cabon, CTFC

Miquel De Cáceres Ainsa, CREAF
