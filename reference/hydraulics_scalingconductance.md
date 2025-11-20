# Scaling from conductivity to conductance

Functions used to scale from tissue conductivity to conductance of
different elements of the continuum.

## Usage

``` r
hydraulics_maximumSoilPlantConductance(krhizomax, krootmax, kstemmax, kleafmax)

hydraulics_soilPlantResistancesSigmoid(
  psiSoil,
  psiRhizo,
  psiStem,
  PLCstem,
  psiLeaf,
  PLCleaf,
  krhizomax,
  n,
  alpha,
  krootmax,
  root_P50,
  root_slope,
  kstemmax,
  stem_P50,
  stem_slope,
  kleafmax,
  leaf_P50,
  leaf_slope
)

hydraulics_soilPlantResistancesWeibull(
  psiSoil,
  psiRhizo,
  psiStem,
  PLCstem,
  psiLeaf,
  PLCleaf,
  krhizomax,
  n,
  alpha,
  krootmax,
  rootc,
  rootd,
  kstemmax,
  stemc,
  stemd,
  kleafmax,
  leafc,
  leafd
)

hydraulics_averageRhizosphereResistancePercent(
  krhizomax,
  n,
  alpha,
  krootmax,
  rootc,
  rootd,
  kstemmax,
  stemc,
  stemd,
  kleafmax,
  leafc,
  leafd,
  psiStep = -0.01
)

hydraulics_findRhizosphereMaximumConductance(
  averageResistancePercent,
  n,
  alpha,
  krootmax,
  rootc,
  rootd,
  kstemmax,
  stemc,
  stemd,
  kleafmax,
  leafc,
  leafd,
  initialValue = 0
)

hydraulics_taperFactorSavage(height)

hydraulics_terminalConduitRadius(height)

hydraulics_referenceConductivityHeightFactor(refheight, height)

hydraulics_maximumStemHydraulicConductance(
  xylemConductivity,
  refheight,
  Al2As,
  height,
  taper = FALSE
)

hydraulics_rootxylemConductanceProportions(L, V)
```

## Arguments

- krhizomax:

  Maximum rhizosphere hydraulic conductance (defined as flow per leaf
  surface unit and per pressure drop).

- krootmax:

  Maximum root xylem hydraulic conductance (defined as flow per leaf
  surface unit and per pressure drop).

- kstemmax:

  Maximum stem xylem hydraulic conductance (defined as flow per leaf
  surface unit and per pressure drop).

- kleafmax:

  Maximum leaf hydraulic conductance (defined as flow per leaf surface
  unit and per pressure drop).

- psiSoil:

  Soil water potential (in MPa). A scalar or a vector depending on the
  function.

- psiRhizo:

  Water potential (in MPa) in the rhizosphere (root surface).

- psiStem:

  Water potential (in MPa) in the stem.

- PLCstem:

  Percent loss of conductance (in %) in the stem.

- psiLeaf:

  Water potential (in MPa) in the leaf.

- n, alpha:

  Parameters of the Van Genuchten function (rhizosphere vulnerability
  curve).

- root_P50, root_slope:

  Parameters of the Sigmoid function for roots (root xylem vulnerability
  curve).

- stem_P50, stem_slope:

  Parameters of the Sigmoid function for stems (stem xylem vulnerability
  curve).

- leaf_P50, leaf_slope:

  Parameters of the Sigmoid function for leaves (leaf vulnerability
  curve).

- rootc, rootd:

  Parameters of the Weibull function for roots (root xylem vulnerability
  curve).

- stemc, stemd:

  Parameters of the Weibull function for stems (stem xylem vulnerability
  curve).

- leafc, leafd:

  Parameters of the Weibull function for leaves (leaf vulnerability
  curve).

- psiStep:

  Water potential precision (in MPa).

- averageResistancePercent:

  Average (across water potential values) resistance percent of the
  rhizosphere, with respect to total resistance (rhizosphere + root
  xylem + stem xylem).

- initialValue:

  Initial value of rhizosphere conductance.

- height:

  Plant height (in cm).

- refheight:

  Reference plant height of measurement of xylem conductivity (in cm).

- xylemConductivity:

  Xylem conductivity as flow per length of conduit and pressure drop (in
  kg·m-1·s-1·MPa-1).

- Al2As:

  Leaf area to sapwood area (in m2·m-2).

- taper:

  A boolean flag to indicate correction by taper of xylem conduits
  (Christoffersen et al. 2017).

- L:

  Vector with the length of coarse roots (mm) for each soil layer.

- V:

  Vector with the proportion \[0-1\] of fine roots within each soil
  layer.

## Value

Values returned for each function are:

- `hydraulics_maximumSoilPlantConductance`: The maximum soil-plant
  conductance, in the same units as the input segment conductances.

- `hydraulics_averageRhizosphereResistancePercent`: The average
  percentage of resistance due to the rhizosphere, calculated across
  water potential values.

- `hydraulics_findRhizosphereMaximumConductance`: The maximum
  rhizosphere conductance value given an average rhizosphere resistance
  and the vulnerability curves of rhizosphere, root and stem elements.

- `hydraulics_taperFactorSavage`: Taper factor according to Savage et
  al. (2010).

## Details

Details of the hydraulic model are given in the medfate book

## References

Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R.
Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A.
Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and
P. Meir. 2016. Linking hydraulic traits to tropical forest function in a
size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific
Model Development Discussions 9: 4227–4255.

Savage, V. M., L. P. Bentley, B. J. Enquist, J. S. Sperry, D. D. Smith,
P. B. Reich, and E. I. von Allmen. 2010. Hydraulic trade-offs and space
filling enable better predictions of vascular structure and function in
plants. Proceedings of the National Academy of Sciences of the United
States of America 107:22722–7.

Olson, M.E., Anfodillo, T., Rosell, J.A., Petit, G., Crivellaro, A.,
Isnard, S., León-Gómez, C., Alvarado-Cárdenas, L.O., and Castorena, M.
2014. Universal hydraulics of the flowering plants: Vessel diameter
scales with stem length across angiosperm lineages, habits and climates.
Ecology Letters 17: 988–997.

## See also

[`hydraulics_psi2K`](https://emf-creaf.github.io/medfate/reference/hydraulics_conductancefunctions.md),
[`hydraulics_supplyFunctionPlot`](https://emf-creaf.github.io/medfate/reference/hydraulics_supplyfunctions.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)

## Author

Miquel De Cáceres Ainsa, CREAF
