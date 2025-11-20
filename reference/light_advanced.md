# Advanced radiation transfer functions

Functions `light_layerDirectIrradianceFraction` and
`light_layerDiffuseIrradianceFraction` calculate the fraction of
above-canopy direct and diffuse radiation reaching each vegetation
layer. Function `light_layerSunlitFraction` calculates the proportion of
sunlit leaves in each vegetation layer. Function
`light_cohortSunlitShadeAbsorbedRadiation` calculates the amount of
radiation absorbed by cohort and vegetation layers, while
differentiating between sunlit and shade leaves.

## Usage

``` r
light_leafAngleCDF(leafAngle, p, q)

light_leafAngleBetaParameters(leafAngle, leafAngleSD)

light_directionalExtinctionCoefficient(p, q, solarElevation)

light_layerDirectIrradianceFraction(
  LAIme,
  LAImd,
  LAImx,
  kb,
  ClumpingIndex,
  alpha,
  gamma,
  trunkExtinctionFraction = 0.1
)

light_layerDiffuseIrradianceFraction(
  LAIme,
  LAImd,
  LAImx,
  K,
  ClumpingIndex,
  ZF,
  alpha,
  gamma,
  trunkExtinctionFraction = 0.1
)

light_cohortSunlitShadeAbsorbedRadiation(
  Ib0,
  Id0,
  LAIme,
  LAImd,
  LAImx,
  kb,
  K,
  ClumpingIndex,
  ZF,
  alpha,
  gamma,
  trunkExtinctionFraction = 0.1
)

light_layerSunlitFraction(LAIme, LAImd, kb, ClumpingIndex)

light_instantaneousLightExtinctionAbsortion(
  LAIme,
  LAImd,
  LAImx,
  p,
  q,
  ClumpingIndex,
  alphaSWR,
  gammaSWR,
  ddd,
  ntimesteps = 24L,
  trunkExtinctionFraction = 0.1
)

light_longwaveRadiationSHAW(
  LAIme,
  LAImd,
  LAImx,
  LWRatm,
  Tsoil,
  Tair,
  trunkExtinctionFraction = 0.1
)
```

## Arguments

- leafAngle:

  Average leaf inclination angle (in radians).

- p, q:

  Parameters of the beta distribution for leaf angles

- leafAngleSD:

  Standard deviation of leaf inclination angle (in radians).

- solarElevation:

  Solar elevation (in radians).

- LAIme:

  A numeric matrix of live expanded LAI values per vegetation layer
  (row) and cohort (column).

- LAImd:

  A numeric matrix of dead LAI values per vegetation layer (row) and
  cohort (column).

- LAImx:

  A numeric matrix of maximum LAI values per vegetation layer (row) and
  cohort (column).

- kb:

  A vector of direct light extinction coefficients.

- ClumpingIndex:

  The extent to which foliage has a nonrandom spatial distribution.

- alpha:

  A vector of leaf absorbance by species.

- gamma:

  A vector of leaf reflectance values.

- trunkExtinctionFraction:

  Fraction of extinction due to trunks (for winter deciduous forests).

- K:

  A vector of light extinction coefficients.

- ZF:

  Fraction of sky angles.

- Ib0:

  Above-canopy direct incident radiation.

- Id0:

  Above-canopy diffuse incident radiation.

- alphaSWR:

  A vecfor of hort-wave absorbance coefficients for each cohort.

- gammaSWR:

  A vector of short-wave reflectance coefficients (albedo) for each
  cohort.

- ddd:

  A dataframe with direct and diffuse radiation for different subdaily
  time steps (see function `radiation_directDiffuseDay` in package
  meteoland).

- ntimesteps:

  Number of subdaily time steps.

- LWRatm:

  Atmospheric downward long-wave radiation (W/m2).

- Tsoil:

  Soil temperature (Celsius).

- Tair:

  Canopy layer air temperature vector (Celsius).

## Value

Functions `light_layerDirectIrradianceFraction`,
`light_layerDiffuseIrradianceFraction` and `light_layerSunlitFraction`
return a numeric vector of length equal to the number of vegetation
layers.

Function `light_cohortSunlitShadeAbsorbedRadiation` returns a list with
two elements (matrices): `I_sunlit` and `I_shade`.

## Details

Functions for short-wave radiation are adapted from Anten & Bastiaans
(2016), whereas long-wave radiation balance follows Flerchinger et al.
(2009). Vegetation layers are assumed to be ordered from bottom to top.

## References

Anten, N.P.R., Bastiaans, L., 2016. The use of canopy models to analyze
light competition among plants, in: Hikosaka, K., Niinemets, U., Anten,
N.P.R. (Eds.), Canopy Photosynthesis: From Basics to Application.
Springer, pp. 379–398.

Flerchinger, G. N., Xiao, W., Sauer, T. J., Yu, Q. 2009. Simulation of
within-canopy radiation exchange. NJAS - Wageningen Journal of Life
Sciences 57 (1): 5–15. https://doi.org/10.1016/j.njas.2009.07.004.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`light_basic`](https://emf-creaf.github.io/medfate/reference/light_basic.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
solarElevation <- 0.67 # in radians
SWR_direct <- 1100
SWR_diffuse <- 300
PAR_direct <- 550
PAR_diffuse <- 150

LAI <- 2
nlayer <- 10
LAIlayerlive <- matrix(rep(LAI/nlayer,nlayer),nlayer,1)
LAIlayerdead <- matrix(0,nlayer,1)
meanLeafAngle <- 60 # in degrees
sdLeafAngle <- 20

beta <- light_leafAngleBetaParameters(meanLeafAngle*(pi/180), sdLeafAngle*(pi/180))

## Extinction coefficients
kb <- light_directionalExtinctionCoefficient(beta["p"], beta["q"], solarElevation)
kd_PAR <- 0.5
kd_SWR <- kd_PAR/1.35
```
