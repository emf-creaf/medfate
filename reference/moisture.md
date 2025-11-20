# Tissue moisture functions

Set of functions used to calculate tissue moisture from water potential
and viceversa.

## Usage

``` r
moisture_sapwoodWaterCapacity(Al2As, height, V, L, wd)

moisture_leafWaterCapacity(SLA, ld)

moisture_turgorLossPoint(pi0, epsilon)

moisture_symplasticRWC(psiSym, pi0, epsilon)

moisture_symplasticPsi(RWC, pi0, epsilon)

moisture_apoplasticRWC(psiApo, c, d)

moisture_apoplasticPsi(RWC, c, d)

moisture_tissueRWC(psiSym, pi0, epsilon, psiApo, c, d, af)

plant_water(x)

moisture_pressureVolumeCurvePlot(
  x,
  segment = "stem",
  fraction = "all",
  psiVec = seq(-0.1, -8, by = -0.01),
  speciesNames = FALSE
)
```

## Arguments

- Al2As:

  Leaf area to sapwood area (in m2·m-2).

- height:

  Plant height (in cm).

- V:

  Vector with the proportion \[0-1\] of fine roots within each soil
  layer.

- L:

  Vector with the length of coarse roots (mm) for each soil layer.

- wd:

  Wood density (g·cm-3).

- SLA:

  Specific leaf area (mm2·mg-1).

- ld:

  Leaf tissue density (g·cm-3).

- pi0:

  Full turgor osmotic potential (MPa).

- epsilon:

  Bulk modulus of elasticity (MPa).

- psiSym, psiApo:

  Symplastic or apoplastic water potential (MPa).

- RWC:

  Relative water content \[0-1\].

- c, d:

  Parameters of the xylem vulnerability curve.

- af:

  Apoplastic fraction (proportion) in the segment (e.g. leaf or stem).

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  or
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- segment:

  Segment whose relative water content curve to plot, either `"stem"` or
  `"leaf"` (the latter only available if `transpirationMode = "Sperry"`
  or `transpirationMode = "Sureau"`).

- fraction:

  Tissue fraction, either `"symplastic"`, `"apoplastic"` or `"all"`.

- psiVec:

  Vector of water potential values to evaluate for the pressure-volume
  curve.

- speciesNames:

  A flag to indicate the use of species names instead of cohort names in
  plots.

## Value

Values returned for each function are:

- `moisture_symplasticRWC`: Relative water content \[0-1\] of the
  symplastic fraction.

- `moisture_apoplasticRWC`: Relative water content \[0-1\] of the
  apoplastic fraction.

- `moisture_symplasticWaterPotential`: Water potential (in MPa) of the
  symplastic fraction.

- `moisture_apoplasticWaterPotential`: Water potential (in MPa) of the
  apoplastic fraction.

- `moisture_turgorLossPoint`: Water potential (in MPa) corresponding to
  turgor loss point.

- `moisture_segmentRWC`: Segment relative water content \[0-1\].

- `water_plant`: A vector of water content (mm) per plant cohort.

## References

Bartlett, M.K., Scoffoni, C., Sack, L. 2012. The determinants of leaf
turgor loss point and prediction of drought tolerance of species and
biomes: a global meta-analysis. Ecology Letters 15: 393–405.

Hölttä, T., Cochard, H., Nikinmaa, E., Mencuccini, M. 2009. Capacitive
effect of cavitation in xylem conduits: Results from a dynamic model.
Plant, Cell and Environment 32: 10–21.

Martin-StPaul, N., Delzon, S., Cochard, H. 2017. Plant resistance to
drought depends on timely stomatal closure. Ecology Letters 20:
1437–1447.

## See also

[`hydraulics_psi2K`](https://emf-creaf.github.io/medfate/reference/hydraulics_conductancefunctions.md),
[`hydraulics_supplyFunctionPlot`](https://emf-creaf.github.io/medfate/reference/hydraulics_supplyfunctions.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
psi = seq(-10,0, by=0.1)
rwc_s = rep(NA, length(psi))
for(i in 1:length(psi)) rwc_s[i] = moisture_symplasticRWC(psi[i],-3,12)
plot(psi, rwc_s, type="l", xlab="Water potential (MPa)", ylab = "Symplasmic RWC")

```
