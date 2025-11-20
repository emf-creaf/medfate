# Photosynthesis submodel functions

Set of functions used in the calculation of photosynthesis

## Usage

``` r
photo_GammaTemp(Tleaf)

photo_KmTemp(Tleaf, Oi = 209)

photo_VmaxTemp(Vmax298, Tleaf)

photo_JmaxTemp(Jmax298, Tleaf)

photo_electronLimitedPhotosynthesis(Q, Ci, GT, Jmax)

photo_rubiscoLimitedPhotosynthesis(Ci, GT, Km, Vmax)

photo_photosynthesis(Q, Catm, Gc, Tleaf, Vmax298, Jmax298, verbose = FALSE)

photo_photosynthesisBaldocchi(
  Q,
  Catm,
  Tleaf,
  u,
  Vmax298,
  Jmax298,
  leafWidth,
  Gsw_AC_slope,
  Gsw_AC_intercept
)

photo_leafPhotosynthesisFunction(
  E,
  psiLeaf,
  Catm,
  Patm,
  Tair,
  vpa,
  u,
  absRad,
  Q,
  Vmax298,
  Jmax298,
  leafWidth = 1,
  refLeafArea = 1,
  verbose = FALSE
)

photo_leafPhotosynthesisFunction2(
  E,
  psiLeaf,
  Catm,
  Patm,
  Tair,
  vpa,
  u,
  SWRabs,
  LWRnet,
  Q,
  Vmax298,
  Jmax298,
  leafWidth = 1,
  refLeafArea = 1,
  verbose = FALSE
)

photo_sunshadePhotosynthesisFunction(
  E,
  psiLeaf,
  Catm,
  Patm,
  Tair,
  vpa,
  SLarea,
  SHarea,
  u,
  absRadSL,
  absRadSH,
  QSL,
  QSH,
  Vmax298SL,
  Vmax298SH,
  Jmax298SL,
  Jmax298SH,
  leafWidth = 1,
  verbose = FALSE
)

photo_multilayerPhotosynthesisFunction(
  E,
  psiLeaf,
  Catm,
  Patm,
  Tair,
  vpa,
  SLarea,
  SHarea,
  u,
  absRadSL,
  absRadSH,
  QSL,
  QSH,
  Vmax298,
  Jmax298,
  leafWidth = 1,
  verbose = FALSE
)
```

## Arguments

- Tleaf:

  Leaf temperature (in ºC).

- Oi:

  Oxigen concentration (mmol\*mol-1).

- Vmax298, Vmax298SL, Vmax298SH:

  Maximum Rubisco carboxylation rate per leaf area at 298ºK (i.e. 25 ºC)
  (micromol*s-1*m-2) (for each canopy layer in the case of
  `photo_multilayerPhotosynthesisFunction`). 'SH' stands for shade
  leaves, whereas 'SL' stands for sunlit leaves.

- Jmax298, Jmax298SL, Jmax298SH:

  Maximum electron transport rate per leaf area at 298ºK (i.e. 25 ºC)
  (micromol*s-1*m-2) (for each canopy layer in the case of
  `photo_multilayerPhotosynthesisFunction`). 'SH' stands for shade
  leaves, whereas 'SL' stands for sunlit leaves.

- Q:

  Active photon flux density (micromol \* s-1 \* m-2).

- Ci:

  CO2 internal concentration (micromol \* mol-1).

- GT:

  CO2 saturation point corrected by temperature (micromol \* mol-1).

- Jmax:

  Maximum electron transport rate per leaf area (micromol*s-1*m-2).

- Km:

  Km = Kc\*(1.0+(Oi/Ko)) - Michaelis-Menten term corrected by
  temperature (in micromol \* mol-1).

- Vmax:

  Maximum Rubisco carboxylation rate per leaf area (micromol*s-1*m-2).

- Catm:

  CO2 air concentration (micromol \* mol-1).

- Gc:

  CO2 leaf (stomatal) conductance (mol \* s-1 \* m-2).

- verbose:

  Boolean flag to indicate console output.

- u:

  Wind speed above the leaf boundary (in m/s) (for each canopy layer in
  the case of `photo_multilayerPhotosynthesisFunction`).

- leafWidth:

  Leaf width (in cm).

- Gsw_AC_slope:

  Slope of the An/C vs Gsw relationship

- Gsw_AC_intercept:

  Intercept of the An/C vs Gsw relationship

- E:

  Transpiration flow rate per leaf area (mmol*s-1*m-2).

- psiLeaf:

  Leaf water potential (MPa).

- Patm:

  Atmospheric air pressure (in kPa).

- Tair:

  Air temperature (in ºC).

- vpa:

  Vapour pressure deficit (in kPa).

- absRad:

  Absorbed long- and short-wave radiation (in W\*m^-2).

- refLeafArea:

  Leaf reference area.

- SWRabs:

  Absorbed short-wave radiation (in W·m-2).

- LWRnet:

  Net long-wave radiation balance (in W·m-2).

- SLarea, SHarea:

  Leaf area index of sunlit/shade leaves (for each canopy layer in the
  case of `photo_multilayerPhotosynthesisFunction`).

- absRadSL, absRadSH:

  Instantaneous absorbed radiation (W·m-2) per unit of sunlit/shade leaf
  area (for each canopy layer in the case of
  `photo_multilayerPhotosynthesisFunction`).

- QSL, QSH:

  Active photon flux density (micromol \* s-1 \* m-2) per unit of
  sunlit/shade leaf area (for each canopy layer in the case of
  `photo_multilayerPhotosynthesisFunction`).

## Value

Values returned for each function are:

- `photo_GammaTemp`: CO2 compensation concentration (micromol \* mol-1).

- `photo_KmTemp`: Michaelis-Menten coefficients of Rubisco for Carbon
  (micromol \* mol-1) and Oxigen (mmol \* mol-1).

- `photo_VmaxTemp`: Temperature correction of Vmax298.

- `photo_JmaxTemp`: Temperature correction of Jmax298.

- `photo_electronLimitedPhotosynthesis`: Electron-limited photosynthesis
  (micromol*s-1*m-2) following Farquhar et al. (1980).

- `photo_rubiscoLimitedPhotosynthesis`: Rubisco-limited photosynthesis
  (micromol*s-1*m-2) following Farquhar et al. (1980).

- `photo_photosynthesis`: Calculates gross photosynthesis
  (micromol*s-1*m-2) following (Farquhar et al. (1980) and Collatz et al
  (1991).

- `photo_leafPhotosynthesisFunction`: Returns a data frame with the
  following columns:

  - `LeafTemperature`: Leaf temperature (ºC).

  - `LeafVPD`: Leaf vapor pressure deficit (kPa).

  - `LeafCi`: Internal CO2 concentration (micromol \* mol-1).

  - `Gsw`: Leaf stomatal conductance to water vapor (mol \* s-1 \* m-2).

  - `GrossPhotosynthesis`: Gross photosynthesis (micromol*s-1*m-2).

  - `NetPhotosynthesis`: Net photosynthesis, after discounting
    autotrophic respiration (micromol*s-1*m-2).

- `photo_sunshadePhotosynthesisFunction`: Returns a data frame with the
  following columns:

  - `GrossPhotosynthesis`: Gross photosynthesis (micromol*s-1*m-2).

  - `NetPhotosynthesis`: Net photosynthesis, after discounting
    autotrophic respiration (micromol*s-1*m-2).

  - `LeafCiSL`: Sunlit leaf internal CO2 concentration (micromol \*
    mol-1).

  - `LeafCiSH`: Shade leaf internal CO2 concentration (micromol \*
    mol-1).

  - `LeafTempSL`: Sunlit leaf temperature (ºC).

  - `LeafTempSH`: Shade leaf temperature (ºC).

  - `LeafVPDSL`: Sunlit leaf vapor pressure deficit (kPa).

  - `LeafVPDSH`: Shade leaf vapor pressure deficit (kPa).

- `photo_multilayerPhotosynthesisFunction`: Return a data frame with the
  following columns:

  - `GrossPhotosynthesis`: Gross photosynthesis (micromol*s-1*m-2).

  - `NetPhotosynthesis`: Net photosynthesis, after discounting
    autotrophic respiration (micromol*s-1*m-2).

## Details

Details of the photosynthesis submodel are given in the medfate book

## References

Bernacchi, C. J., E. L. Singsaas, C. Pimentel, A. R. Portis, and S. P.
Long. 2001. Improved temperature response functions for models of
Rubisco-limited photosynthesis. Plant, Cell and Environment 24:253–259.

Collatz, G. J., J. T. Ball, C. Grivet, and J. A. Berry. 1991.
Physiological and environmental regulation of stomatal conductance,
photosynthesis and transpiration: a model that includes a laminar
boundary layer. Agricultural and Forest Meteorology 54:107–136.

Farquhar, G. D., S. von Caemmerer, and J. A. Berry. 1980. A biochemical
model of photosynthetic CO2 assimilation in leaves of C3 species. Planta
149:78–90.

Leuning, R. 2002. Temperature dependence of two parameters in a
photosynthesis model. Plant, Cell and Environment 25:1205–1210.

Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S.
Mackay, Y. Wang, and D. M. Love. 2016. Predicting stomatal responses to
the environment from the optimization of photosynthetic gain and
hydraulic cost. Plant Cell and Environment.

## See also

[`hydraulics_supplyFunctionNetwork`](https://emf-creaf.github.io/medfate/reference/hydraulics_supplyfunctions.md),
[`biophysics_leafTemperature`](https://emf-creaf.github.io/medfate/reference/biophysics.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)

## Author

Miquel De Cáceres Ainsa, CREAF
