# Physical and biophysical utility functions

Internal utility functions for the calculation of biophysical variables.

## Usage

``` r
biophysics_radiationDiurnalPattern(t, daylength)

biophysics_temperatureDiurnalPattern(
  t,
  tmin,
  tmax,
  tminPrev,
  tmaxPrev,
  tminNext,
  daylength
)

biophysics_leafTemperature(absRad, airTemperature, u, E, leafWidth = 1)

biophysics_leafTemperature2(
  SWRabs,
  LWRnet,
  airTemperature,
  u,
  E,
  leafWidth = 1
)

biophysics_leafVapourPressure(leafTemp, leafPsi)

biophysics_irradianceToPhotonFlux(I, lambda = 546.6507)

biophysics_waterDynamicViscosity(temp)
```

## Arguments

- t:

  Time of the day (in seconds).

- daylength:

  Day length (in seconds).

- tmin, tmax:

  Minimum and maximum daily temperature (ºC).

- tminPrev, tmaxPrev, tminNext:

  Maximum and minimum daily temperatures of the previous and following
  day (ºC).

- absRad:

  Absorbed long- and short-wave radiation (in W·m-2).

- airTemperature:

  Air temperature (in ºC).

- u:

  Wind speed above the leaf boundary layer (in m/s).

- E:

  Transpiration flow (in mmol H20·m-2·s-1) per one sided leaf area
  basis.

- leafWidth:

  Leaf width (in cm).

- SWRabs:

  Absorbed short-wave radiation (in W·m-2).

- LWRnet:

  Net long-wave radiation balance (in W·m-2).

- leafTemp:

  Leaf temperature (ºC).

- leafPsi:

  Leaf water potential (MPa).

- I:

  Irradiance (in W\*m-2).

- lambda:

  Wavelength (in nm).

- temp:

  Temperature (ºC).

## Value

Values returned for each function are:

- `biophysics_leafTemperature` and `biophysics_leafTemperature2`: leaf
  temperature (in ºC)

- `biophysics_leafVapourPressure`: leaf vapour pressure (in kPa)

- `biophysics_radiationDiurnalPattern`: the proportion of daily
  radiation corresponding to the input time in seconds after sunrise.

- `biophysics_temperatureDiurnalPattern`: diurnal pattern of
  temperature.

- `biophysics_waterDynamicViscosity`: Water dynamic viscosity relative
  to 20ºC.

## Details

Functions `biophysics_leafTemperature` and `biophysics_leafTemperature2`
calculate leaf temperature according to energy balance equation given in
Campbell and Norman (1988).

Function `biophysics_radiationDiurnalPattern` follows the equations
given in Liu and Jordan (1960).

Function `biophysics_temperatureDiurnalPattern` determines diurnal
temperature pattern assuming a sinusoidal pattern with T = Tmin at
sunrise and T = (Tmin+Tmax)/2 at sunset and a linear change in
temperature between sunset and Tmin of the day after (McMurtrie et al.
1990).

Function `biophysics_waterDynamicViscosity` calculates water dynamic
viscosity following the Vogel (1921) equation.

## References

Campbell, G. S., and J. M. Norman. 1998. An introduction to
environmental biophysics: 2nd edition. (eqns. 14.1 & 14.3)

B. Y. H. Liu and R. C. Jordan, “The interrelationship and characteristic
distribution of direct, diffuse and total solar radiation,” Solar
Energy, vol. 4, no. 3, pp. 1–19, 1960.

McMurtrie, R. E., D. A. Rook, and F. M. Kelliher. 1990. Modelling the
yield of Pinus radiata on a site limited by water and nitrogen. Forest
Ecology and Management 30:381–413.

H. Vogel, "Das Temperaturabhangigkeitsgesetz der Viskositat von
Flussigkeiten", Physikalische Zeitschrift, vol. 22, pp. 645–646, 1921.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)

## Author

Miquel De Cáceres Ainsa, CREAF
