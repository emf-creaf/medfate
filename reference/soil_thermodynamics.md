# Soil thermodynamic functions

Functions `soil_thermalConductivity` and `soil_thermalCapacity`
calculate thermal conductivity and thermal capacity for each soil layer,
given its texture and water content. Functions
`soil_temperatureGradient` and `soil_temperatureChange` are used to
calculate soil temperature gradients (in ºC/m) and temporal temperature
change (in ºC/s) given soil layer texture and water content (and
possibly including heat flux from above).

## Usage

``` r
soil_thermalCapacity(soil, model = "SX")

soil_thermalConductivity(soil, model = "SX")

soil_temperatureGradient(widths, Temp)

soil_temperatureChange(
  widths,
  Temp,
  sand,
  clay,
  W,
  Theta_SAT,
  Theta_FC,
  Gdown,
  tstep
)
```

## Arguments

- soil:

  Soil object (returned by function
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)).

- model:

  Either 'SX' or 'VG' for Saxton's or Van Genuchten's pedotransfer
  models.

- widths:

  Width of soil layers (in mm).

- Temp:

  Temperature (in ºC) for each soil layer.

- sand:

  Percentage of sand (in percent weight) for each layer.

- clay:

  Percentage of clay (in percent weight) for each layer.

- W:

  Soil moisture (in percent of field capacity) for each layer.

- Theta_SAT:

  Relative water content (in percent volume) at saturation for each
  layer.

- Theta_FC:

  Relative water content (in percent volume) at field capacity for each
  layer.

- Gdown:

  Downward heat flux from canopy to soil (in W·m-2).

- tstep:

  Time step (interval) in seconds.

## Value

Function `soil_thermalConductivity` returns a vector with values of
thermal conductivity (W/m/ºK) for each soil layer.

Function `soil_thermalCapacity` returns a vector with values of heat
storage capacity (J/m3/ºK) for each soil layer.

Function `soil_temperatureGradient` returns a vector with values of
temperature gradient between consecutive soil layers.

Function `soil_temperatureChange` returns a vector with values of
instantaneous temperature change (ºC/s) for each soil layer.

## References

Cox, P.M., Betts, R.A., Bunton, C.B., Essery, R.L.H., Rowntree, P.R.,
and Smith, J. 1999. The impact of new land surface physics on the GCM
simulation of climate and climate sensitivity. Climate Dynamics 15:
183–203.

Dharssi, I., Vidale, P.L., Verhoef, A., MacPherson, B., Jones, C., and
Best, M. 2009. New soil physical properties implemented in the Unified
Model at PS18. 9–12.

## See also

[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Define soil and complete parameters
examplesoil = soil(defaultSoilParams(4))

soil_thermalConductivity(examplesoil)
#> [1] 2.2 2.2 2.2 2.2
soil_thermalCapacity(examplesoil)
#> [1] 2486214 2486214 2486214 2486214

#Values change when altering water content (drier layers have lower conductivity and capacity)
examplesoil$W = c(0.1, 0.4, 0.7, 1.0)
soil_thermalConductivity(examplesoil)
#> [1] 0.4266171 1.4943000 1.9252995 2.2000000
soil_thermalCapacity(examplesoil)
#> [1] 1342121 1723486 2104850 2486214
```
