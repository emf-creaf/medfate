# Fire severity functions

Functions to estimate fire effects on foliage, buds and cambium, based
on the model by Michaletz & Johnson (2008)

## Usage

``` r
fire_plumeTemperature(Ib_surf, z, T_air = 25, rho_air = 1.169)

fire_barkThermalDiffusivity(fmc_bark, rho_bark = 500, T_air = 25)

fire_radialBoleNecrosis(
  Ib_surf,
  t_res,
  bark_diffusivity,
  T_air = 25,
  rho_air = 1.169,
  T_necrosis = 60
)

fire_leafThermalFactor(SLA, h = 130, c = 2500)

fire_necrosisCriticalTemperature(
  t_res,
  thermal_factor,
  T_air = 25,
  T_necrosis = 60
)

fire_necrosisHeight(
  Ib_surf,
  t_res,
  thermal_factor,
  T_air = 25,
  rho_air = 1.169,
  T_necrosis = 60
)
```

## Arguments

- Ib_surf:

  Surface fireline intensity (kW/m).

- z:

  height (m).

- T_air:

  Air temperature (degrees Celsius).

- rho_air:

  Air density (kg/m3).

- fmc_bark:

  Bark moisture content (% dry weight).

- rho_bark:

  Bark density (kg/m3).

- t_res:

  fire residence time (seconds).

- bark_diffusivity:

  Bark thermal diffusivity (m2/s).

- T_necrosis:

  Temperature of tissue necrosis (degrees Celsius).

- SLA:

  Specific leaf area (m2/kg).

- h:

  Heat transfer coefficient

- c:

  Specific heat capacity

- thermal_factor:

  Tissue thermal factor.

## Value

- Function `fire_plumeTemperature` returns the plume temperature at a
  given height.

- Function `fire_barkThermalDiffusivity` returns the bark thermal
  diffusivity given a bark moisture value.

- Function `fire_radialBoleNecrosis` returns the depth of radial bole
  necrosis in cm.

- Function `fire_leafThermalFactor` returns the thermal factor of leaves
  as a function of specific leaf area.

- Function `fire_necrosisCriticalTemperature` returns the (plume)
  temperature yielding necrosis for a given residence time and tissue
  thermal factor.

- Function `fire_necrosisHeight` returns the height (in m) of necrosis
  for tissues with given thermal factor.

## References

Michaletz, S.T., and Johnson, E.A. 2006. A heat transfer model of crown
scorch in forest fires. Can. J. For. Res. 36: 2839–2851.
doi:10.1139/X06-158.

Michaletz ST, Johnson EA. 2008. A biophysical process model of tree
mortality in surface fires. Canadian Journal of Forest Research 38:
2013–2029.
