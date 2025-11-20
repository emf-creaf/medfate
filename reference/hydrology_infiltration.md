# Soil infiltration

Soil infiltration functions:

- Function `hydrology_infiltrationBoughton` calculates the amount of
  water that infiltrates into the topsoil, according to the USDA SCS
  curve number method (Boughton 1989).

- Function `hydrology_infiltrationGreenAmpt` calculates the amount of
  water that infiltrates into the topsoil, according to the model by
  Green & Ampt (1911).

- Function `hydrology_infiltrationAmount` uses either Green &
  Ampt (1911) or Boughton (1989) to estimate infiltration.

- Function `hydrology_infiltrationRepartition` distributes infiltration
  among soil layers depending on macroporosity.

## Usage

``` r
hydrology_infiltrationBoughton(input, Ssoil)

hydrology_infiltrationGreenAmpt(t, psi_w, Ksat, theta_sat, theta_dry)

hydrology_infiltrationRepartition(I, widths, macro, a = -0.005, b = 3)

hydrology_infiltrationAmount(
  rainfallInput,
  rainfallIntensity,
  soil,
  soilFunctions,
  model = "GreenAmpt1911",
  K_correction = 1
)
```

## Arguments

- input:

  A numeric vector of (daily) water input (in mm of water).

- Ssoil:

  Soil water storage capacity (can be referred to topsoil) (in mm of
  water).

- t:

  Time of the infiltration event

- psi_w:

  Matric potential at the wetting front

- Ksat:

  hydraulic conductivity at saturation

- theta_sat:

  volumetric content at saturation

- theta_dry:

  volumetric content at the dry side of the wetting front

- I:

  Soil infiltration (in mm of water).

- widths:

  Width of soil layers (in mm).

- macro:

  Macroporosity of soil layers (in %).

- a, b:

  Parameters of the extinction function used for water infiltration.

- rainfallInput:

  Water from the rainfall event reaching the soil surface (mm)

- rainfallIntensity:

  rainfall intensity rate (mm/h)

- soil:

  A list containing the description of the soil (see
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)).

- soilFunctions:

  Soil water retention curve and conductivity functions, either 'SX'
  (for Saxton) or 'VG' (for Van Genuchten).

- model:

  Infiltration model, either "GreenAmpt1911" or "Boughton1989"

- K_correction:

  Correction for saturated conductivity, to account for increased
  infiltration due to macropore presence

## Value

Functions `hydrology_infiltrationBoughton`,
`hydrology_infiltrationGreenAmpt` and `hydrology_infiltrationAmount`
return the daily amount of water that infiltrates into the soil (in mm
of water).

Function `hydrology_infiltrationRepartition` returns the amount of
infiltrated water that reaches each soil layer.

## Details

When using function `hydrology_infiltrationGreenAmpt`, the units of
`Ksat`, `t` and `psi_wat` have to be in the same system (e.g. cm/h, h
and cm).

## References

Boughton (1989). A review of the USDA SCS curve number method. -
Australian Journal of Soil Research 27: 511-523.

Green, W.H. and Ampt, G.A. (1911) Studies on Soil Physics, 1: The Flow
of Air and Water through Soils. The Journal of Agricultural Science, 4,
1-24.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`hydrology_waterInputs`](https://emf-creaf.github.io/medfate/reference/hydrology_verticalInputs.md)

## Author

Miquel De Cáceres Ainsa, CREAF
