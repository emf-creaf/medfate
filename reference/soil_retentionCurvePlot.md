# Soil water retention and conductivity plots

Functions to display water retention curves and conductivity curves.

## Usage

``` r
soil_retentionCurvePlot(
  soil,
  model = "SX",
  layer = 1,
  psi = seq(0, -6, by = -0.01),
  relative = TRUE,
  to = "SAT"
)

soil_conductivityCurvePlot(
  soil,
  model = "SX",
  layer = 1,
  psi = seq(0, -6, by = -0.01),
  relative = TRUE,
  to = "SAT",
  log = TRUE,
  mmol = TRUE
)
```

## Arguments

- soil:

  Initialized soil object (returned by function
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)).

- model:

  model Either 'SX' or 'VG' for Saxton's or Van Genuchten's water
  retention models; or 'both' to plot both retention models.

- layer:

  Soil layer to be plotted.

- psi:

  A numeric vector specifying a sequence of water potential values.

- relative:

  Boolean flag to indicate that retention curve should be relative to
  field capacity or saturation.

- to:

  Either 'SAT' (saturation) or 'FC' (field capacity).

- log:

  Boolean to display the y-axis in logarithm units

- mmol:

  Boolean flag to indicate that saturated conductivity units should be
  returned in mmol/m/s/MPa. If `mmol = FALSE` then units are cm/day.

## Value

An object of class ggplot.

## Details

- `soil_retentionCurvePlot()` allows plotting the water retention curve
  of a given soil layer.

- `soil_conductivityCurvePlot()` allows plotting the conductivity curve
  of a given soil layer.

## See also

[`soil_texture`](https://emf-creaf.github.io/medfate/reference/soil_texture.md)

## Author

Miquel De Cáceres Ainsa, CREAF
