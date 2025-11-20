# Example observed data

Example (fake) data set of variables measured in a plot.

## Format

A data frame containing daily 'observed' values for year 2001:

- `dates`:

  Measurement dates.

- `SWC`:

  Soil moisture content (in m3/m3).

- `ETR`:

  Total evapotranspiration (mm).

- `E_T1_148`:

  Transpiration of Pinus halepensis cohort 'T1_148' (L/m2 of leaf area).

- `E_T2_168`:

  Transpiration of Quercus ilex cohort 'T2_168' (L/m2 of leaf area).

- `FMC_T1_148`:

  Fuel moisture content of Pinus halepensis cohort 'T1_148' (in
  percent).

- `FMC_T2_168`:

  Fuel moisture content of Quercus ilex cohort 'T2_168' (in percent).

- `BAI_T1_148`:

  Basal area increment for Pinus halepensis cohort 'T1_148' (in cm2).

- `BAI_T2_168`:

  Basal area increment for Quercus ilex cohort 'T2_168' (in cm2).

- `DI_T1_148`:

  Diameter increment for Pinus halepensis cohort 'T1_148' (in cm).

- `DI_T2_168`:

  Diameter increment for Quercus ilex cohort 'T2_168' (in cm).

## Source

This data set was actually created by running a simulation and adding
some gaussian error to the outputs.

## See also

[`evaluation`](https://emf-creaf.github.io/medfate/reference/evaluation.md)

## Examples

``` r
data(exampleobs)
```
