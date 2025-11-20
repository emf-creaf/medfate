# Data tables with species parameter definitions and values

A data sets of species parameter definition and values, the latter
resulting from existing databases, fit to empirical data or expert-based
guesses.

## Format

- Data frame `SpParamsDefinition` has parameters in rows and columns
  'ParameterName', 'ParameterGroup', 'Definition', 'Type' and 'Units'.

- Data frames `SpParamsMED` has species or genus as rows and column
  names equal to parameter names in `SpParamsDefinition`.

## Details

`SpParamsMED` was the official species parameter for package versions up
to v.4.0.0, but will not be maintained in the future. Additional species
parameter tables for different countries are distributed via package
[traits4models](https://emf-creaf.github.io/traits4models/).

## Examples

``` r
data(SpParamsDefinition)
data(SpParamsMED)
```
