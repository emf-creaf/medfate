# Example daily meteorology data

Example data set of meteorological input.

## Format

A data frame containing daily meteorology of a location in Catalonia
(Spain) for year 2001:

- `dates`:

  Vector of [`Date`](https://rdrr.io/r/base/Dates.html) objects.

- `MinTemperature`:

  Minimum daily temperature (in degrees Celsius).

- `MaxTemperature`:

  Maximum daily temperature (in degrees Celsius).

- `Precipitation`:

  Daily precipitation (in mm of water).

- `MinRelativeHumidity`:

  Minimum daily relative humidity (in percent).

- `MaxRelativeHumidity`:

  Maximum daily relative humidity (in percent).

- `Radiation`:

  Incoming radiation (in MJ/m2).

- `WindSpeed`:

  Wind speed (in m/s).

## Source

Interpolated from weather station data (Spanish and Catalan meteorology
agencies) using package 'meteoland'.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)

## Examples

``` r
 data(examplemeteo)
```
