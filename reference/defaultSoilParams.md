# Default soil parameters

Creates a data frame with default soil physical description for model
functions

## Usage

``` r
defaultSoilParams(n = 4)
```

## Arguments

- n:

  An integer with the number of soil layers (between two and five).

## Value

A data frame with layers in rows and the following columns (and default
values):

- `widths (= c(300,700,1000,2000)`: Width of soil layers (in mm).

- `clay (= 25)`: Clay percentage for each layer (in %).

- `sand (= 25)`: Sand percentage for each layer (in %).

- `om (= NA)`: Organic matter percentage for each layer (in %)
  (optional).

- `nitrogen (= NA)`: Sum of total nitrogen (ammonia, organic and reduced
  nitrogen) for each layer (in g/kg) (optional).

- `bd (= 1.5)`: Bulk density for each layer (in g/cm3).

- `rfc (= c(20,40,60,85))`: Percentage of rock fragment content (volume
  basis) for each layer.

## Details

The function returns a data frame with default physical soil
description, with soil layers in rows. Users can change those that need
to be set to other values and use the list as input for function
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md).

## Note

While this function is limited to five soil layers, user defined data
frames can discretize soils using an unlimited number of soil layers.

## See also

[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
[`soil_redefineLayers`](https://emf-creaf.github.io/medfate/reference/soil_redefineLayers.md),
[`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md),
[`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
defaultSoilParams(4)
#>   widths clay sand om nitrogen  bd rfc
#> 1    300   25   25 NA       NA 1.5  25
#> 2    700   25   25 NA       NA 1.5  45
#> 3   1000   25   25 NA       NA 1.5  75
#> 4   2000   25   25 NA       NA 1.5  95
```
