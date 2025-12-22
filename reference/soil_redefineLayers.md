# Redefine soil layer widths

Allows redefining soil layer widths of an input data frame of soil
parameters.

## Usage

``` r
soil_redefineLayers(x, widths = c(300, 700, 1000, 2000))
```

## Arguments

- x:

  A data frame of soil parameters (see an example in
  [`defaultSoilParams`](https://emf-creaf.github.io/medfate/reference/defaultSoilParams.md))
  or an object of class
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md).

- widths:

  A numeric vector indicating the desired layer widths, in mm.

## Value

A data frame or
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md) object
with soil parameters, depending on the class of `x`.

## Details

If an initialized
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md) is
supplied, its hydraulic parameters will be recalculated and the value of
state variables will be lost.

## See also

[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
[`defaultSoilParams`](https://emf-creaf.github.io/medfate/reference/defaultSoilParams.md)

## Author

Víctor Granda, EMF-CREAF

Miquel De Cáceres Ainsa, EMF-CREAF

## Examples

``` r
# Define initial soil with 5 layers
spar <- defaultSoilParams(5)
spar
#>   widths clay sand om nitrogen ph  bd rfc
#> 1    300   25   25 NA       NA NA 1.5  25
#> 2    700   25   25 NA       NA NA 1.5  45
#> 3   1000   25   25 NA       NA NA 1.5  75
#> 4   2000   25   25 NA       NA NA 1.5  95
#> 5   4000   25   25 NA       NA NA 1.5  98

# Redefine to four layers
soil_redefineLayers(spar)
#>   widths clay sand om nitrogen ph  bd rfc
#> 1    300   25   25 NA       NA NA 1.5  25
#> 2    700   25   25 NA       NA NA 1.5  45
#> 3   1000   25   25 NA       NA NA 1.5  75
#> 4   2000   25   25 NA       NA NA 1.5  95

# Same but after soil parameter initialization
examplesoil <- soil(spar)
examplesoil
#>   widths sand clay      usda om nitrogen ph  bd rfc  macro     Ksat VG_alpha
#> 1    300   25   25 Silt loam NA       NA NA 1.5  25 0.0485 5401.471 89.16112
#> 2    700   25   25 Silt loam NA       NA NA 1.5  45 0.0485 5401.471 89.16112
#> 3   1000   25   25 Silt loam NA       NA NA 1.5  75 0.0485 5401.471 89.16112
#> 4   2000   25   25 Silt loam NA       NA NA 1.5  95 0.0485 5401.471 89.16112
#> 5   4000   25   25 Silt loam NA       NA NA 1.5  98 0.0485 5401.471 89.16112
#>       VG_n VG_theta_res VG_theta_sat W Temp
#> 1 1.303861        0.041     0.423715 1   NA
#> 2 1.303861        0.041     0.423715 1   NA
#> 3 1.303861        0.041     0.423715 1   NA
#> 4 1.303861        0.041     0.423715 1   NA
#> 5 1.303861        0.041     0.423715 1   NA

soil_redefineLayers(examplesoil)
#>   widths sand clay      usda om nitrogen ph  bd rfc  macro     Ksat VG_alpha
#> 1    300   25   25 Silt loam NA       NA NA 1.5  25 0.0485 5401.471 89.16112
#> 2    700   25   25 Silt loam NA       NA NA 1.5  45 0.0485 5401.471 89.16112
#> 3   1000   25   25 Silt loam NA       NA NA 1.5  75 0.0485 5401.471 89.16112
#> 4   2000   25   25 Silt loam NA       NA NA 1.5  95 0.0485 5401.471 89.16112
#>       VG_n VG_theta_res VG_theta_sat W Temp
#> 1 1.303861        0.041     0.423715 1   NA
#> 2 1.303861        0.041     0.423715 1   NA
#> 3 1.303861        0.041     0.423715 1   NA
#> 4 1.303861        0.041     0.423715 1   NA
```
