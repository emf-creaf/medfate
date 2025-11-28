# Soil initialization

Initializes soil parameters and state variables for its use in
simulations.

## Usage

``` r
soil(x, VG_PTF = "Toth")

# S3 method for class 'soil'
summary(object, model = "SX", ...)
```

## Arguments

- x:

  A data frame of soil parameters (see an example in
  [`defaultSoilParams`](https://emf-creaf.github.io/medfate/reference/defaultSoilParams.md)).

- VG_PTF:

  Pedotransfer functions to obtain parameters for the van
  Genuchten-Mualem equations. Either `"Carsel"` (Carsel and
  Parrish 1988) or `"Toth"` (Toth et al. 2015).

- object:

  An object of class `soil`.

- model:

  Either 'SX' or 'VG' for Saxton or Van Genuchten pedotransfer models.

- ...:

  Additional parameters to `summary`.

## Value

Function `soil` returns a data frame of class `soil` with the following
columns:

- `widths`: Width of soil layers (in mm).

- `sand`: Sand percentage for each layer (in percent volume).

- `clay`: Clay percentage for each layer (in percent volume).

- `om`: Organic matter percentage for each layer (in percent volume).

- `nitrogen`: Sum of total nitrogen (ammonia, organic and reduced
  nitrogen) for each layer (in g/kg).

- `rfc`: Percentage of rock fragment content for each layer.

- `macro`: Macroporosity for each layer (estimated using Stolf et al.
  2011).

- `Ksat`: Saturated soil conductivity for each layer (in
  mmol·m-1·s-1·MPa-1, estimated using function
  [`soil_saturatedConductivitySX`](https://emf-creaf.github.io/medfate/reference/soil_texture.md).

- `VG_alpha`, `VG_n`, `VG_theta_res`, `VG_theta_sat`: Parameters for van
  Genuchten's pedotransfer functions, for each layer, corresponding to
  the USDA texture type.

- `W`: State variable with relative water content of each layer (in as
  proportion relative to FC).

- `Temp`: State variable with temperature (in ºC) of each layer.

## Details

Function `summary` prompts a description of soil characteristics and
state variables (water content and temperature) according to a water
retention curve (either Saxton's or Van Genuchten's). Volume at field
capacity is calculated assuming a soil water potential equal to -0.033
MPa. Parameter `Temp` is initialized as missing for all soil layers.

If available, the user can specify columns `VG_alpha`, `VG_n`,
`VG_theta_res`, `VG_theta_sat` and `K_sat`, to override Van Genuchten
parameters an saturated conductivity estimated from pedotransfer
functions when calling function `soil`.

## References

Carsel, R.F., and Parrish, R.S. 1988. Developing joint probability
distributions of soil water retention characteristics. Water Resources
Research 24: 755–769.

Tóth, B., Weynants, M., Nemes, A., Makó, A., Bilas, G., and Tóth, G.
2015. New generation of hydraulic pedotransfer functions for Europe.
European Journal of Soil Science 66: 226–238.

Stolf, R., Thurler, A., Oliveira, O., Bacchi, S., Reichardt, K., 2011.
Method to estimate soil macroporosity and microporosity based on sand
content and bulk density. Rev. Bras. Ciencias do Solo 35, 447–459.

## See also

[`soil_redefineLayers`](https://emf-creaf.github.io/medfate/reference/soil_redefineLayers.md),
[`soil_psi2thetaSX`](https://emf-creaf.github.io/medfate/reference/soil_texture.md),
[`soil_psi2thetaVG`](https://emf-creaf.github.io/medfate/reference/soil_texture.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`defaultSoilParams`](https://emf-creaf.github.io/medfate/reference/defaultSoilParams.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
# Default parameters
df_soil <- defaultSoilParams()

# Initializes soil
s = soil(df_soil)
s
#>   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha     VG_n VG_theta_res
#> 1    300   25   25 Silt loam NA       NA 1.5  25 0.0485 5401.471 89.16112 1.303861        0.041
#> 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112 1.303861        0.041
#> 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112 1.303861        0.041
#> 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112 1.303861        0.041
#>   VG_theta_sat W Temp
#> 1     0.423715 1   NA
#> 2     0.423715 1   NA
#> 3     0.423715 1   NA
#> 4     0.423715 1   NA

# Prints soil characteristics according to Saxton's water retention curve
summary(s, model="SX")
#> Soil depth (mm): 4000 
#> 
#> Layer  1  [ 0  to  300 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 25 Macroporosity (%): 5 
#>     Theta WP (%): 14 Theta FC (%): 30 Theta SAT (%): 49 Theta current (%) 30 
#>     Vol. WP (mm): 32 Vol. FC (mm): 68 Vol. SAT (mm): 111 Vol. current (mm): 68 
#>     Temperature (Celsius): NA 
#> 
#> Layer  2  [ 300  to  1000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 45 Macroporosity (%): 5 
#>     Theta WP (%): 14 Theta FC (%): 30 Theta SAT (%): 49 Theta current (%) 30 
#>     Vol. WP (mm): 55 Vol. FC (mm): 117 Vol. SAT (mm): 190 Vol. current (mm): 117 
#>     Temperature (Celsius): NA 
#> 
#> Layer  3  [ 1000  to  2000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 75 Macroporosity (%): 5 
#>     Theta WP (%): 14 Theta FC (%): 30 Theta SAT (%): 49 Theta current (%) 30 
#>     Vol. WP (mm): 36 Vol. FC (mm): 76 Vol. SAT (mm): 123 Vol. current (mm): 76 
#>     Temperature (Celsius): NA 
#> 
#> Layer  4  [ 2000  to  4000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 95 Macroporosity (%): 5 
#>     Theta WP (%): 14 Theta FC (%): 30 Theta SAT (%): 49 Theta current (%) 30 
#>     Vol. WP (mm): 14 Vol. FC (mm): 30 Vol. SAT (mm): 49 Vol. current (mm): 30 
#>     Temperature (Celsius): NA 
#> 
#> Total soil saturated capacity (mm): 473 
#> Total soil water holding capacity (mm): 291 
#> Total soil extractable water (mm): 183 
#> Total soil current Volume (mm): 291 
#> Saturated water depth (mm): NA 
#> 

# Prints soil characteristics according to Van Genuchten's water retention curve
summary(s, model="VG")
#> Soil depth (mm): 4000 
#> 
#> Layer  1  [ 0  to  300 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 25 Macroporosity (%): 5 
#>     Theta WP (%): 13 Theta FC (%): 30 Theta SAT (%): 42 Theta current (%) 30 
#>     Vol. WP (mm): 29 Vol. FC (mm): 68 Vol. SAT (mm): 95 Vol. current (mm): 68 
#>     Temperature (Celsius): NA 
#> 
#> Layer  2  [ 300  to  1000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 45 Macroporosity (%): 5 
#>     Theta WP (%): 13 Theta FC (%): 30 Theta SAT (%): 42 Theta current (%) 30 
#>     Vol. WP (mm): 49 Vol. FC (mm): 117 Vol. SAT (mm): 163 Vol. current (mm): 117 
#>     Temperature (Celsius): NA 
#> 
#> Layer  3  [ 1000  to  2000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 75 Macroporosity (%): 5 
#>     Theta WP (%): 13 Theta FC (%): 30 Theta SAT (%): 42 Theta current (%) 30 
#>     Vol. WP (mm): 32 Vol. FC (mm): 76 Vol. SAT (mm): 106 Vol. current (mm): 76 
#>     Temperature (Celsius): NA 
#> 
#> Layer  4  [ 2000  to  4000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 95 Macroporosity (%): 5 
#>     Theta WP (%): 13 Theta FC (%): 30 Theta SAT (%): 42 Theta current (%) 30 
#>     Vol. WP (mm): 13 Vol. FC (mm): 30 Vol. SAT (mm): 42 Vol. current (mm): 30 
#>     Temperature (Celsius): NA 
#> 
#> Total soil saturated capacity (mm): 407 
#> Total soil water holding capacity (mm): 291 
#> Total soil extractable water (mm): 194 
#> Total soil current Volume (mm): 291 
#> Saturated water depth (mm): NA 
#> 

# Add columns 'VG_theta_sat' and 'VG_theta_res' with custom values
df_soil$VG_theta_sat <- 0.400 
df_soil$VG_theta_res <- 0.040 

# Reinitialize soil (should override estimations)
s2 = soil(df_soil)
s2
#>   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha     VG_n VG_theta_res
#> 1    300   25   25 Silt loam NA       NA 1.5  25 0.0485 5401.471 89.16112 1.303861         0.04
#> 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112 1.303861         0.04
#> 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112 1.303861         0.04
#> 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112 1.303861         0.04
#>   VG_theta_sat W Temp
#> 1          0.4 1   NA
#> 2          0.4 1   NA
#> 3          0.4 1   NA
#> 4          0.4 1   NA
summary(s2, model="VG")
#> Soil depth (mm): 4000 
#> 
#> Layer  1  [ 0  to  300 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 25 Macroporosity (%): 5 
#>     Theta WP (%): 12 Theta FC (%): 29 Theta SAT (%): 40 Theta current (%) 29 
#>     Vol. WP (mm): 27 Vol. FC (mm): 64 Vol. SAT (mm): 90 Vol. current (mm): 64 
#>     Temperature (Celsius): NA 
#> 
#> Layer  2  [ 300  to  1000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 45 Macroporosity (%): 5 
#>     Theta WP (%): 12 Theta FC (%): 29 Theta SAT (%): 40 Theta current (%) 29 
#>     Vol. WP (mm): 47 Vol. FC (mm): 110 Vol. SAT (mm): 154 Vol. current (mm): 110 
#>     Temperature (Celsius): NA 
#> 
#> Layer  3  [ 1000  to  2000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 75 Macroporosity (%): 5 
#>     Theta WP (%): 12 Theta FC (%): 29 Theta SAT (%): 40 Theta current (%) 29 
#>     Vol. WP (mm): 30 Vol. FC (mm): 72 Vol. SAT (mm): 100 Vol. current (mm): 72 
#>     Temperature (Celsius): NA 
#> 
#> Layer  4  [ 2000  to  4000 mm ] 
#>     clay (%): 25 silt (%): 50 sand (%): 25 organic matter (%): NA [ Silt loam ]
#>     Rock fragment content (%): 95 Macroporosity (%): 5 
#>     Theta WP (%): 12 Theta FC (%): 29 Theta SAT (%): 40 Theta current (%) 29 
#>     Vol. WP (mm): 12 Vol. FC (mm): 29 Vol. SAT (mm): 40 Vol. current (mm): 29 
#>     Temperature (Celsius): NA 
#> 
#> Total soil saturated capacity (mm): 384 
#> Total soil water holding capacity (mm): 275 
#> Total soil extractable water (mm): 182 
#> Total soil current Volume (mm): 275 
#> Saturated water depth (mm): NA 
#> 
```
