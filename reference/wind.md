# Models for canopy turbulence

Models for canopy turbulence by Katul et al (2004).

## Usage

``` r
wind_canopyTurbulenceModel(zm, Cx, hm, d0, z0, model = "k-epsilon")

wind_canopyTurbulence(
  zmid,
  LAD,
  canopyHeight,
  u,
  windMeasurementHeight = 200,
  model = "k-epsilon"
)
```

## Arguments

- zm:

  A numeric vector with height values (m).

- Cx:

  Effective drag = Cd x leaf area density.

- hm:

  Canopy height (m).

- d0:

  Zero displacement height (m).

- z0:

  Momentum roughness height (m).

- model:

  Closure model.

- zmid:

  A numeric vector of mid-point heights (in cm) for canopy layers.

- LAD:

  A numeric vector of leaf area density values (m3/m2).

- canopyHeight:

  Canopy height (in cm).

- u:

  Measured wind speed (m/s).

- windMeasurementHeight:

  Height of wind speed measurement with respect to canopy height (cm).

## Value

Function `wind_canopyTurbulenceModel` returns a data frame of vertical
profiles for variables:

- `z1`: Height values.

- `U1`: U/u\*, where U is mean velocity and u\* is friction velocity.

- `dU1`: dUdz/u\*, where dUdz is mean velocity gradient and u\* is
  friction velocity.

- `epsilon1`: epsilon/(u^3/h) where epsilon is the turbulent kinetic
  dissipation rate, u\* is friction velocity and h is canopy height.

- `k1`: k/(u\*^2), where k is the turbulent kinetic energy and u\* is
  friction velocity.

- `uw1`: uw/(u\*^2), where uw is the Reynolds stress and u\* is friction
  velocity.

- `Lmix1`: Mixing length.

Function `wind_canopyTurbulence` returns a data frame of vertical
profiles for transformed variables:

- `zmid`: Input mid-point heights (in cm) for canopy layers.

- `u`: Wind speed (m/s).

- `du`: Mean velocity gradient (1/s).

- `epsilon`: Turbulent kinetic dissipation rate.

- `k`: Turbulent kinetic energy.

- `uw`: Reynolds stress.

## Details

Implementation in Rcpp of the K-epsilon canopy turbulence models by
Katul et al (2004) originally in Matlab code
(https://nicholas.duke.edu/people/faculty/katul/k_epsilon_model.htm).

## References

Katul GG, Mahrt L, Poggi D, Sanz C (2004) One- and two-equation models
for canopy turbulence. Boundary-Layer Meteorol 113:81–109.
https://doi.org/10.1023/B:BOUN.0000037333.48760.e5

## See also

[`vprofile_windExtinction`](https://emf-creaf.github.io/medfate/reference/vprofile_leafAreaDensity.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Default species parameterization
data(SpParamsMED)

#Load example plot plant data
data(exampleforest)

#Canopy height (in m)
h= max(exampleforest$treeData$Height/100) 
d0 = 0.67*h
z0 = 0.08*h

#Height values (cm)
z = seq(50,1000, by=50)
zm = z/100 # (in m)

# Leaf area density
lad = vprofile_leafAreaDensity(exampleforest, SpParamsMED, draw = FALSE,
                               z = c(0,z))
  
# Effective drag
Cd = 0.2
Cx = Cd*lad
  
# canopy turbulence model
wind_canopyTurbulenceModel(zm, Cx,h,d0,z0)
#>      z1        U1        dU1    epsilon1        k1          uw1 Lmix1
#> 1   0.5 0.8771067 0.01790035 0.003310369 0.2262022 -0.003700077 1.056
#> 2   1.0 0.8851619 0.01790035 0.004185830 0.2353521 -0.003700077 1.056
#> 3   1.5 0.9005842 0.03236945 0.005264621 0.2516399 -0.006929133 1.056
#> 4   2.0 0.9156470 0.03163474 0.006341795 0.2734986 -0.007068251 1.056
#> 5   2.5 0.9303537 0.03091913 0.007407269 0.3020996 -0.007267271 1.056
#> 6   3.0 0.9467057 0.03415407 0.008544316 0.3387575 -0.008506022 1.056
#> 7   3.5 0.9812685 0.07008201 0.010940058 0.3929857 -0.018799114 1.056
#> 8   4.0 1.0486878 0.13501612 0.015991296 0.4768353 -0.039860220 1.056
#> 9   4.5 1.1588413 0.21967239 0.025382494 0.6041078 -0.072879682 1.056
#> 10  5.0 1.3211630 0.32330652 0.041874250 0.7942168 -0.122740408 1.056
#> 11  5.5 1.5432474 0.44238140 0.069450240 1.0706074 -0.194603677 1.056
#> 12  6.0 1.8292514 0.57010444 0.112737942 1.4558017 -0.291958707 1.056
#> 13  6.5 2.1786093 0.69701881 0.175244947 1.9628728 -0.413982255 1.056
#> 14  7.0 2.5856994 0.81291075 0.256638622 2.5862302 -0.553774257 1.056
#> 15  7.5 3.0227532 0.87351990 0.342893758 3.2541967 -0.667250674 1.056
#> 16  8.0 3.4832552 0.92085788 0.425900167 3.9396751 -0.773826995 1.056
#> 17  8.5 3.9621123 0.95782453 0.494742074 4.6128006 -1.035841159 1.256
#> 18  9.0 4.3475754 0.77160633 0.519606183 5.0536114 -1.012585285 1.456
#> 19  9.5 4.6719318 0.64955414 0.530416195 5.3755536 -1.000000000 1.656
#> 20 10.0 4.9525037 0.56204626 0.538793103 5.6312500 -1.000000000 1.856
```
