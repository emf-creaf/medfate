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
#> 1   0.5 0.7506036 0.01611031 0.002184125 0.1697988 -0.002669998 1.056
#> 2   1.0 0.7578532 0.01611031 0.002859429 0.1780338 -0.002669998 1.056
#> 3   1.5 0.7707260 0.02713846 0.003724363 0.1914893 -0.004672562 1.056
#> 4   2.0 0.7833098 0.02655020 0.004517893 0.2098024 -0.004791237 1.056
#> 5   2.5 0.7956052 0.02597262 0.005214342 0.2339515 -0.004954432 1.056
#> 6   3.0 0.8091031 0.02833618 0.005845081 0.2646464 -0.005752844 1.056
#> 7   3.5 0.8383893 0.05946380 0.007365186 0.3090357 -0.013045314 1.056
#> 8   4.0 0.8984498 0.12024976 0.010973156 0.3780062 -0.029145794 1.056
#> 9   4.5 0.9988562 0.20013962 0.018023315 0.4839961 -0.054784431 1.056
#> 10  5.0 1.1495492 0.29999743 0.030962972 0.6458571 -0.094631728 1.056
#> 11  5.5 1.3597602 0.41854076 0.053751538 0.8887614 -0.154499332 1.056
#> 12  6.0 1.6364288 0.55128712 0.091766039 1.2410790 -0.239983876 1.056
#> 13  6.5 1.9824795 0.69022606 0.150479788 1.7265353 -0.353859041 1.056
#> 14  7.0 2.3956389 0.82488644 0.232484440 2.3525175 -0.493164134 1.056
#> 15  7.5 2.8536979 0.91540151 0.327260616 3.0632852 -0.624197144 1.056
#> 16  8.0 3.3497381 0.99188362 0.424989344 3.8257576 -0.755673218 1.056
#> 17  8.5 3.8762609 1.05319400 0.509242413 4.5922802 -1.045511476 1.256
#> 18  9.0 4.2961137 0.84049330 0.534498864 5.0647502 -1.015857598 1.456
#> 19  9.5 4.6482246 0.70512898 0.539495761 5.3894102 -1.000000000 1.656
#> 20 10.0 4.9525037 0.60947693 0.538793103 5.6312500 -1.000000000 1.856
```
