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
#> 1   0.5 0.8550725 0.01611031 0.003070609 0.2077187 -0.003156593 1.056
#> 2   1.0 0.8623221 0.01611031 0.003843443 0.2159537 -0.003156593 1.056
#> 3   1.5 0.8738023 0.02441644 0.004755993 0.2297953 -0.004942526 1.056
#> 4   2.0 0.8850883 0.02401959 0.005692790 0.2493185 -0.005070773 1.056
#> 5   2.5 0.8961739 0.02361695 0.006639946 0.2755342 -0.005246501 1.056
#> 6   3.0 0.9070512 0.02320348 0.007547697 0.3086430 -0.005459465 1.056
#> 7   3.5 0.9349807 0.05684244 0.009551333 0.3575494 -0.014396237 1.056
#> 8   4.0 0.9965002 0.12321936 0.014064599 0.4349114 -0.034388602 1.056
#> 9   4.5 1.1024788 0.21128611 0.022773360 0.5553624 -0.066520792 1.056
#> 10  5.0 1.2634226 0.32048042 0.038565340 0.7400501 -0.116225228 1.056
#> 11  5.5 1.4877301 0.44675159 0.065766621 1.0150209 -0.189347190 1.056
#> 12  6.0 1.7796965 0.58198655 0.109530957 1.4054661 -0.289757506 1.056
#> 13  6.5 2.1380296 0.71498586 0.173799260 1.9255898 -0.416170718 1.056
#> 14  7.0 2.5557528 0.83425304 0.258167007 2.5681284 -0.560379831 1.056
#> 15  7.5 3.0008802 0.88979336 0.346410777 3.2493215 -0.672082116 1.056
#> 16  8.0 3.4676143 0.93344426 0.430499881 3.9433012 -0.776599671 1.056
#> 17  8.5 3.9520287 0.96904099 0.499662705 4.6222636 -1.038149018 1.256
#> 18  9.0 4.3415186 0.77968801 0.523301827 5.0622764 -1.013427511 1.456
#> 19  9.5 4.6691390 0.65607082 0.532374145 5.3806174 -1.000000000 1.656
#> 20 10.0 4.9525037 0.56759594 0.538793103 5.6312500 -1.000000000 1.856
```
