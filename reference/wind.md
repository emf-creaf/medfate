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
#> 1   0.5 0.8976867 0.01790035 0.003555201 0.2358795 -0.003841021 1.056
#> 2   1.0 0.9057419 0.01790035 0.004430662 0.2450295 -0.003841021 1.056
#> 3   1.5 0.9214767 0.03298998 0.005536175 0.2616747 -0.007326270 1.056
#> 4   2.0 0.9368416 0.03223347 0.006658577 0.2840095 -0.007466124 1.056
#> 5   2.5 0.9518399 0.03149592 0.007792155 0.3132268 -0.007668203 1.056
#> 6   3.0 0.9686140 0.03498894 0.009034259 0.3507736 -0.009020231 1.056
#> 7   3.5 1.0042684 0.07224693 0.011620628 0.4067173 -0.020055793 1.056
#> 8   4.0 1.0733309 0.13828607 0.016975507 0.4933194 -0.042242799 1.056
#> 9   4.5 1.1857708 0.22423342 0.026852532 0.6247334 -0.076963013 1.056
#> 10  5.0 1.3509392 0.32900058 0.044080646 0.8205695 -0.129166186 1.056
#> 11  5.5 1.5761026 0.44855879 0.072653086 1.1039961 -0.203877374 1.056
#> 12  6.0 1.8648340 0.57559926 0.117047950 1.4963973 -0.304103719 1.056
#> 13  6.5 2.2158190 0.70032639 0.180381672 2.0087127 -0.428195837 1.056
#> 14  7.0 2.6227064 0.81255924 0.261755803 2.6327623 -0.568369762 1.056
#> 15  7.5 3.0561575 0.86636314 0.346315586 3.2927366 -0.677487705 1.056
#> 16  8.0 3.5098584 0.90728710 0.426441922 3.9633296 -0.778275536 1.056
#> 17  8.5 3.9793793 0.93916414 0.492375633 4.6186505 -1.034345477 1.256
#> 18  9.0 4.3579571 0.75781934 0.517034233 5.0527930 -1.012062132 1.456
#> 19  9.5 4.6767184 0.63834791 0.528812661 5.3735399 -1.000000000 1.656
#> 20 10.0 4.9525037 0.55246219 0.538793103 5.6312500 -1.000000000 1.856
```
