# Fire behaviour functions

Function `fire_FCCS()` implements a modification of the fire behavior
models described for the Fuel Characteristics Classification System
(FCCS) in Prichard et al. (2013). Function `fire_Rothermel()` implements
Rothermel's (1972) fire behaviour model (modified from package
'Rothermel' (Giorgio Vacchiano, Davide Ascoli)).

## Usage

``` r
fire_FCCS(
  FCCSpropsSI,
  MliveSI = as.numeric(c(90, 90, 60)),
  MdeadSI = as.numeric(c(6, 6, 6, 6, 6)),
  slope = 0,
  windSpeedSI = 11
)

fire_Rothermel(
  modeltype,
  wSI,
  sSI,
  delta,
  mx_dead,
  hSI,
  mSI,
  u,
  windDir,
  slope,
  aspect
)
```

## Arguments

- FCCSpropsSI:

  A data frame describing the properties of five fuel strata (canopy,
  shrub, herbs, dead woody and litter) returned by
  [`fuel_FCCS`](https://emf-creaf.github.io/medfate/reference/fuel_properties.md).

- MliveSI:

  Moisture of live fuels (in percent of dry weight) for canopy, shrub,
  and herb strata. Live moisture values are drawn from column `ActFCM`
  in `FCCSpropsSI` if available (see
  [`fuel_FCCS`](https://emf-creaf.github.io/medfate/reference/fuel_properties.md)).
  Otherwise, moisture values supplied for `MliveSI` are used.

- MdeadSI:

  Moisture of dead fuels (in percent of dry weight) for canopy, shrub,
  herb, woody and litter strata.

- slope:

  Slope (in degrees).

- windSpeedSI:

  Wind speed (in m/s) at 20 ft (6 m) over vegetation (default 11 m/s =
  40 km/h)

- modeltype:

  'S'(tatic) or 'D'(ynamic)

- wSI:

  A vector of fuel load (t/ha) for five fuel classes.

- sSI:

  A vector of surface-to-volume ratio (m2/m3) for five fuel classes.

- delta:

  A value of fuel bed depth (cm).

- mx_dead:

  A value of dead fuel moisture of extinction (percent).

- hSI:

  A vector of heat content (kJ/kg) for five fuel classes.

- mSI:

  A vector of percent moisture on a dry weight basis (percent) for five
  fuel classes.

- u:

  A value of windspeed (m/s) at midflame height.

- windDir:

  Wind direction (in degrees from north). North means blowing from north
  to south.

- aspect:

  Aspect (in degrees from north).

## Value

Both functions return list with fire behavior variables.

In the case of `fire_FCCS`, the function returns the variables in three
blocks (lists `SurfaceFire`, `CrownFire` and `FirePotentials`), and the
values are:

- `` SurfaceFire$`midflame_WindSpeed [m/s]` ``: Midflame wind speed in
  the surface fire.

- `SurfaceFire$phi_wind`: Spread rate modifier due to wind.

- `SurfaceFire$phi_slope`: Spread rate modifier due to slope.

- `` SurfaceFire$`I_R_surf [kJ/m2/min]` ``: Intensity of the surface
  fire reaction.

- `` SurfaceFire$`I_R_litter [kJ/m2/min]` ``: Intensity of the litter
  fire reaction.

- `` SurfaceFire$`q_surf [kJ/m2]` ``: Heat sink of the surface fire.

- `` SurfaceFire$`q_litter [kJ/m2]` ``: Heat sink of the litter fire.

- `SurfaceFire$xi_surf`: Propagating flux ratio of the surface fire.

- `SurfaceFire$xi_litter`: Propagating flux ratio of the litter fire.

- `` SurfaceFire$`ROS_surf [m/min]` ``: Spread rate of the surface
  fire(without accounting for faster spread in the litter layer).

- `` SurfaceFire$`ROS_litter [m/min]` ``: Spread rate of the litter
  fire.

- `` SurfaceFire$`ROS_windslopecap [m/min]` ``: Maximum surface fire
  spread rate according to wind speed.

- `` SurfaceFire$`ROS [m/min]` ``: Final spread rate of the surface
  fire.

- `` SurfaceFire$`I_b [kW/m]` ``: Fireline intensity of the surface
  fire.

- `` SurfaceFire$`FL [m]` ``: Flame length of the surface fire.

- `` CrownFire$`I_R_canopy [kJ/m2/min]` ``: Intensity of the canopy fire
  reaction.

- `` CrownFire$`I_R_crown [kJ/m2/min]` ``: Intensity of the crown fire
  reaction (adding surface and canopy reactions).

- `` CrownFire$`q_canopy [kJ/m2]` ``: Heat sink of the canopy fire.

- `` CrownFire$`q_crown [kJ/m2]` ``: Heat sink of the crown fire (adding
  surface and canopy heat sinks).

- `CrownFire$xi_surf`: Propagating flux ratio of the crown fire.

- `` CrownFire$`canopy_WindSpeed [m/s]` ``: Wind speed in the canopy
  fire (canopy top wind speed).

- `CrownFire$WAF`: Wind speed adjustment factor for crown fires.

- `` CrownFire$`ROS [m/min]` ``: Spread rate of the crown fire.

- `CrownFire$Ic_ratio`: Crown initiation ratio.

- `` CrownFire$`I_b [kW/m]` ``: Fireline intensity of the crown fire.

- `` CrownFire$`FL [m]` ``: Flame length of the crown fire.

- `FirePotentials$RP`: Surface fire reaction potential (\[0-9\]).

- `FirePotentials$SP`: Surface fire spread rate potential (\[0-9\]).

- `FirePotentials$FP`: Surface fire flame length potential (\[0-9\]).

- `FirePotentials$SFP`: Surface fire potential (\[0-9\]).

- `FirePotentials$IC`: Crown initiation potential (\[0-9\]).

- `FirePotentials$TC`: Crown-to-crown transmission potential (\[0-9\]).

- `FirePotentials$RC`: Crown fire spread rate potential (\[0-9\]).

- `FirePotentials$CFC`: Crown fire potential (\[0-9\]).

## Details

Default moisture, slope and windspeed values are benchmark conditions
used to calculate fire potentials (Sandberg et al. 2007) and map
vulnerability to fire.

## Note

Default moisture, slope and windspeed values are benchmark conditions
used to calculate fire potentials (Sandberg et al. 2007) and map
vulnerability to fire.

## References

Albini, F. A. (1976). Computer-based models of wildland fire behavior: A
users' manual. Ogden, UT: US Department of Agriculture, Forest Service,
Intermountain Forest and Range Experiment Station.

Rothermel, R. C. 1972. A mathematical model for predicting fire spread
in wildland fuels. USDA Forest Service Research Paper INT USA.

Prichard, S. J., D. V Sandberg, R. D. Ottmar, E. Eberhardt, A. Andreu,
P. Eagle, and K. Swedin. 2013. Classification System Version 3.0:
Technical Documentation.

## See also

[`fuel_FCCS`](https://emf-creaf.github.io/medfate/reference/fuel_properties.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Calculate fuel properties according to FCCS
fccs <- fuel_FCCS(exampleforest, SpParamsMED)
fccs
#>                 w  cover hbc htc habc hatc       delta        rhob     rhop
#> canopy 0.52550038 100.00 2.8 7.1  2.4  7.9 4.833189136  0.10872746 605.5066
#> shrub  0.01688118   3.75 0.0 0.1  0.1  0.8 0.642625347  0.02626908 628.1644
#> herb   0.00000000   0.00 0.0  NA  0.0   NA 0.000000000  0.00000000 400.0000
#> woody  0.19252124     NA 0.0  NA  0.0   NA 0.007284194 26.43000000 730.0000
#> litter 0.20052115     NA 0.0  NA  0.0   NA 0.009766190 20.53217836 366.9009
#>                  PV         beta   betarel etabetarel     sigma        pDead
#> canopy 9.109147e-04 1.884707e-04 0.1262340  0.3024463  5284.915 0.0004081897
#> shrub  2.687382e-05 4.181879e-05 0.2024229  0.4494103  4141.000 0.0006800000
#> herb   0.000000e+00 0.000000e+00 1.6654521  0.8561108 11483.000 0.0000000000
#> woody  2.637277e-04 3.620548e-02 1.6654521  0.8561108  1601.050 1.0000000000
#> litter 5.465268e-04 5.596111e-02 9.0541457  0.1485569  7298.666 1.0000000000
#>              FAI        h           RV   MinFMC    MaxFMC ActFMC
#> canopy 4.9967361 21059.75 8.678690e-04 55.62584 117.21661     NA
#> shrub  0.1112845 20000.00 2.687382e-05 62.68714  98.39177     NA
#> herb   0.0000000 18608.00 0.000000e+00       NA        NA     NA
#> woody  0.4222413 18608.00 2.637277e-04       NA        NA     NA
#> litter 3.9889167 18608.00 5.465268e-04       NA        NA     NA
  
#Calculate fire behavior according to FCCS
fire_FCCS(fccs)
#> $SurfaceFire
#> $SurfaceFire$`midflame_WindSpeed [m/s]`
#> [1] 2.258073
#> 
#> $SurfaceFire$phi_wind
#> [1] 21.76496
#> 
#> $SurfaceFire$phi_slope
#> [1] 0
#> 
#> $SurfaceFire$`I_R_surf [kJ/m2/min]`
#> [1] 10957.84
#> 
#> $SurfaceFire$`I_R_litter [kJ/m2/min]`
#> [1] 2178.891
#> 
#> $SurfaceFire$`q_surf [kJ/m2]`
#> [1] 15003.88
#> 
#> $SurfaceFire$`q_litter [kJ/m2]`
#> [1] 1773.69
#> 
#> $SurfaceFire$xi_surf
#> [1] 0.101426
#> 
#> $SurfaceFire$xi_litter
#> [1] 0.1699028
#> 
#> $SurfaceFire$`ROS_surf [m/min]`
#> [1] 1.687518
#> 
#> $SurfaceFire$`ROS_litter [m/min]`
#> [1] 4.754843
#> 
#> $SurfaceFire$`ROS_windslopecap [m/min]`
#> [1] 135.4844
#> 
#> $SurfaceFire$`ROS [m/min]`
#> [1] 4.754843
#> 
#> $SurfaceFire$`I_b [kW/m]`
#> [1] 382.867
#> 
#> $SurfaceFire$`t_r [s]`
#> [1] 26.43495
#> 
#> $SurfaceFire$`FL [m]`
#> [1] 1.195005
#> 
#> 
#> $CrownFire
#> $CrownFire$`I_R_canopy [kJ/m2/min]`
#> [1] 21202.28
#> 
#> $CrownFire$`I_R_crown [kJ/m2/min]`
#> [1] 32160.12
#> 
#> $CrownFire$`q_canopy [kJ/m2]`
#> [1] 345.4915
#> 
#> $CrownFire$`q_crown [kJ/m2]`
#> [1] 15349.37
#> 
#> $CrownFire$xi_crown
#> [1] 0.07575536
#> 
#> $CrownFire$`canopy_WindSpeed [m/s]`
#> [1] 5.006563
#> 
#> $CrownFire$WAF
#> [1] 2.027292
#> 
#> $CrownFire$`ROS_crown [m/min]`
#> [1] 14.30609
#> 
#> $CrownFire$`I_b_crown [kW/m]`
#> [1] 1829.264
#> 
#> $CrownFire$`t_r_crown [s]`
#> [1] 14.30308
#> 
#> $CrownFire$Ic_ratio
#> [1] 0.5852784
#> 
#> $CrownFire$`FL_crown [m]`
#> [1] 6.003662
#> 
#> 
#> $FirePotentials
#> $FirePotentials$RP
#> [1] 2.485908
#> 
#> $FirePotentials$SP
#> [1] 9
#> 
#> $FirePotentials$FP
#> [1] 4.950139
#> 
#> $FirePotentials$SFP
#> [1] 9
#> 
#> $FirePotentials$IC
#> [1] 3.593623
#> 
#> $FirePotentials$TC
#> [1] 9
#> 
#> $FirePotentials$RC
#> [1] 6.850984
#> 
#> $FirePotentials$CFP
#> [1] 4.96184
#> 
#> 
  
 
```
