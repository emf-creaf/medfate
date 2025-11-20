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
#> canopy 0.52550038 100.00 2.7 7.1  2.6  7.9 4.791658510  0.10966983 592.0044
#> shrub  0.01407945   3.75 0.0 0.1  0.1  0.8 0.642625347  0.02190927 412.0091
#> herb   0.01929299  10.00 0.0  NA  0.0   NA 0.200000000  0.09646495 400.0000
#> woody  0.16542073     NA 0.0  NA  0.0   NA 0.006258824 26.43000000 438.9106
#> litter 0.23060466     NA 0.0  NA  0.0   NA 0.011699765 19.71019565 370.9679
#>                  PV         beta   betarel etabetarel     sigma        pDead
#> canopy 9.181138e-04 1.916067e-04 0.1276082  0.3053187  5284.915 0.0004081897
#> shrub  3.417267e-05 5.317666e-05 0.2856939  0.5836066  4141.000 0.1448400000
#> herb   4.823248e-05 2.411624e-04 0.6924824  0.9418071 11483.000 0.0000000000
#> woody  3.768894e-04 6.021728e-02 0.6924824  0.9418071  1601.050 1.0000000000
#> litter 6.216297e-04 5.313181e-02 9.1968815  0.1441747  7401.336 1.0000000000
#>              FAI        h           RV   MinFMC    MaxFMC ActFMC
#> canopy 5.0076821 21059.75 8.876630e-04 75.21455 113.45355     NA
#> shrub  0.1415090 20117.67 3.417267e-05 63.64891  96.53441     NA
#> herb   0.5538535 18608.00 4.823248e-05       NA        NA     NA
#> woody  0.6034187 18608.00 3.768894e-04       NA        NA     NA
#> litter 4.6008905 18608.00 6.216297e-04       NA        NA     NA
  
#Calculate fire behavior according to FCCS
fire_FCCS(fccs)
#> $SurfaceFire
#> $SurfaceFire$`midflame_WindSpeed [m/s]`
#> [1] 2.232265
#> 
#> $SurfaceFire$phi_wind
#> [1] 17.58488
#> 
#> $SurfaceFire$phi_slope
#> [1] 0
#> 
#> $SurfaceFire$`I_R_surf [kJ/m2/min]`
#> [1] 15797.37
#> 
#> $SurfaceFire$`I_R_litter [kJ/m2/min]`
#> [1] 2431.865
#> 
#> $SurfaceFire$`q_surf [kJ/m2]`
#> [1] 16391
#> 
#> $SurfaceFire$`q_litter [kJ/m2]`
#> [1] 1652.455
#> 
#> $SurfaceFire$xi_surf
#> [1] 0.1008317
#> 
#> $SurfaceFire$xi_litter
#> [1] 0.1628295
#> 
#> $SurfaceFire$`ROS_surf [m/min]`
#> [1] 1.80737
#> 
#> $SurfaceFire$`ROS_litter [m/min]`
#> [1] 4.456703
#> 
#> $SurfaceFire$`ROS_windslopecap [m/min]`
#> [1] 133.9359
#> 
#> $SurfaceFire$`ROS [m/min]`
#> [1] 4.456703
#> 
#> $SurfaceFire$`I_b [kW/m]`
#> [1] 385.2773
#> 
#> $SurfaceFire$`t_r [s]`
#> [1] 19.68641
#> 
#> $SurfaceFire$`FL [m]`
#> [1] 1.19846
#> 
#> 
#> $CrownFire
#> $CrownFire$`I_R_canopy [kJ/m2/min]`
#> [1] 21279.7
#> 
#> $CrownFire$`I_R_crown [kJ/m2/min]`
#> [1] 37077.07
#> 
#> $CrownFire$`q_canopy [kJ/m2]`
#> [1] 341.4614
#> 
#> $CrownFire$`q_crown [kJ/m2]`
#> [1] 16732.46
#> 
#> $CrownFire$xi_crown
#> [1] 0.07654697
#> 
#> $CrownFire$`canopy_WindSpeed [m/s]`
#> [1] 5.006563
#> 
#> $CrownFire$WAF
#> [1] 2.027292
#> 
#> $CrownFire$`ROS_crown [m/min]`
#> [1] 16.86238
#> 
#> $CrownFire$`I_b_crown [kW/m]`
#> [1] 2485.778
#> 
#> $CrownFire$`t_r_crown [s]`
#> [1] 14.30308
#> 
#> $CrownFire$Ic_ratio
#> [1] 0.6232663
#> 
#> $CrownFire$`FL_crown [m]`
#> [1] 6.375402
#> 
#> 
#> $FirePotentials
#> $FirePotentials$RP
#> [1] 2.984799
#> 
#> $FirePotentials$SP
#> [1] 9
#> 
#> $FirePotentials$FP
#> [1] 4.957289
#> 
#> $FirePotentials$SFP
#> [1] 9
#> 
#> $FirePotentials$IC
#> [1] 3.639106
#> 
#> $FirePotentials$TC
#> [1] 9
#> 
#> $FirePotentials$RC
#> [1] 7.437928
#> 
#> $FirePotentials$CFP
#> [1] 5.202635
#> 
#> 
  
 
```
