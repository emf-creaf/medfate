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
#>                w  cover hbc htc habc hatc       delta        rhob     rhop
#> canopy 0.4718550 100.00 2.8 6.9  2.4  7.9 4.826614951  0.09776106 607.3280
#> shrub  0.0170672   3.75 0.0 0.1  0.1  0.8 0.642625347  0.02655855 628.1644
#> herb   0.0000000   0.00 0.0  NA  0.0   NA 0.000000000  0.00000000 400.0000
#> woody  0.1740125     NA 0.0  NA  0.0   NA 0.006583900 26.43000000 730.0000
#> litter 0.1811697     NA 0.0  NA  0.0   NA 0.008875105 20.41324407 366.9009
#>                  PV         beta   betarel etabetarel     sigma        pDead
#> canopy 8.157094e-04 1.690024e-04 0.1129683  0.2742772  5272.152 0.0004066583
#> shrub  2.716995e-05 4.227961e-05 0.1704513  0.3907228  4141.000 0.0006800000
#> herb   0.000000e+00 0.000000e+00 1.6654521  0.8561108 11483.000 0.0000000000
#> woody  2.383732e-04 3.620548e-02 1.6654521  0.8561108  1601.050 1.0000000000
#> litter 4.937838e-04 5.563695e-02 9.0747984  0.1479229  7313.522 1.0000000000
#>              FAI        h           RV   MinFMC    MaxFMC ActFMC
#> canopy 4.4655429 21041.57 7.769360e-04 55.65721 117.13407     NA
#> shrub  0.1125108 20000.00 2.716995e-05 62.68714  98.39177     NA
#> herb   0.0000000 18608.00 0.000000e+00       NA        NA     NA
#> woody  0.3816475 18608.00 2.383732e-04       NA        NA     NA
#> litter 3.6112983 18608.00 4.937838e-04       NA        NA     NA
  
#Calculate fire behavior according to FCCS
fire_FCCS(fccs)
#> $SurfaceFire
#> $SurfaceFire$`midflame_WindSpeed [m/s]`
#> [1] 2.291788
#> 
#> $SurfaceFire$phi_wind
#> [1] 24.00006
#> 
#> $SurfaceFire$phi_slope
#> [1] 0
#> 
#> $SurfaceFire$`I_R_surf [kJ/m2/min]`
#> [1] 9944.641
#> 
#> $SurfaceFire$`I_R_litter [kJ/m2/min]`
#> [1] 1960.213
#> 
#> $SurfaceFire$`q_surf [kJ/m2]`
#> [1] 14977.29
#> 
#> $SurfaceFire$`q_litter [kJ/m2]`
#> [1] 1755.889
#> 
#> $SurfaceFire$xi_surf
#> [1] 0.09156883
#> 
#> $SurfaceFire$xi_litter
#> [1] 0.1690924
#> 
#> $SurfaceFire$`ROS_surf [m/min]`
#> [1] 1.521091
#> 
#> $SurfaceFire$`ROS_litter [m/min]`
#> [1] 4.722611
#> 
#> $SurfaceFire$`ROS_windslopecap [m/min]`
#> [1] 137.5073
#> 
#> $SurfaceFire$`ROS [m/min]`
#> [1] 4.722611
#> 
#> $SurfaceFire$`I_b [kW/m]`
#> [1] 344.7192
#> 
#> $SurfaceFire$`t_r [s]`
#> [1] 26.40498
#> 
#> $SurfaceFire$`FL [m]`
#> [1] 1.13868
#> 
#> 
#> $CrownFire
#> $CrownFire$`I_R_canopy [kJ/m2/min]`
#> [1] 18298.27
#> 
#> $CrownFire$`I_R_crown [kJ/m2/min]`
#> [1] 28242.92
#> 
#> $CrownFire$`q_canopy [kJ/m2]`
#> [1] 310.8646
#> 
#> $CrownFire$`q_crown [kJ/m2]`
#> [1] 15288.16
#> 
#> $CrownFire$xi_crown
#> [1] 0.06807188
#> 
#> $CrownFire$`canopy_WindSpeed [m/s]`
#> [1] 4.96174
#> 
#> $CrownFire$WAF
#> [1] 2.018978
#> 
#> $CrownFire$`ROS_crown [m/min]`
#> [1] 12.49535
#> 
#> $CrownFire$`I_b_crown [kW/m]`
#> [1] 1406.52
#> 
#> $CrownFire$`t_r_crown [s]`
#> [1] 14.3377
#> 
#> $CrownFire$Ic_ratio
#> [1] 0.5269629
#> 
#> $CrownFire$`FL_crown [m]`
#> [1] 5.624276
#> 
#> 
#> $FirePotentials
#> $FirePotentials$RP
#> [1] 2.368194
#> 
#> $FirePotentials$SP
#> [1] 9
#> 
#> $FirePotentials$FP
#> [1] 4.832072
#> 
#> $FirePotentials$SFP
#> [1] 9
#> 
#> $FirePotentials$IC
#> [1] 3.518974
#> 
#> $FirePotentials$TC
#> [1] 9
#> 
#> $FirePotentials$RC
#> [1] 6.402753
#> 
#> $FirePotentials$CFP
#> [1] 4.746695
#> 
#> 
  
 
```
