# Fuel stratification and fuel characteristics

Function `fuel_stratification` provides a stratification of the stand
into understory and canopy strata. Function `fuel_FCCS` calculates fuel
characteristics from a `forest` object following an adaptation of the
protocols described for the Fuel Characteristics Classification System
(Prichard et al. 2013).

## Usage

``` r
fuel_stratification(
  object,
  SpParams,
  gdd = NA_real_,
  heightProfileStep = 10,
  maxHeightProfile = 5000,
  bulkDensityThreshold = 0.05
)

fuel_FCCS(
  object,
  SpParams,
  cohortFMC = as.numeric(c()),
  loadingOffset = as.numeric(c(0, 0, 0, 0, 0)),
  gdd = NA_real_,
  heightProfileStep = 10,
  maxHeightProfile = 5000,
  bulkDensityThreshold = 0.05,
  depthMode = "crownaverage"
)
```

## Arguments

- object:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- gdd:

  Growth degree-days.

- heightProfileStep:

  Precision for the fuel bulk density profile.

- maxHeightProfile:

  Maximum height for the fuel bulk density profile.

- bulkDensityThreshold:

  Minimum fuel bulk density to delimit fuel strata.

- cohortFMC:

  A numeric vector of (actual) fuel moisture content by cohort.

- loadingOffset:

  A vector of length five with fine fuel loading values (canopy, shrub,
  herb, woody and litter) to be added to loading estimations from
  `forest`.

- depthMode:

  Specifies how fuel depth (and therefore canopy and understory bulk
  density) should be estimated:

  - `"crownaverage"`: As weighed average of crown lengths using loadings
    as weights.

  - `"profile"`: As the difference of base and top heights in bulk
    density profiles.

  - `"absoluteprofile"`: As the difference of absolute base and absolute
    top heights in bulk density profiles.

## Value

Function `fuel_FCCS` returns a data frame with five rows corresponding
to fuel layers: `canopy`, `shrub`, `herb`, `woody` and `litter`. Columns
correspond fuel properties:

- `w`: Fine fuel loading (in kg/m2).

- `cover`: Percent cover.

- `hbc`: Height to base of crowns (in m).

- `htc`: Height to top of crowns (in m).

- `delta`: Fuel depth (in m).

- `rhob`: Fuel bulk density (in kg/m3).

- `rhop`: Fuel particle density (in kg/m3).

- `PV`: Particle volume (in m3/m2).

- `beta`: Packing ratio (unitless).

- `betarel`: Relative packing ratio (unitless).

- `etabetarel`: Reaction efficiency (unitless).

- `sigma`: Surface area-to-volume ratio (m2/m3).

- `pDead`: Proportion of dead fuels.

- `FAI`: Fuel area index (unitless).

- `h`: High heat content (in kJ/kg).

- `RV`: Reactive volume (in m3/m2).

- `MinFMC`: Minimum fuel moisture content (as percent over dry weight).

- `MaxFMC`: Maximum fuel moisture content (as percent over dry weight).

- `ActFMC`: Actual fuel moisture content (as percent over dry weight).
  These are set to `NA` if parameter `cohortFMC` is empty.

Function `fuel_stratification` returns a list with the following items:

- `surfaceLayerBaseHeight`: Base height of crowns of shrubs in the
  surface layer (in cm).

- `surfaceLayerTopHeight`: Top height of crowns of shrubs in the surface
  layer (in cm).

- `understoryLAI`: Cumulated LAI of the understory layer (i.e. leaf area
  comprised between surface layer base and top heights).

- `canopyBaseHeight`: Base height of tree crowns in the canopy (in cm).

- `canopyTopHeight`: Top height of tree crowns in the canopy (in cm).

- `canopyLAI`: Cumulated LAI of the canopy (i.e. leaf area comprised
  between canopy base and top heights).

## References

Prichard, S. J., D. V Sandberg, R. D. Ottmar, E. Eberhardt, A. Andreu,
P. Eagle, and K. Swedin. 2013. Classification System Version 3.0:
Technical Documentation.

Reinhardt, E., D. Lutes, and J. Scott. 2006. FuelCalc: A method for
estimating fuel characteristics. Pages 273–282.

## See also

[`fire_FCCS`](https://emf-creaf.github.io/medfate/reference/fire_behaviour.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Show stratification of fuels
fuel_stratification(exampleforest, SpParamsMED)
#> $surfaceLayerBaseHeight
#> [1] 0
#> 
#> $surfaceLayerTopHeight
#> [1] 10
#> 
#> $surfaceLayerAbsoluteBaseHeight
#> [1] 10
#> 
#> $surfaceLayerAbsoluteTopHeight
#> [1] 80
#> 
#> $understoryLAI
#> [1] 0
#> 
#> $canopyBaseHeight
#> [1] 270
#> 
#> $canopyTopHeight
#> [1] 710
#> 
#> $canopyAbsoluteBaseHeight
#> [1] 260
#> 
#> $canopyAbsoluteTopHeight
#> [1] 790
#> 
#> $canopyLAI
#> [1] 1.453648
#> 
  
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

```
