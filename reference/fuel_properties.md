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
#> [1] 280
#> 
#> $canopyTopHeight
#> [1] 760
#> 
#> $canopyAbsoluteBaseHeight
#> [1] 240
#> 
#> $canopyAbsoluteTopHeight
#> [1] 790
#> 
#> $canopyLAI
#> [1] 1.695727
#> 
  
#Calculate fuel properties according to FCCS
fccs <- fuel_FCCS(exampleforest, SpParamsMED)
fccs
#>                 w  cover hbc htc habc hatc       delta        rhob     rhop
#> canopy 0.62413684 100.00 2.8 7.6  2.4  7.9 4.960217784  0.12582851 570.3129
#> shrub  0.01561312   3.75 0.0 0.1  0.1  0.8 0.642625347  0.02429584 628.1644
#> herb   0.00000000   0.00 0.0  NA  0.0   NA 0.000000000  0.00000000 400.0000
#> woody  0.22731400     NA 0.0  NA  0.0   NA 0.008600605 26.43000000 730.0000
#> litter 0.23084538     NA 0.0  NA  0.0   NA 0.010374497 22.25123541 366.9009
#>                  PV         beta   betarel etabetarel     sigma        pDead
#> canopy 1.138507e-03 2.295276e-04 0.1593580  0.3693686  5531.511 0.0004377814
#> shrub  2.485515e-05 3.867752e-05 0.2788513  0.5735398  4141.000 0.0006800000
#> herb   0.000000e+00 0.000000e+00 1.6654521  0.8561108 11483.000 0.0000000000
#> woody  3.113890e-04 3.620548e-02 1.6654521  0.8561108  1601.050 1.0000000000
#> litter 6.291765e-04 6.064645e-02 8.7556347  0.1577217  7083.948 1.0000000000
#>              FAI        h           RV   MinFMC    MaxFMC ActFMC
#> canopy 6.4740133 21411.15 1.094376e-03 55.01978 118.81135     NA
#> shrub  0.1029252 20000.00 2.485515e-05 62.68714  98.39177     NA
#> herb   0.0000000 18608.00 0.000000e+00       NA        NA     NA
#> woody  0.4985494 18608.00 3.113890e-04       NA        NA     NA
#> litter 4.4570531 18608.00 6.291765e-04       NA        NA     NA

```
