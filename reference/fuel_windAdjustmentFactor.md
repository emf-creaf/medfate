# Wind adjustment factor for Rothermel's model

Function fuel_windAdjustmentFactor determines the adjustment factor of
wind for surface fires, according to Andrews (2012).

## Usage

``` r
fuel_windAdjustmentFactor(
  topShrubHeight,
  bottomCanopyHeight,
  topCanopyHeight,
  canopyCover
)
```

## Arguments

- topShrubHeight:

  Shrub stratum top height (in m).

- bottomCanopyHeight:

  Canopy base height (in m).

- topCanopyHeight:

  Canopy top height (in m).

- canopyCover:

  Canopy percent cover.

## Value

A scalar value between 0 and 1

## References

Andrews, P. L. 2012. Modeling wind adjustment factor and midflame wind
speed for Rothermel’s surface fire spread model. USDA Forest Service -
General Technical Report RMRS-GTR:1–39.

## Examples

``` r
#Load example plot plant data
 data(exampleforest)
  
#Default species parameterization
data(SpParamsMED)

#Calculate fuel properties according to FCCS
fccs <- fuel_FCCS(exampleforest, SpParamsMED)

# Estimate wind adjustment factor
fuel_windAdjustmentFactor(fccs$htc[2], fccs$hbc[1], fccs$htc[1], fccs$cover[1])
#> [1] 0.2029332
```
