# Modify parameters

Routines to modify species parameter table or model input objects

## Usage

``` r
modifySpParams(SpParams, customParams, subsetSpecies = TRUE)

modifyCohortParams(x, customParams, verbose = TRUE)

modifyInputParams(x, customParams, verbose = TRUE)
```

## Arguments

- SpParams:

  A species parameter data frame, typically
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md).

- customParams:

  A data frame or a named vector with new parameter values (see
  details).

- subsetSpecies:

  A logical flag to indicate that the output data frame should include
  only those species mentioned in `customParams`.

- x:

  A model input object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  or
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- verbose:

  A logical flag to indicate that messages should be printed on the
  console.

## Value

Function `modifySpParams` returns a modified species parameter data
frame.

Functions `modifyCohortParams` and `modifyInputParams` return a modified
[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
or
[`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
object. Note that modifications may affect other parameters beyond those
indicated in `customParams`, as a result of parameter dependencies (see
examples).

## Details

When calling function `modifySpParams`, `customParams` should be a data
frame with as many rows as species and as many columns as parameters to
modify, plus a column called 'Name' or 'Species' to match species names
between the two tables. In both cases, the function will match input
strings with column 'Name' of `x`. Alternatively, `customParams` can
contain a column 'SpIndex' for matching of species indices, but this is
deprecated.

When calling `modifyCohortParams`, `customParams` can be a data frame
with as many rows as cohorts and as many columns as parameters to
modify, plus a column called 'Cohort' which will be matched with the
cohort names given by
[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
or
[`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).
Alternatively, `customParams` can be a named list or named numeric
vector as for `modifyInputParams`.

When calling `modifyInputParams`, `customParams` must be either a named
list or a named numeric vector. Cohort parameters are specified using
the syntax "\[cohortName\]/\[paramName\]" for names (e.g. "T2_176/Z50"
to modify parameter 'Z50' of cohort 'T2_176'). Soil layer parameters are
specified using the syntax "\[paramName\]@#layer" for names, where
\#layer is the layer index (e.g. "rfc@1" will modify the rock fragment
content of soil layer 1). Control parameters are specified using either
"\[paramName\]" (e.g "phloemConductanceFactor") or
"\[paramName\]\$\[subParamName\]" (e.g
"maximumRelativeGrowthRates\$leaf"). It may seem unnecessary to modify
soil or control parameters via a function, but `modifyInputParams` is
called from optimization functions (see
[`optimization`](https://emf-creaf.github.io/medfate/reference/optimization.md)).

## See also

[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md),
[`optimization`](https://emf-creaf.github.io/medfate/reference/optimization.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

#Initialize control parameters
control <- defaultControl("Granier")

#Initialize input
x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

# Cohort name for Pinus halepensis
PH_coh <- paste0("T1_", SpParamsMED$SpIndex[SpParamsMED$Name=="Pinus halepensis"])
PH_coh 
#> [1] "T1_148"

# Modify Z50 and Z95 of Pinus halepensis cohort 
customParams <- c(200,2000)
names(customParams) <- paste0(PH_coh,c("/Z50", "/Z95"))
x1m <- modifyInputParams(x1, customParams)

# Inspect original and modified objects 
x1$below
#>        Z50  Z95 Z100
#> T1_148 100  600   NA
#> T2_168 300 1000   NA
#> S1_165 200 1000   NA
x1m$below
#>        Z50  Z95 Z100 fineRootBiomass coarseRootSoilVolume
#> T1_148 200 2000   NA              NA                    0
#> T2_168 300 1000   NA              NA                    0
#> S1_165 200 1000   NA              NA                    0

# Inspect dependencies: fine root distribution across soil layers
x1$belowLayers$V
#>                1         2          3           4
#> T1_148 0.8604899 0.1194556 0.01511005 0.004944476
#> T2_168 0.5008953 0.4505941 0.04064831 0.007862284
#> S1_165 0.6799879 0.2737911 0.03567632 0.010544678
x1m$belowLayers$V
#>                1         2          3           4
#> T1_148 0.6402830 0.2655064 0.06472163 0.029488953
#> T2_168 0.5008953 0.4505941 0.04064831 0.007862284
#> S1_165 0.6799879 0.2737911 0.03567632 0.010544678

# Modify rock fragment content and sand proportion of soil layer 1
x1s <- modifyInputParams(x1, c("rfc@1" = 5, "sand@1" = 10))

# Inspect original and modified soils 
x1$soil
#>   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha
#> 1    300   25   25 Silt loam NA       NA 1.5  25 0.0485 5401.471 89.16112
#> 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112
#> 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112
#> 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112
#>       VG_n VG_theta_res VG_theta_sat W Temp
#> 1 1.303861        0.041     0.423715 1   NA
#> 2 1.303861        0.041     0.423715 1   NA
#> 3 1.303861        0.041     0.423715 1   NA
#> 4 1.303861        0.041     0.423715 1   NA
x1s$soil
#>   widths sand clay      usda om nitrogen  bd rfc  macro     Ksat VG_alpha
#> 1    300   10   25 Silt loam NA       NA 1.5   5 0.0167 7046.523 97.10141
#> 2    700   25   25 Silt loam NA       NA 1.5  45 0.0485 5401.471 89.16112
#> 3   1000   25   25 Silt loam NA       NA 1.5  75 0.0485 5401.471 89.16112
#> 4   2000   25   25 Silt loam NA       NA 1.5  95 0.0485 5401.471 89.16112
#>       VG_n VG_theta_res VG_theta_sat W Temp
#> 1 1.240345        0.041     0.426520 1   NA
#> 2 1.303861        0.041     0.423715 1   NA
#> 3 1.303861        0.041     0.423715 1   NA
#> 4 1.303861        0.041     0.423715 1   NA

# When modifying growth input objects dependencies increase
x1 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
customParams <- c(2000,2)
names(customParams) <- paste0(PH_coh,c("/Al2As", "/LAI_live"))
x1m <- modifyInputParams(x1, customParams)
```
