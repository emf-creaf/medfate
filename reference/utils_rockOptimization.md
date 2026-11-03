# Optimization of rock fragment content

Function `utils_rockOptimization` finds optimum rock fragment content in
the soil corresponding to given vegetation, weather and target percent
loss of conductance (PLC), following a modification of the method
proposed by Druel et al. (2023).

## Usage

``` r
utils_rockOptimization(
  x,
  soil,
  SpParams,
  control,
  meteo,
  PLCquantile = 0.9,
  qPLC_target = 12,
  qPLC_tol = 0.5,
  sew_min = 30,
  max_rocks = 99,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- soil:

  An object of class
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
  containing soil parameters per soil layer.

- SpParams:

  A data frame with species parameters (see
  [`SpParamsDefinition`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  and
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- control:

  A list with default control parameters (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

- meteo:

  A data frame with daily meteorological data series (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- PLCquantile:

  Maximum PLC quantile to be calculated across years.

- qPLC_target:

  Target PLC to be achieved (by default 12%).

- qPLC_tol:

  Tolerance of PLC difference to target accepted when finding solution.

- sew_min:

  Minimum soil extractable water (mm) for rock exploration.

- max_rocks:

  Maximum content in coarse fragments allowed for any soil layer.

- verbose:

  A logical value. Print the internal messages of the function?

- ...:

  Additional parameters to function
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md).

## Value

Function `utils_rockOptimization` returns a list containing:

- `RFC`: A vector with the estimated rock fragment content for each soil
  layer.

- `SEW`: Soil extractable water (mm).

- `runs`: Number of simulations performed.

- `message`: Text message indicating whether optimization could be done
  (OK) or not.

## Details

The function performs a model inversion based on an ecohydrological
assumption, consisting in that forest leaf area index is in equilibrium
with a low embolism rate under normal conditions. This is translated in
that the (by default 90%) interannual quantile of the maximum annual
percent loss of conductance (PLC), averaged over plant cohorts, should
be close to a target PLC value (by default 12%).

The algorithm first determines the PLC corresponding to the minimum and
maximum soil extractable water (SEW). The minimum SEW (SEW_min) is an
input parameter, whereas the maximum SEW (SEW_max) corresponds to no
rock fragments in the soil.

Then three situations are distinguished:

1.  If `PLC(SEW_min) < qPLC_target` and `PLC(SEW_max) < qPLC_target`,
    the function will use
    [`uniroot`](https://rdrr.io/r/stats/uniroot.html) to find the root
    of the function `f(x) = PLC(x) - qPLC_target`, where `x` is SEW,
    which corresponds to a factor that multiplies the original rock
    fragment content.

2.  If both `PLC(SEW_min) < qPLC_target` and
    `PLC(SEW_max) < qPLC_target`, the function cannot find an optimum,
    because PLC is always too low, and will return the original rock
    fragment content

3.  Analogously, if both `PLC(SEW_min) > qPLC_target` and
    `PLC(SEW_max) > qPLC_target`, the function cannot find an optimum,
    because PLC is always too large, and will return the original rock
    fragment content

## References

Druel, A., Martins, N., Cochard, H., De Caceres, M., Delzon, S.,
Mencuccini, M., Torres-Ruiz, J., and Ruffault, J.: European forest
vulnerability to hydraulic failure: an ecohydrological approach, EGU
General Assembly 2023, Vienna, Austria, 24–28 Apr 2023, EGU23-17068,
https://doi.org/10.5194/egusphere-egu23-17068, 2023.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
[`utils_ldrOptimization`](https://emf-creaf.github.io/medfate/reference/utils_ldrOptimization.md)

## Author

Arsène Druel, URFM-INRAE

Nicolas Martin-StPaul, URFM-INRAE

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
# \donttest{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with two layers
examplesoil <- defaultSoilParams(4)

#Rock fragment content optimization (Granier)
utils_rockOptimization(exampleforest, soil = examplesoil,
                       SpParams = SpParamsMED, meteo = examplemeteo,
                       control = defaultControl("Granier"),
                       elevation = 100, latitude = 41.82592)
#> $RFC
#> [1] 25 45 75 95
#> 
#> $SEW
#> [1] 193.9525
#> 
#> $runs
#> [1] 2
#> 
#> $message
#> [1] "PLC lower than target for the whole range of SEW. Returning original SEW."
#> 
# }
```
