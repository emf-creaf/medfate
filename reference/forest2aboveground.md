# Input for simulation models (deprecated)

Functions `forest2spwbInput()` and `forest2growthInput()` take an object
of class
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md) and
a soil data input to create input objects for simulation functions
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) (or
[`pwb`](https://emf-creaf.github.io/medfate/reference/pwb.md)) and
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
respectively. Function `forest2aboveground()` calculates aboveground
variables such as leaf area index. Function `forest2belowground()`
calculates belowground variables such as fine root distribution.

## Usage

``` r
forest2aboveground(x, SpParams, gdd = NA_real_, loading = FALSE)

forest2belowground(x, soil, SpParams)

forest2spwbInput(x, soil, SpParams, control)

forest2growthInput(x, soil, SpParams, control)
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md).

- SpParams:

  A data frame with species parameters (see
  [`SpParamsDefinition`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  and
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)).

- gdd:

  Growth degree days to account for leaf phenology effects (in Celsius).
  This should be left `NA` in most applications.

- loading:

  A logical flag to indicate that fuel loading should be included (for
  fire hazard calculations).

- soil:

  An object of class
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
  containing soil parameters per soil layer.

- control:

  A list with default control parameters (see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

## Value

Function `forest2aboveground()` returns a data frame with the following
columns (rows are identified as specified by function
[`plant_ID`](https://emf-creaf.github.io/medfate/reference/plant_values.md)):

- `SP`: Species identity (an integer) (first species is 0).

- `N`: Cohort density (ind/ha) (see function
  [`plant_density`](https://emf-creaf.github.io/medfate/reference/plant_values.md)).

- `DBH`: Tree diameter at breast height (cm).

- `H`: Plant total height (cm).

- `CR`: Crown ratio (crown length to total height) (between 0 and 1).

- `LAI_live`: Live leaf area index (m2/m2) (one-side leaf area relative
  to plot area), includes leaves in winter dormant buds.

- `LAI_expanded`: Leaf area index of expanded leaves (m2/m2) (one-side
  leaf area relative to plot area).

- `LAI_dead`: Dead leaf area index (m2/m2) (one-side leaf area relative
  to plot area).

- `LAI_nocomp`: Leaf area index (m2/m2) (one-side leaf area relative to
  plot area) assuming no aboveground competition.

- `Loading`: Fine fuel loading (kg/m2), only if `loading = TRUE`.

- `Age`: A numeric vector indicating age of cohorts in years. Used to
  track cohort age in simulations with
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

- `ObsID`: A string identifying plant cohorts at the stage of forest
  sampling. Used to track the fate of particular plant cohorts in
  simulations with
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

## Details

Function `forest2aboveground()` extracts height and species identity
from plant cohorts of `x`, and calculate leaf area index and crown
ratio.

*IMPORTANT NOTE*: Function names `forest2spwbInput()` and
`forest2growthInput()` are now internal and deprecated, but they can
still be used for back-compatibility. They correspond to functions
[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
and
[`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)

## See also

[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md),
[`defaultSoilParams`](https://emf-creaf.github.io/medfate/reference/defaultSoilParams.md),
[`plant_ID`](https://emf-creaf.github.io/medfate/reference/plant_values.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

# Aboveground parameters
forest2aboveground(exampleforest, SpParamsMED)
#>         SP        N   DBH Cover   H        CR   LAI_live LAI_expanded LAI_dead
#> T1_148 148 168.0000 37.55    NA 800 0.6605196 0.84874773   0.84874773        0
#> T2_168 168 384.0000 14.60    NA 660 0.6055642 0.70557382   0.70557382        0
#> S1_165 165 749.4923    NA  3.75  80 0.8032817 0.03062604   0.03062604        0
#>        LAI_nocomp Age ObsID
#> T1_148 1.29720268  NA  <NA>
#> T2_168 1.01943205  NA  <NA>
#> S1_165 0.04412896  NA  <NA>

# Example of aboveground parameters taken from a forest
# described using LAI and crown ratio
data(exampleforest2)
forest2aboveground(exampleforest2, SpParamsMED)
#>         SP  N DBH Cover   H   CR LAI_live LAI_expanded LAI_dead LAI_nocomp Age
#> T1_148 148 NA  NA    NA 800 0.66     0.80         0.80        0       0.80  NA
#> T2_168 168 NA  NA    NA 660 0.60     0.50         0.50        0       0.50  NA
#> S1_165 165 NA  NA    NA  80 0.80     0.03         0.03        0       0.03  NA
#>        ObsID
#> T1_148  <NA>
#> T2_168  <NA>
#> S1_165  <NA>

# Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

# Bewowground parameters (distribution of fine roots)
forest2belowground(exampleforest, examplesoil, SpParamsMED)
#>                1         2          3           4
#> T1_148 0.8604899 0.1194556 0.01511005 0.004944476
#> T2_168 0.5008953 0.4505941 0.04064831 0.007862284
#> S1_165 0.6799879 0.2737911 0.03567632 0.010544678

```
