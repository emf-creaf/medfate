# Extracts model outputs

Function `extract()` extracts daily or subdaily output and returns it as
a tidy data frame.

## Usage

``` r
extract(
  x,
  level = "forest",
  output = NULL,
  vars = NULL,
  dates = NULL,
  subdaily = FALSE,
  addunits = FALSE
)
```

## Arguments

- x:

  An object returned by simulation functions
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
  [`aspwb`](https://emf-creaf.github.io/medfate/reference/aspwb.md),
  [`pwb`](https://emf-creaf.github.io/medfate/reference/pwb.md),
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md) or
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

- level:

  Level of simulation output, either "forest" (stand-level results),
  "soillayer" (soil layer-level results), "cohort" (cohort-level
  results), "sunlitleaf" or "shadeleaf" (leaf-level results)

- output:

  Section of the model output to be explored. See details.

- vars:

  Variables to be extracted (by default, all of them).

- dates:

  A date vector indicating the subset of simulated days for which output
  is desired.

- subdaily:

  A flag to indicate that subdaily values are desired (see details).

- addunits:

  A flag to indicate that variable units should be added whenever
  possible.

## Value

Function `extract()` returns a data frame:

- If `level = "forest"`, columns are "date" and variable names.

- If `level = "soillayer"`, columns are "date", "soillayer" and variable
  names.

- If `level = "cohort"`, `level = "sunlitleaf"` or
  `level = "shadeleaf"`, columns are "date", "cohorts", "species" and
  variable names.

- If `subdaily = TRUE`, columns are "datetime", "cohorts", "species" and
  variable names.

## Details

When `subdaily = FALSE`, parameter `output` is used to restrict the
section in `x` where variables are located. For example
`output = "Plants"` will correspond to variables "LAI", "LAIlive",
"Transpiration", "StemPLC",... as returned by a call `names(x$Plants)`.

Option `subdaily = TRUE` only works when simulations have been carried
using control option 'subdailyResults = TRUE' (see
[`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).
When using `subdaily = TRUE`, parameter `output` is not taken into
account, and options for parameter `vars` are the following:

- Variables for `level = "forest"` or `level = "soillayer"`: Not
  allowed. An error is raised.

- Variables for `level = "cohort"`:
  "E","Ag","An","dEdP","RootPsi","StemPsi","LeafPsi","StemPLC","StemRWC","LeafRWC","StemSympRWC","LeafSympRWC","PWB".

- Variables for `level = "shadeleaf"` and `level="sunlitleaf"`:
  "Abs_SWR","Abs_PAR","Net_LWR","E","Ag","An","Ci","Gsw","VPD","Temp","Psi","iWUE".

## See also

[`summary.spwb`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md)

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
x <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function (ten days)
S1<-spwb(x, examplemeteo[1:10, ], latitude = 41.82592, elevation = 100)
#> Initial plant water content (mm): 4.69853
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  [Year 2001]:
#> 
#> Final plant water content (mm): 4.69847
#> Final soil water content (mm): 294.834
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -6.55784e-05
#> Plant water balance result (mm): -6.55784e-05
#> Change in soil water content (mm): 3.95915
#> Soil water balance result (mm): 3.95915
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 0
#> Water balance components:
#>   Precipitation (mm) 39 Rain (mm) 34 Snow (mm) 5
#>   Interception (mm) 8 Net rainfall (mm) 26
#>   Infiltration (mm) 31 Infiltration excess (mm) 0 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 3  Herbaceous transpiration (mm) 0 Woody plant transpiration (mm) 3
#>   Plant extraction from soil (mm) 3  Plant water balance (mm) -0 Hydraulic redistribution (mm) 0
#>   Runoff (mm) 0 Deep drainage (mm) 22

#Extracts daily forest-level output as a data frame
extract(S1, level = "forest", addunits = TRUE)
#> Error in UseMethod("extract"): no applicable method for 'extract' applied to an object of class "c('spwb', 'list')"
#> Error in extract(S1, level = "forest", addunits = TRUE): Arguments in `...` must be used.
#> ✖ Problematic arguments:
#> • level = "forest"
#> • addunits = TRUE
#> ℹ Did you misspell an argument name?

#Extracts daily soil layer-level output as a data frame
extract(S1, level = "soillayer", addunits = TRUE)
#> Error in UseMethod("extract"): no applicable method for 'extract' applied to an object of class "c('spwb', 'list')"
#> Error in extract(S1, level = "soillayer", addunits = TRUE): Arguments in `...` must be used.
#> ✖ Problematic arguments:
#> • level = "soillayer"
#> • addunits = TRUE
#> ℹ Did you misspell an argument name?

#Extracts daily cohort-level output as a data frame
extract(S1, level = "cohort", addunits = TRUE)
#> Error in UseMethod("extract"): no applicable method for 'extract' applied to an object of class "c('spwb', 'list')"
#> Error in extract(S1, level = "cohort", addunits = TRUE): Arguments in `...` must be used.
#> ✖ Problematic arguments:
#> • level = "cohort"
#> • addunits = TRUE
#> ℹ Did you misspell an argument name?

#Select the output tables/variables to be extracted
extract(S1, level ="cohort", output="Plants", vars = c("PlantStress", "StemPLC"))
#> Error in UseMethod("extract"): no applicable method for 'extract' applied to an object of class "c('spwb', 'list')"
#> Error in extract(S1, level = "cohort", output = "Plants", vars = c("PlantStress",     "StemPLC")): Arguments in `...` must be used.
#> ✖ Problematic arguments:
#> • level = "cohort"
#> • output = "Plants"
#> • vars = c("PlantStress", "StemPLC")
#> ℹ Did you misspell an argument name?
```
