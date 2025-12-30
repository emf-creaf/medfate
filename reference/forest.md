# Description of a forest stand.

`exampleforest` is an example of forest stand description, whereas
`exampleforest2` is an alternative forest description where leaf area
index and crown ratio are supplied instead of structural (density, DBH
and cover) parameters.

## Format

An object of class `forest` contains the description of the woody (tree
or shrub) cohorts of a forest patch. It has the following structure (see
additional columns in details):

- `treeData`: A data frame of tree cohorts (in rows) and the following
  columns:

  - `Species`: String with species (taxon) name or a non-negative
    integer for tree species identity (i.e., 0,1,2,...) matching
    SpParams.

  - `Height`: Total tree height (in cm).

  - `DBH`: Tree diameter at breast height (in cm).

  - `N`: Density (number of individuals/hectare) that the measured tree
    represents.

  - `Z50`: Depth (in mm) corresponding to 50% of fine roots.

  - `Z95`: Depth (in mm) corresponding to 95% of fine roots.

- `shrubData`: A data frame of shrub cohorts (in rows) and the following
  columns:

  - `Species`: String with species (taxon) name or a non-negative
    integer for shrub species identity (i.e., 0,1,2,...) matching
    SpParams.

  - `Height`: Average total height of plants (in cm).

  - `Cover`: Percent cover.

  - `Z50`: Depth (in mm) corresponding to 50% of fine roots.

  - `Z95`: Depth (in mm) corresponding to 95% of fine roots.

Additionally, it can contain information about other elements: the
herbaceous layer, seedling/sapling bank, seed bank, snags, litter or
soil carbon:

- `herbData`: A data frame of herbaceous cohorts (in rows) and the
  following columns (see an alternative way of defining the herbaceous
  layer in details):

  - `Species`: String with species (taxon) name or a non-negative
    integer for herbaceous species identity (i.e., 0,1,2,...) matching
    SpParams.

  - `Height`: Average total height of plants (in cm).

  - `Cover`: Percent cover.

  - `Z50`: Depth (in mm) corresponding to 50% of fine roots.

  - `Z95`: Depth (in mm) corresponding to 95% of fine roots.

- `seedlingBank`: An optional data frame containing seedling/sapling
  information with the following columns:

  - `Species`: String with species (taxon) name or a non-negative
    integer for seedling/sapling species identity (i.e., 0,1,2,...)
    matching SpParams.

  - `Percent`: Amount of seedling in relation to full seedling bank (in
    %).

  - `Age`: A numeric vector indicating age of seedlings/saplings in
    years.

  - `Z50`: Depth (in mm) corresponding to 50% of fine roots.

  - `Z95`: Depth (in mm) corresponding to 95% of fine roots.

- `seedBank`: An optional data frame containing seed bank information
  with the following columns:

  - `Species`: String with species (taxon) name or a non-negative
    integer for seed species identity (i.e., 0,1,2,...) matching
    SpParams.

  - `Percent`: Amount of seeds in relation to full seed bank (in %).

- `snagData`: An optional data frame containing information for snag
  cohorts, with the following columns:

  - `Species`: String with species (taxon) name.

  - `DBH`: Diameter at breast height (in cm), for tree snags only.

  - `Height`: Height of the snag cohort (tree or shrub) (in cm).

  - `DeathAge`: Years passed since plant death (in yrs).

  - `SmallBranches`: Mass of (standing) small dead branches (in g C/m2).

  - `LargeWood`: Mass of (standing) large dead wood (in g C/m2).

- `litterData`: An optional data frame containing litter information
  with the following columns:

  - `Species`: String with species (taxon) name.

  - `Leaves`: Mass of leaf litter (in g C/m2)

  - `Twigs`: Mass of twig litter (in g C/m2)

  - `SmallBranches`: Mass of (downed) small dead branches (in g C/m2)

  - `LargeWood`: Mass of (downed) large dead wood (in g C/m2)

  - `CoarseRoots`: Mass of dead coarse roots (in g C/m2)

  - `FineRoots`: Mass of fine root litter (in g C/m2)

- `SOCData`: An optional named numeric vector containing CENTURY soil
  organic carbon pools (normally from previous simulations):

  - `SurfaceMetabolic`: Metabolic carbon in the surface (in g C/m2).

  - `SoilMetabolic`: Metabolic carbon in the soil (in g C/m2)

  - `SurfaceActive`: Active decomposition carbon pool in the surface (in
    g C/m2)

  - `SoilActive`: Active decomposition carbon pool in the soil (in g
    C/m2)

  - `SurfaceSlow`: Slow decomposition carbon pool in the surface (in g
    C/m2)

  - `SoilSlow`: Slow decomposition carbon pool in the soil (in g C/m2)

  - `SoilPassive`: Passive decomposition carbon pool in the soil (in g
    C/m2)

## Source

DGCN (2005). Tercer Inventario Forestal Nacional (1997-2007): Catalunya.
Dirección General de Conservación de la Naturaleza, Ministerio de Medio
Ambiente, Madrid.

## Details

The structure presented above for `forest` objects corresponds to the
main (required and optional) data elements and their required columns. A
`forest` object can contain additional information when this is
available. Data frames `treeData`, `shrubData` and `herbData` can
contain additional columns (`Z100` is also commonly used for data frame
`seedlingBank`):

- `Z100`: Depth (in mm) corresponding to 100% of fine roots (to specify
  a truncated root distribution).

- `Age`: A numeric vector indicating age of cohorts in years. Used to
  track cohort age in simulations with
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

- `LAI`: Leaf area index (m2/m2).

- `FoliarBiomass`: Standing dry biomass of leaves (kg/m2).

- `FuelLoading`: Fine fuel loading (kg/m2).

- `CrownRatio`: The ratio between crown length and total height (between
  0 and 1)

- `ObsID`: A string identifying plant cohorts at the stage of forest
  sampling. Used to track the fate of particular plant cohorts through
  simulations. In
  [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md),
  the use of these identifiers can be combined with the control option
  `keepCohortsWithObsID` so that these cohorts are not merged or removed
  during simulations.

Leaf area index (`LAI`), foliar biomass (`FoliarBiomass`), fuel loading
(`FuelLoading`) and `CrownRatio` are used to override allometry-based
estimates of those variables when initializing inputs for functions
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
but they will have no effect in simulations of
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
[`growth_day`](https://emf-creaf.github.io/medfate/reference/growth_day.md)
or [`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).
Note that leaf area index, foliar biomass and fuel loading are related
entities, and they are treated as such in medfate. Therefore, users are
expected to supply one or the other, and not all of them at the same
time.

Instead of defining `herbData`, the herbaceous layer can be defined
collectively (avoiding the need to specify species parameters for
herbaceous species), using the following list items:

- `herbCover`: Percent cover of the herb layer.

- `herbHeight`: Mean height (in cm) of the herb layer.

- `herbLAI`: Leaf area index (m2/m2) of the herb layer.

- `herbFoliarBiomass`: Dry biomass (kg/m2) of the herb layer.

- `herbFuelLoading`: Fine fuel loading (kg/m2) of the herb layer.

## See also

[`summary.forest`](https://emf-creaf.github.io/medfate/reference/summary.forest.md),
[`emptyforest`](https://emf-creaf.github.io/medfate/reference/emptyforest.md),
[`plot.forest`](https://emf-creaf.github.io/medfate/reference/plot.forest.md),
[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)

## Examples

``` r
data(exampleforest)
exampleforest
#> $treeData
#>            Species   DBH Height   N Z50  Z95
#> 1 Pinus halepensis 37.55    800 168 100  300
#> 2     Quercus ilex 14.60    660 384 300 1000
#> 
#> $shrubData
#>             Species Height Cover Z50  Z95
#> 1 Quercus coccifera     80  3.75 200 1000
#> 
#> attr(,"class")
#> [1] "forest" "list"  

data(exampleforest2)
exampleforest2
#> $treeData
#>            Species DBH Height  N Z50  Z95 LAI CrownRatio
#> 1 Pinus halepensis  NA    800 NA 100  300 0.8       0.66
#> 2     Quercus ilex  NA    660 NA 300 1000 0.5       0.60
#> 
#> $shrubData
#>             Species Height Cover Z50  Z95  LAI CrownRatio
#> 1 Quercus coccifera     80    NA 200 1000 0.03        0.8
#> 
#> $herbLAI
#> [1] 0.25
#> 
#> attr(,"class")
#> [1] "forest" "list"  
```
