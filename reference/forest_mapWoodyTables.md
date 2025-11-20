# Map forest plot data

Mapping functions to facilitate building forest objects from forest plot
data

## Usage

``` r
forest_mapTreeTable(x, mapping_x, SpParams, plot_size_x = NULL)

forest_mapShrubTable(y, mapping_y, SpParams, plot_size_y = NULL)

forest_mapWoodyTables(
  x = NULL,
  y = NULL,
  mapping_x = NULL,
  mapping_y = NULL,
  SpParams,
  plot_size_x = NULL,
  plot_size_y = NULL
)
```

## Arguments

- x:

  A data frame with tree records in rows and attributes in columns. Tree
  records can correspond to individual trees or groups of trees with an
  associated density.

- mapping_x:

  A named character vector to specify mappings of columns in `x` into
  attributes of `treeData` data frames. Accepted names (and the
  corresponding specifications for the columns in `x`) are:

- SpParams:

  A data frame with species parameters (see
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md))
  from which valid species names are drawn.

- plot_size_x:

  The size of tree plot sampled area (in m2). Alternatively,
  'plot_size_x' can be a column in `x` and specified in `mapping_x` to
  indicate that trees have been measured in different subplots and,
  therefore, they represent different densities per hectare.

- y:

  A data frame with shrub records in rows and attributes in columns.
  Records can correspond to individual shrubs (with crown dimensions and
  height) or groups of shrubs with an associated cover estimate.

- mapping_y:

  A named character vector to specify mappings of columns in `y` into
  attributes of `shrubData` data frames. Accepted names (and the
  corresponding specifications for the columns in `y`) are:

  - `Species`: Species code (should follow codes in `SpParams`).

  - `Species.name`: Species name. In this case, the species code will be
    drawn by matching names with species names in `SpParams`.

  - `N`: Tree density (in individuals per ha).

  - `Cover`: Shrub cover (in percent).

  - `D1`: Shrub largest crown diameter (in cm).

  - `D2`: Shrub crown diameter orthogonal to the largest one (in cm).

  - `plot.size`: Plot size (in m2) to which each record refers to. This
    is used to calculate tree density (stems per hectare) when not
    supplied or shrub cover when shrub data is given at the individual
    level.

  - `DBH`: Diameter at breast height (in cm).

  - `Height`: Tree or shrub height (in cm).

  - `Z50`: Depth (in mm) corresponding to 50 percent of fine roots.

  - `Z95`: Depth (in mm) corresponding to 95 percent of fine roots.

  - `Z100`: Depth (in mm) corresponding to 100 percent of fine roots.

- plot_size_y:

  The size of shrub plot sampled area (in m2). Alternatively,
  'plot_size_y' can be a column in `y` and specified in `mapping_y` to
  indicate that shrubs have been measured in different subplots and,
  therefore, they represent different cover values.

## Value

Functions `forest_mapTreeTable` and `forest_mapShrubTable` return a data
frame with the structure of `treeData` and `shrubData` from
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
objects. Function `forest_mapWoodyTable` returns directly a
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
object.

## See also

[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`poblet_trees`](https://emf-creaf.github.io/medfate/reference/poblet_trees.md),
[`forest_mergeTrees`](https://emf-creaf.github.io/medfate/reference/forest_simplification.md),
[`tree2forest`](https://emf-creaf.github.io/medfate/reference/tree2forest.md)

## Author

Miquel De Cáceres Ainsa, EMF-CREAF

## Examples

``` r
# Load species parameters
data(SpParamsMED)

# Create an empty forest object
f <- emptyforest()

# (1) Mapping tree data
# Load Poblet tree data
data(poblet_trees)

# Subset control plot
x <- subset(poblet_trees, Plot.Code=="POBL_CTL")

# Estimate sampled area (15-m radius plot)
sampled_area <- pi*15^2

# Define mapping
mapping_x <- c("Species.name" = "Species", "DBH" = "Diameter.cm")

# Map tree data for plot 'POBL_CTL'
f$treeData <- forest_mapTreeTable(x,
                    mapping_x = mapping_x, SpParams = SpParamsMED,
                    plot_size_x = sampled_area)

# (2) Mapping shrub individual data
#
# Create the individual shrub data frame
species <- c("Erica arborea","Cistus albidus", "Erica arborea", "Cistus albidus", "Cistus albidus")
H <- c(200,50,100,40,30)
D1 <- c(140,40,100, 35,30)
D2 <- D1
y <- data.frame(species, H, D1, D2)

# Define mapping (D1 and D2 map to variables with the same name)
mapping_y <- c("Species.name"= "species", "Height" ="H", "D1", "D2")

# Map individual shrub data to cover data (here each individual becomes a cohort)
# assuming that the sampled area was 4 m2
f$shrubData <- forest_mapShrubTable(y,
                     mapping_y = mapping_y, SpParams = SpParamsMED,
                     plot_size_y = 4)

# (3) Print forest attributes
summary(f, SpParamsMED)
#> Tree BA (m2/ha): 42.6957047  adult trees: 42.6957047  saplings: 0 
#> Density (ind/ha) adult trees: 3777.277316  saplings: 0  shrubs (estimated): 19051.5105038 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 65.4334845  herbs: 0 
#> LAI (m2/m2) total: 6.0900572  adult trees: 5.6770407  saplings: 0  shrubs: 0.4130165  herbs: 0 
#> Fuel loading (kg/m2) total: 1.5959112  adult trees: 1.493419  saplings: 0  shrubs: 0.1024922  herbs: 0 
#> PAR ground (%): NA  SWR ground (%): NA 

# (4) Forest initialization in a single step
f <- forest_mapWoodyTables(x, y,
                           mapping_x = mapping_x, mapping_y = mapping_y,
                           SpParams = SpParamsMED,
                           plot_size_x = sampled_area, plot_size_y = 4)
summary(f, SpParamsMED)
#> Tree BA (m2/ha): 42.6957047  adult trees: 42.6957047  saplings: 0 
#> Density (ind/ha) adult trees: 3777.277316  saplings: 0  shrubs (estimated): 19051.5105038 
#> Cover (%) adult trees: 100  saplings: 0  shrubs: 65.4334845  herbs: 0 
#> LAI (m2/m2) total: 6.0900572  adult trees: 5.6770407  saplings: 0  shrubs: 0.4130165  herbs: 0 
#> Fuel loading (kg/m2) total: 1.5959112  adult trees: 1.493419  saplings: 0  shrubs: 0.1024922  herbs: 0 
#> PAR ground (%): NA  SWR ground (%): NA 
```
