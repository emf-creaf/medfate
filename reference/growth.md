# Forest growth

Function `growth` is a process-based model that performs energy, water
and carbon balances; and determines changes in water/carbon pools,
functional variables (leaf area, sapwood area, root area) and structural
ones (tree diameter, tree height, shrub cover) for woody plant cohorts
in a given forest stand during a period specified in the input climatic
data.

## Usage

``` r
growth(
  x,
  meteo,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  CO2ByYear = numeric(0),
  waterTableDepth = NA_real_
)
```

## Arguments

- x:

  An object of class
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- meteo:

  A data frame with daily meteorological data series (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- latitude:

  Latitude (in degrees).

- elevation, slope, aspect:

  Elevation above sea level (in m), slope (in degrees) and aspect (in
  degrees from North).

- CO2ByYear:

  A named numeric vector with years as names and atmospheric CO2
  concentration (in ppm) as values. Used to specify annual changes in
  CO2 concentration along the simulation (as an alternative to
  specifying daily values in `meteo`).

- waterTableDepth:

  Water table depth (in mm). When not missing, capillarity rise will be
  allowed if lower than total soil depth.

## Value

Function `growth` returns a list of class 'growth'. Since lists are
difficult to handle, we recommend using function
[`extract`](https://emf-creaf.github.io/medfate/reference/extract.md) to
reshape simulation results (including their units) from those objects.

List elements are as follows:

- `"latitude"`: Latitude (in degrees) given as input.

- `"topography"`: Vector with elevation, slope and aspect given as
  input.

- `"weather"`: A copy of the input weather data frame.

- `"growthInput"`: A copy of the object `x` of class
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  given as input.

- `"growthOutput"`: An copy of the final state of the object `x` of
  class
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- `"WaterBalance"`: A data frame where different water balance variables
  (see [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- `"EnergyBalance"`: A data frame with the daily values of energy
  balance components for the soil and the canopy (only for
  `transpirationMode = "Sperry"` or `transpirationMode = "Sureau"`; see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- `"CarbonBalance"`: A data frame where different stand-level carbon
  balance components (gross primary production, maintenance respiration,
  synthesis respiration and net primary production), all in g C · m-2.

- `"BiomassBalance"`: A data frame with the daily values of stand
  biomass balance components (in g dry · m-2.

- `"Temperature"`: A data frame with the daily values of
  minimum/mean/maximum temperatures for the atmosphere (input), canopy
  and soil (only for `transpirationMode = "Sperry"` or
  `transpirationMode = "Sureau"`; see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- `"Soil"`: A data frame where different soil variables (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- `"Stand"`: A data frame where different stand-level variables (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- `"Plants"`: A list of daily results for plant cohorts (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- `"SunlitLeaves"` and `"ShadeLeaves"`: A list with daily results for
  sunlit and shade leaves (only for `transpirationMode = "Sperry"` or
  `transpirationMode = "Sureau"`; see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- `"LabileCarbonBalance"`: A list of daily labile carbon balance results
  for plant cohorts, with elements:

  - `"GrossPhotosynthesis"`: Daily gross photosynthesis per dry weight
    of living biomass (g gluc · g dry-1).

  - `"MaintentanceRespiration"`: Daily maintenance respiration per dry
    weight of living biomass (g gluc · g dry-1).

  - `"GrowthCosts"`: Daily growth costs per dry weight of living biomass
    (g gluc · g dry-1).

  - `"RootExudation"`: Root exudation per dry weight of living biomass
    (g gluc · g dry-1).

  - `"LabileCarbonBalance"`: Daily labile carbon balance
    (photosynthesis - maintenance respiration - growth costs - root
    exudation) per dry weight of living biomass (g gluc · g dry-1).

  - `"SugarLeaf"`: Sugar concentration (mol·l-1) in leaves.

  - `"StarchLeaf"`: Starch concentration (mol·l-1) in leaves.

  - `"SugarSapwood"`: Sugar concentration (mol·l-1) in sapwood.

  - `"StarchSapwood"`: Starch concentration (mol·l-1) in sapwood.

  - `"SugarTransport"`: Average instantaneous rate of carbon transferred
    between leaves and stem compartments via floem (mol gluc·s-1).

- `"PlantBiomassBalance"`: A list of daily plant biomass balance results
  for plant cohorts, with elements:

  - `"StructuralBiomassBalance"`: Daily structural biomass balance (g
    dry · m-2).

  - `"LabileBiomassBalance"`: Daily labile biomass balance (g dry ·
    m-2).

  - `"PlantBiomassBalance"`: Daily plant biomass balance, i.e. labile
    change + structural change (g dry · m-2).

  - `"MortalityBiomassLoss"`: Biomass loss due to mortality (g dry ·
    m-2).

  - `"CohortBiomassBalance"`: Daily cohort biomass balance (including
    mortality) (g dry · m-2).

- `"PlantStructure"`: A list of daily area and biomass values for
  compartments of plant cohorts, with elements:

  - `"LeafBiomass"`: Daily amount of leaf structural biomass (in g dry)
    for an average individual of each plant cohort.

  - `"SapwoodBiomass"`: Daily amount of sapwood structural biomass (in g
    dry) for an average individual of each plant cohort.

  - `"FineRootBiomass"`: Daily amount of fine root biomass (in g dry)
    for an average individual of each plant cohort.

  - `"LeafArea"`: Daily amount of leaf area (in m2) for an average
    individual of each plant cohort.

  - `"SapwoodArea"`: Daily amount of sapwood area (in cm2) for an
    average individual of each plant cohort.

  - `"FineRootArea"`: Daily amount of fine root area (in m2) for an
    average individual of each plant cohort.

  - `"HuberValue"`: The ratio of sapwood area to (target) leaf area (in
    cm2/m2).

  - `"RootAreaLeafArea"`: The ratio of fine root area to (target) leaf
    area (in m2/m2).

  - `"DBH"`: Diameter at breast height (in cm) for an average individual
    of each plant cohort.

  - `"Height"`: Height (in cm) for an average individual of each plant
    cohort.

- `"GrowthMortality"`: A list of daily growth and mortality rates for
  plant cohorts, with elements:

  - `"LAgrowth"`: Leaf area growth (in m2·day-1) for an average
    individual of each plant cohort.

  - `"SAgrowth"`: Sapwood area growth rate (in cm2·day-1) for an average
    individual of each plant cohort.

  - `"FRAgrowth"`: Fine root area growth (in m2·day-1) for an average
    individual of each plant cohort.

  - `"StarvationRate"`: Daily mortality rate from starvation (ind/d-1).

  - `"DessicationRate"`: Daily mortality rate from dessication
    (ind/d-1).

  - `"MortalityRate"`: Daily mortality rate (any cause) (ind/d-1).

- `"subdaily"`: A list of objects of class
  [`growth_day`](https://emf-creaf.github.io/medfate/reference/growth_day.md),
  one per day simulated (only if required in `control` parameters, see
  [`defaultControl`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)).

## Details

Detailed model description is available in the medfate book.

Forest growth simulations allow using different sub-models for bulk soil
water flows and different sub-models of transpiration and photosynthesis
(see details in
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

## References

De Cáceres M, Molowny-Horas R, Cabon A, Martínez-Vilalta J, Mencuccini
M, García-Valdés R, Nadal-Sala D, Sabaté S, Martin-StPaul N, Morin X,
D'Adamo F, Batllori E, Améztegui A (2023) MEDFATE 2.9.3: A trait-enabled
model to simulate Mediterranean forest function and dynamics at regional
scales. Geoscientific Model Development 16: 3165-3201
(https://doi.org/10.5194/gmd-16-3165-2023).

## See also

[`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`growth_day`](https://emf-creaf.github.io/medfate/reference/growth_day.md),
[`plot.growth`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md)

## Author

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
  
#Initialize control parameters
control <- defaultControl("Granier")
  
#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

#Initialize model input
x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control)

#Call simulation function
G1 <- growth(x1, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant cohort biomass (g/m2): 5041.87
#> Initial plant water content (mm): 4.69853
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  Year 2001:............
#> 
#> Final plant biomass (g/m2): 5210.49
#> Change in plant biomass (g/m2): 168.622
#> Plant biomass balance result (g/m2): 168.622
#> Plant biomass balance components:
#>   Structural balance (g/m2) 87 Labile balance (g/m2) 89
#>   Plant individual balance (g/m2) 176 Mortality loss (g/m2) 8
#> Final plant water content (mm): 4.706
#> Final soil water content (mm): 275.66
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): 0.00746396
#> Plant water balance result (mm): -0.00136322
#> Change in soil water content (mm): -15.2151
#> Soil water balance result (mm): -15.2151
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 0
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 83 Net rainfall (mm) 379
#>   Infiltration (mm) 409 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 26  Herbaceous transpiration (mm) 0 Woody plant transpiration (mm) 247
#>   Plant extraction from soil (mm) 247  Plant water balance (mm) -0 Hydraulic redistribution (mm) 2
#>   Runoff (mm) 21 Deep drainage (mm) 151
 
#Switch to 'Sperry' transpiration mode
control <- defaultControl("Sperry")

#Initialize model input
x2 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
G2 <-growth(x2, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant cohort biomass (g/m2): 6177.36
#> Initial plant water content (mm): 6.71035
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  Year 2001:............
#> 
#> Final plant biomass (g/m2): 6323.12
#> Change in plant biomass (g/m2): 145.761
#> Plant biomass balance result (g/m2): 145.761
#> Plant biomass balance components:
#>   Structural balance (g/m2) 59 Labile balance (g/m2) 97
#>   Plant individual balance (g/m2) 155 Mortality loss (g/m2) 10
#> Final plant water content (mm): 6.72222
#> Final soil water content (mm): 274.718
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): 0.0118714
#> Plant water balance result (mm): 4.21353e-16
#> Change in soil water content (mm): -16.1567
#> Soil water balance result (mm): -16.1567
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 83 Net rainfall (mm) 379
#>   Infiltration (mm) 411 Infiltration excess (mm) 19 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 22  Herbaceous transpiration (mm) 0 Woody plant transpiration (mm) 243
#>   Plant extraction from soil (mm) 243  Plant water balance (mm) 0 Hydraulic redistribution (mm) 4
#>   Runoff (mm) 19 Deep drainage (mm) 163

#Switch to 'Sureau' transpiration mode
control <- defaultControl("Sureau")

#Initialize model input
x3 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
G3 <-growth(x3, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant cohort biomass (g/m2): 7776.21
#> Initial plant water content (mm): 6.71035
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  Year 2001:............
#> 
#> Final plant biomass (g/m2): 7641.04
#> Change in plant biomass (g/m2): -135.171
#> Plant biomass balance result (g/m2): -136.855
#> Plant biomass balance components:
#>   Structural balance (g/m2) 26 Labile balance (g/m2) -65
#>   Plant individual balance (g/m2) -40 Mortality loss (g/m2) 97
#> Final plant water content (mm): 6.6505
#> Final soil water content (mm): 279.247
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.059855
#> Plant water balance result (mm): -0.134463
#> Change in soil water content (mm): -11.6279
#> Soil water balance result (mm): -11.6279
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 0
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 82 Net rainfall (mm) 380
#>   Infiltration (mm) 410 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 31  Herbaceous transpiration (mm) 0 Woody plant transpiration (mm) 185
#>   Plant extraction from soil (mm) 185  Plant water balance (mm) -0 Hydraulic redistribution (mm) 0
#>   Runoff (mm) 21 Deep drainage (mm) 206
# }
      
```
