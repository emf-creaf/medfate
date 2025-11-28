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

A list of class 'growth' with the following elements:

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
#> Initial plant cohort biomass (g/m2): 5068.34
#> Initial plant water content (mm): 4.73001
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  Year 2001:............
#> 
#> Final plant biomass (g/m2): 5256.04
#> Change in plant biomass (g/m2): 187.703
#> Plant biomass balance result (g/m2): 187.703
#> Plant biomass balance components:
#>   Structural balance (g/m2) 104 Labile balance (g/m2) 92
#>   Plant individual balance (g/m2) 196 Mortality loss (g/m2) 8
#> Final plant water content (mm): 4.73772
#> Final soil water content (mm): 274.838
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): 0.00771124
#> Plant water balance result (mm): -0.00124258
#> Change in soil water content (mm): -16.0372
#> Soil water balance result (mm): -16.0372
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 0
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 370
#>   Infiltration (mm) 400 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 26  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 247
#>   Plant extraction from soil (mm) 247  Plant water balance (mm) -0 Hydraulic redistribution (mm) 3
#>   Runoff (mm) 21 Deep drainage (mm) 130
 
#Switch to 'Sperry' transpiration mode
control <- defaultControl("Sperry")

#Initialize model input
x2 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
G2 <-growth(x2, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant cohort biomass (g/m2): 6245.67
#> Initial plant water content (mm): 6.78662
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  Year 2001:............
#> 
#> Final plant biomass (g/m2): 6385.8
#> Change in plant biomass (g/m2): 140.134
#> Plant biomass balance result (g/m2): 140.134
#> Plant biomass balance components:
#>   Structural balance (g/m2) 54 Labile balance (g/m2) 96
#>   Plant individual balance (g/m2) 150 Mortality loss (g/m2) 10
#> Final plant water content (mm): 6.79935
#> Final soil water content (mm): 274.014
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): 0.012731
#> Plant water balance result (mm): 4.21944e-16
#> Change in soil water content (mm): -16.8613
#> Soil water balance result (mm): -16.8613
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): 7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 370
#>   Infiltration (mm) 402 Infiltration excess (mm) 19 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 23  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 240
#>   Plant extraction from soil (mm) 240  Plant water balance (mm) 0 Hydraulic redistribution (mm) 4
#>   Runoff (mm) 19 Deep drainage (mm) 143

#Switch to 'Sureau' transpiration mode
control <- defaultControl("Sureau")

#Initialize model input
x3 <- growthInput(exampleforest,examplesoil, SpParamsMED, control)

#Call simulation function
G3 <-growth(x3, examplemeteo, latitude = 41.82592, elevation = 100)
#> Initial plant cohort biomass (g/m2): 7826.12
#> Initial plant water content (mm): 6.78662
#> Initial soil water content (mm): 290.875
#> Initial snowpack content (mm): 0
#> Performing daily simulations
#> 
#>  Year 2001:............
#> 
#> Final plant biomass (g/m2): 7693.63
#> Change in plant biomass (g/m2): -132.488
#> Plant biomass balance result (g/m2): -133.783
#> Plant biomass balance components:
#>   Structural balance (g/m2) 28 Labile balance (g/m2) -66
#>   Plant individual balance (g/m2) -37 Mortality loss (g/m2) 96
#> Final plant water content (mm): 6.72728
#> Final soil water content (mm): 278.431
#> Final snowpack content (mm): 0
#> Change in plant water content (mm): -0.0593373
#> Plant water balance result (mm): -0.119326
#> Change in soil water content (mm): -12.4441
#> Soil water balance result (mm): -12.4441
#> Change in snowpack water content (mm): 0
#> Snowpack water balance result (mm): -7.10543e-15
#> Water balance components:
#>   Precipitation (mm) 513 Rain (mm) 462 Snow (mm) 51
#>   Interception (mm) 92 Net rainfall (mm) 371
#>   Infiltration (mm) 401 Infiltration excess (mm) 21 Saturation excess (mm) 0 Capillarity rise (mm) 0
#>   Soil evaporation (mm) 31  Herbaceous transpiration (mm) 14 Woody plant transpiration (mm) 184
#>   Plant extraction from soil (mm) 184  Plant water balance (mm) -0 Hydraulic redistribution (mm) 0
#>   Runoff (mm) 21 Deep drainage (mm) 185
# }
      
```
