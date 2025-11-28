# Single-day forest growth

Function `growth_day` performs water and carbon balance for a single
day.

## Usage

``` r
growth_day(
  x,
  date,
  meteovec,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  runon = 0,
  lateralFlows = NULL,
  waterTableDepth = NA_real_,
  modifyInput = TRUE
)
```

## Arguments

- x:

  An object of class
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- date:

  Date as string "yyyy-mm-dd".

- meteovec:

  A named numerical vector with weather data. See variable names in
  parameter `meteo` of
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md).

- latitude:

  Latitude (in degrees).

- elevation, slope, aspect:

  Elevation above sea level (in m), slope (in degrees) and aspect (in
  degrees from North).

- runon:

  Surface water amount running on the target area from upslope (in mm).

- lateralFlows:

  Lateral source/sink terms for each soil layer (interflow/to from
  adjacent locations) as mm/day.

- waterTableDepth:

  Water table depth (in mm). When not missing, capillarity rise will be
  allowed if lower than total soil depth.

- modifyInput:

  Boolean flag to indicate that the input `x` object is allowed to be
  modified during the simulation.

## Value

Function `growth_day()` returns a list of class `growth_day` with the
same elements as
[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
and the following:

- `"LabileCarbonBalance"`: A data frame with labile carbon balance
  results for plant cohorts, with elements:

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

- `"PlantBiomassBalance"`: A data frame with plant biomass balance
  results for plant cohorts, with elements:

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

- `"PlantStructure"`: A data frame with area and biomass values for
  compartments of plant cohorts, with elements:

  - `"LeafBiomass"`: Leaf structural biomass (in g dry) for an average
    individual of each plant cohort.

  - `"SapwoodBiomass"`: Sapwood structural biomass (in g dry) for an
    average individual of each plant cohort.

  - `"FineRootBiomass"`: Fine root biomass (in g dry) for an average
    individual of each plant cohort.

  - `"LeafArea"`: Leaf area (in m2) for an average individual of each
    plant cohort.

  - `"SapwoodArea"`: Sapwood area (in cm2) for an average individual of
    each plant cohort.

  - `"FineRootArea"`: Fine root area (in m2) for an average individual
    of each plant cohort.

  - `"HuberValue"`: Sapwood area to (target) leaf area (in cm2/m2).

  - `"RootAreaLeafArea"`: The ratio of fine root area to (target) leaf
    area (in m2/m2).

  - `"DBH"`: Diameter at breast height (in cm) for an average individual
    of each plant cohort.

  - `"Height"`: Height (in cm) for an average individual of each plant
    cohort.

- `"GrowthMortality"`: A data frame with growth and mortality rates for
  plant cohorts, with elements:

  - `"LAgrowth"`: Leaf area growth (in m2·day-1) for an average
    individual of each plant cohort.

  - `"SAgrowth"`: Sapwood area growth rate (in cm2·day-1) for an average
    individual of each plant cohort.

  - `"FRAgrowth"`: Fine root area growth (in m2·day-1) for an average
    individual of each plant cohort.

  - `"StarvationRate"`: Mortality rate from starvation (ind/d-1).

  - `"DessicationRate"`: Mortality rate from dessication (ind/d-1).

  - `"MortalityRate"`: Mortality rate (any cause) (ind/d-1).

## Details

Detailed model description is available in the medfate book.

Forest growth simulations allow using different sub-models for bulk soil
water flows and different sub-models of transpiration and photosynthesis
(see details in
[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)).

## References

De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
forest inventory data to predict drought stress: the role of forest
structural changes vs. climate changes. Agricultural and Forest
Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).

De Cáceres M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L,
Poyatos R, Cabon A, Granda V, Forner A, Valladares F, Martínez-Vilalta J
(2021) Unravelling the effect of species mixing on water use and drought
stress in holm oak forests: a modelling approach. Agricultural and
Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).

Granier A, Bréda N, Biron P, Villette S (1999) A lumped water balance
model to evaluate duration and intensity of drought constraints in
forest stands. Ecol Modell 116:269–283.
https://doi.org/10.1016/S0304-3800(98)00205-1.

Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022)
SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations
of plant water status and drought-induced mortality at the ecosystem
level. Geoscientific Model Development 15, 5593-5626
(doi:10.5194/gmd-15-5593-2022).

Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S.
Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to
the environment from the optimization of photosynthetic gain and
hydraulic cost. Plant Cell and Environment 40, 816-830 (doi:
10.1111/pce.12852).

## See also

[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
[`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md),
[`plot.growth_day`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md)

## Author

- Miquel De Cáceres Ainsa, CREAF

- Nicolas Martin-StPaul, URFM-INRAE

## Examples

``` r
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Define soil parameters
examplesoil <- defaultSoilParams(4)

# Day to be simulated
d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

#Simulate water and carbon balance for one day only (Granier mode)
control <- defaultControl("Granier")
x4  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
sd4 <- growth_day(x4, date, meteovec,
                latitude = 41.82592, elevation = 100, slope=0, aspect=0)

#Simulate water and carbon balance for one day only (Sperry mode)
control <- defaultControl("Sperry")
x5  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
sd5 <- growth_day(x5, date, meteovec,
                latitude = 41.82592, elevation = 100, slope=0, aspect=0)

#Simulate water and carbon balance for one day only (Sureau mode)
control <- defaultControl("Sureau")
x6  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
sd6 <- growth_day(x6, date, meteovec,
                latitude = 41.82592, elevation = 100, slope=0, aspect=0)
```
