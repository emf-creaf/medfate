# Carbon-related functions

Low-level functions used in the calculation of carbon balance.

## Usage

``` r
carbon_belowgroundSapwoodStructuralBiomass(SA, L, V, woodDensity)

carbon_sapwoodStructuralBiomass(SA, H, L, V, woodDensity)

carbon_belowgroundHeartwoodStructuralBiomass(DBH, SA, L, V, woodDensity)

carbon_heartwoodStructuralBiomass(DBH, SA, H, L, V, woodDensity)

carbon_sapwoodStarchCapacity(SA, H, L, V, woodDensity, conduit2sapwood)

carbon_carbonCompartments(x, biomassUnits = "g_m2")

carbon_sugarStarchDynamicsLeaf(sugarConc, starchConc, eqSugarConc)

carbon_sugarStarchDynamicsStem(sugarConc, starchConc, eqSugarConc)

carbon_osmoticWaterPotential(sugarConc, temp, nonSugarConc)

carbon_sugarConcentration(osmoticWP, temp, nonSugarConc)

carbon_relativeSapViscosity(sugarConc, temp)

carbon_leafStructuralBiomass(LAI, N, SLA)

carbon_twigStructuralBiomass(LAI, N, SLA, r635)

carbon_leafStarchCapacity(LAI, N, SLA, leafDensity)

carbon_abovegroundSapwoodStructuralBiomass(SA, H, woodDensity)

carbon_abovegroundHeartwoodStructuralBiomass(DBH, SA, H, woodDensity)
```

## Arguments

- SA:

  Sapwood area (cm2).

- L:

  Coarse root length (mm) for each soil layer.

- V:

  Proportion of fine roots in each soil layer.

- woodDensity:

  Wood density (dry weight over volume).

- H:

  Plant height (cm).

- DBH:

  Tree diameter (cm).

- conduit2sapwood:

  Proportion of sapwood corresponding to conducive elements (vessels or
  tracheids) as opposed to parenchymatic tissue.

- x:

  An object of class
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- biomassUnits:

  A string for output biomass units, either "g_ind" (g per individual)
  or "g_m2" (g per square meter).

- sugarConc:

  Concentration of soluble sugars (mol/l).

- starchConc:

  Concentration of starch (mol/l)

- eqSugarConc:

  Equilibrium concentration of soluble sugars (mol/l).

- temp:

  Temperature (degrees Celsius).

- nonSugarConc:

  Concentration of inorganic solutes (mol/l).

- osmoticWP:

  Osmotic water potential (MPa).

- LAI:

  Leaf area index.

- N:

  Density (ind·ha-1).

- SLA:

  Specific leaf area (mm2/mg = m2/kg).

- leafDensity:

  Density of leaf tissue (dry weight over volume).

## Value

Values returned for each function are:

- `carbon_leafStarchCapacity`: Capacity of storing starch in the leaf
  compartment (mol gluc/ind.).

- `carbon_leafStructuralBiomass`: Leaf structural biomass (g dry/ind.)

- `carbon_sapwoodStarchCapacity`: Capacity of storing starch in the
  sapwood compartment (mol gluc/ind.).

- `carbon_sapwoodStructuralBiomass`: Sapwood structural biomass (g
  dry/ind.)

- `carbon_sapwoodStructuralLivingBiomass`: Living sapwood (parenchyma)
  structural biomass (g dry/ind.)

- `carbon_sugarConcentration`: Sugar concentration (mol gluc/l)

- `carbon_osmoticWaterPotential`: Osmotic component of water potential
  (MPa)

- `carbon_relativeSapViscosity`: Relative viscosity of sapwood with
  respect to pure water (according to Forst et al. (2002)).

- `carbon_sugarStarchDynamicsLeaf`: Rate of conversion from sugar to
  starch in leaf (mol gluc/l/s).

- `carbon_sugarStarchDynamicsStem`: Rate of conversion from sugar to
  starch in leaf (mol gluc/l/s).

- `carbon_carbonCompartments`: A data frame with the size of
  compartments for each plant cohort, in the specified units.

## References

Forst P, Wermer F, Delgado A (2002). On the pressure dependence of the
viscosity of aqueous sugar solutions. Rheol Acta 41: 369–374 DOI
10.1007/s00397-002-0238-y

## See also

[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
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

# Estimate carbon compartments
carbon_carbonCompartments(x1)
#>   LeafStorageVolume SapwoodStorageVolume LeafStarchMaximumConcentration
#> 1        39.5098502          20.22759439                       0.925123
#> 2         1.9930757           4.20844949                       0.925123
#> 3         0.1889106           0.06116612                       0.925123
#>   SapwoodStarchMaximumConcentration LeafStarchCapacity SapwoodStarchCapacity
#> 1                          4.625615         36.5514716            93.5650651
#> 2                          4.625615          1.8438402            19.4666673
#> 3                          4.625615          0.1747655             0.2829309
#>   LeafStructuralBiomass TwigStructuralBiomass TwigLivingStructuralBiomass
#> 1            235.258940            226.842893                   16.825556
#> 2             89.575518             72.186379                   26.973382
#> 3              6.819588              8.793534                    3.285819
#>   DeadLeafStructuralBiomass DeadTwigStructuralBiomass SapwoodStructuralBiomass
#> 1                         0                         0               4503.53340
#> 2                         0                         0                911.94903
#> 3                         0                         0                 17.02772
#>   SapwoodLivingStructuralBiomass HeartwoodStructuralBiomass
#> 1                      334.03935                   4624.499
#> 2                      340.76165                   3123.249
#> 3                        6.36263                      0.000
#>   AbovegroundWoodBiomass BelowgroundWoodBiomass FineRootBiomass LabileBiomass
#> 1             8930.18145             197.851192       29.436067    263.758249
#> 2             3776.25018             258.947653       48.165198    129.402075
#> 3               11.69174               5.335979        2.640632      4.705064
#>   StructuralBiomass DeadBiomass TotalLivingBiomass TotalBiomass
#> 1        4995.07130    4624.499          879.31817   9883.32879
#> 2        1121.87613    3123.249          634.87782   4374.52701
#> 3          35.28148       0.000           23.81373     39.98654
```
