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
#> 1        24.6407807          12.61517607                       0.925123
#> 2         2.2605172           4.77316160                       0.925123
#> 3         0.2065041           0.06686261                       0.925123
#>   SapwoodStarchMaximumConcentration LeafStarchCapacity SapwoodStarchCapacity
#> 1                          4.625615         22.7957532            58.3529483
#> 2                          4.625615          2.0912565            22.0788082
#> 3                          4.625615          0.1910417             0.3092807
#>   LeafStructuralBiomass TwigStructuralBiomass TwigLivingStructuralBiomass
#> 1            146.721992            141.473226                   10.493455
#> 2            101.595237             81.872731                   30.592814
#> 3              7.454707              9.612489                    3.591832
#>   DeadLeafStructuralBiomass DeadTwigStructuralBiomass SapwoodStructuralBiomass
#> 1                         0                         0               2808.68133
#> 2                         0                         0               1034.31920
#> 3                         0                         0                 18.61354
#>   SapwoodLivingStructuralBiomass HeartwoodStructuralBiomass
#> 1                     208.327554                   6319.351
#> 2                     386.486860                   3000.879
#> 3                       6.955191                      0.000
#>   AbovegroundWoodBiomass BelowgroundWoodBiomass FineRootBiomass LabileBiomass
#> 1             8930.18145             197.851192       24.832485    164.495920
#> 2             3776.25018             258.947653       54.628260    146.765933
#> 3               12.78061               5.832926        2.886558      5.143253
#>   StructuralBiomass DeadBiomass TotalLivingBiomass TotalBiomass
#> 1        3121.70904    6319.351          554.87141   9605.55626
#> 2        1272.41543    3000.879          720.06910   4420.06000
#> 3          38.56729       0.000           26.03154     43.71054
```
