# Carbon-related functions

Low-level functions used in the calculation of carbon balance.

## Usage

``` r
carbon_sugarStarchDynamicsLeaf(sugarConc, starchConc, eqSugarConc)

carbon_sugarStarchDynamicsStem(sugarConc, starchConc, eqSugarConc)

carbon_osmoticWaterPotential(sugarConc, temp, nonSugarConc)

carbon_sugarConcentration(osmoticWP, temp, nonSugarConc)

carbon_relativeSapViscosity(sugarConc, temp)

carbon_leafStructuralBiomass(LAI, N, SLA)

carbon_twigStructuralBiomass(LAI, N, SLA, r635)

carbon_leafStarchCapacity(LAI, N, SLA, leafDensity)

carbon_abovegroundSapwoodStructuralBiomass(SA, H, woodDensity)

carbon_belowgroundSapwoodStructuralBiomass(SA, L, V, woodDensity)

carbon_sapwoodStructuralBiomass(SA, H, L, V, woodDensity)

carbon_abovegroundHeartwoodStructuralBiomass(DBH, SA, H, woodDensity)

carbon_belowgroundHeartwoodStructuralBiomass(DBH, SA, L, V, woodDensity)

carbon_heartwoodStructuralBiomass(DBH, SA, H, L, V, woodDensity)

carbon_sapwoodStructuralLivingBiomass(
  SA,
  H,
  L,
  V,
  woodDensity,
  conduit2sapwood
)

carbon_sapwoodStarchCapacity(SA, H, L, V, woodDensity, conduit2sapwood)

carbon_carbonCompartments(x, biomassUnits = "g_m2")
```

## Arguments

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

- SA:

  Sapwood area (cm2).

- H:

  Plant height (cm).

- woodDensity:

  Wood density (dry weight over volume).

- L:

  Coarse root length (mm) for each soil layer.

- V:

  Proportion of fine roots in each soil layer.

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
#> 1        26.5664330          14.49487928                       0.925123
#> 2         4.0406733           5.17625488                       0.925123
#> 3         0.1679022           0.03056702                       0.925123
#>   SapwoodStarchMaximumConcentration LeafStarchCapacity SapwoodStarchCapacity
#> 1                          4.625615         24.5772184            67.0477318
#> 2                          4.625615          3.7381199            23.9433625
#> 3                          4.625615          0.1553302             0.1413913
#>   LeafStructuralBiomass TwigStructuralBiomass SapwoodStructuralBiomass
#> 1            165.109206            171.234885              3201.235698
#> 2            111.289246            138.097964              1146.794384
#> 3              6.149703              4.769236                 3.738437
#>   SapwoodLivingStructuralBiomass HeartwoodStructuralBiomass
#> 1                     244.444357                   6043.965
#> 2                     431.409712                   2937.490
#> 3                       1.406353                      0.000
#>   AbovegroundWoodBiomass BelowgroundWoodBiomass FineRootBiomass
#> 1            9044.809878             200.390824      23.2157679
#> 2            3822.186939             262.097660      20.9930167
#> 3               2.566923               1.171515       0.7933894
#>   StructuralBiomass LabileBiomass TotalLivingBiomass TotalBiomass
#> 1        9604.76056    199.895416          632.66475   9804.65598
#> 2        4354.66483    159.784125          723.47610   4514.44895
#> 3          15.45077      2.872282           11.22173     18.32305
```
