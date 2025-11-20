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

carbon_leafStarchCapacity(LAI, N, SLA, leafDensity)

carbon_sapwoodStructuralBiomass(SA, H, L, V, woodDensity)

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

- L:

  Coarse root length (mm) for each soil layer.

- V:

  Proportion of fine roots in each soil layer.

- woodDensity:

  Wood density (dry weight over volume).

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
