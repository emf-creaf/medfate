# Root functions

Functions to calculate properties of fine/coarse roots within the soil,
given root system parameters and soil layer definition.

## Usage

``` r
root_conicDistribution(Zcone, d)

root_ldrDistribution(Z50, Z95, Z100, d)

root_individualRootedGroundArea(VolInd, V, d, rfc)

root_specificRootSurfaceArea(specificRootLength, rootTissueDensity)

root_fineRootRadius(specificRootLength, rootTissueDensity)

root_fineRootHalfDistance(rootLengthDensity)

root_fineRootAreaIndex(
  Ksoil,
  krhizo,
  lai,
  specificRootLength,
  rootTissueDensity,
  rootLengthDensity
)

root_fineRootBiomass(
  Ksoil,
  krhizo,
  lai,
  N,
  specificRootLength,
  rootTissueDensity,
  rootLengthDensity
)

root_rhizosphereMaximumConductance(
  Ksoil,
  fineRootBiomass,
  lai,
  N,
  specificRootLength,
  rootTissueDensity,
  rootLengthDensity
)

root_fineRootSoilVolume(fineRootBiomass, specificRootLength, rootLengthDensity)

root_coarseRootSoilVolumeFromConductance(
  Kmax_rootxylem,
  VCroot_kmax,
  Al2As,
  v,
  d,
  rfc
)

root_coarseRootLengthsFromVolume(VolInd, v, d, rfc)

root_coarseRootLengths(v, d, depthWidthRatio = 1)

root_coarseRootSoilVolume(v, d, depthWidthRatio = 1)

root_horizontalProportions(poolProportions, VolInd, N, V, d, rfc)
```

## Arguments

- Zcone:

  A vector of depths (in mm) corresponding to the root cone tip.

- d:

  The width (in mm) corresponding to each soil layer.

- Z50:

  A vector of depths (in mm) corresponding to 50% of roots.

- Z95:

  A vector of depths (in mm) corresponding to 95% of roots.

- Z100:

  A vector of depths (in mm) corresponding to 100% of roots.

- VolInd:

  Volume of soil (in m3) occupied by coarse roots per individual.

- V:

  Matrix of proportions of fine roots (cohorts x soil layers).

- rfc:

  Percentage of rock fragment content (volume basis) for each layer.

- specificRootLength:

  Specific fine root length (length of fine roots over weight).

- rootTissueDensity:

  Fine root tissue density (weight over volume at turgidity).

- rootLengthDensity:

  Fine root length density (length of fine roots over soil volume;
  cm/cm3)

- Ksoil:

  Soil saturated conductivity (mmol·m-1·s-1·MPa-1).

- krhizo:

  Rhizosphere maximum conductance per leaf area (mmol·m-2·s-1·MPa-1).

- lai:

  Leaf area index.

- N:

  Density of individuals per hectare.

- fineRootBiomass:

  Biomass of fine roots (g).

- Kmax_rootxylem:

  Sapwood-specific hydraulic conductivity of root xylem (in kg
  H2O·s-1·m-1·MPa-1).

- VCroot_kmax:

  Root xylem maximum conductance per leaf area (mmol·m-2·s-1·MPa-1).

- Al2As:

  Leaf area to sapwood area ratio (in m2·m-2).

- v:

  Vector of proportions of fine roots in each soil layer.

- depthWidthRatio:

  Ratio between radius of the soil layer with the largest radius and
  maximum rooting depth.

- poolProportions:

  Division of the stand area among plant cohorts (proportions).

## Value

See details.

## Details

- `root_conicDistribution()` assumes a (vertical) conic distribution of
  fine roots, whereas `root_ldrDistribution()` distributes fine roots
  according to the linear dose response model of Schenck & Jackson
  (2002). Return a matrix of fine root proportions in each layer with as
  many rows as elements in `Z` (or `Z50`) and as many columns as soil
  layers.

- `root_coarseRootLengths()` and `root_coarseRootLengthsFromVolume()`
  estimate the length of coarse roots (mm) for each soil layer,
  including axial and radial lengths.

- `root_coarseRootSoilVolume` estimates the soil volume (m3) occupied by
  coarse roots of an individual.

- `root_coarseRootSoilVolumeFromConductance` estimates the soil volume
  (m3) occupied by coarse roots of an individual from root xylem
  conductance.

- `root_fineRootHalfDistance()` calculates the half distance (cm)
  between neighbouring fine roots.

- `root_fineRootRadius()` calculates the radius of fine roots (cm).

- `root_fineRootAreaIndex()` estimates the fine root area index for a
  given soil conductivity and maximum rhizosphere conductance.

- `root_fineRootBiomass()` estimates the biomass of fine roots (g
  dry/individual) for a given soil conductivity and maximum rhizosphere
  conductance.

- `root_rhizosphereMaximumConductance()` is the inverse of the
  preceeding function, i.e. it estimates rhizosphere conductance from
  soil conductivity and fine root biomass.

- `root_fineRootSoilVolume()` calculates the soil volume (m3) occupied
  with fine roots.

- `root_specificRootSurfaceArea()` returns the specific fine root area
  (cm2/g).

- `root_individualRootedGroundArea()` calculates the area (m2) covered
  by roots of an individual, for each soil layer.

- `root_horizontalProportions()` calculates the (horizontal) proportion
  of roots of each cohort in the water pool corresponding to itself and
  that of other cohorts, for each soil layer. Returns a list (with as
  many elements as cohorts) with each element being a matrix.

## References

Schenk, H., Jackson, R., 2002. The global biogeography of roots. Ecol.
Monogr. 72, 311–328.

Sperry, J. S., Y. Wang, B. T. Wolfe, D. S. Mackay, W. R. L. Anderegg, N.
G. Mcdowell, and W. T. Pockman. 2016. Pragmatic hydraulic theory
predicts stomatal responses to climatic water deficits. New Phytologist
212, 577–589.

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

ntree <- nrow(exampleforest$treeData)

#Initialize soil with default soil params
s <- defaultSoilParams(4)

#Calculate conic root system for trees
V1 <- root_conicDistribution(Z=rep(2000,ntree), s$widths)            
print(V1)
#>          [,1]     [,2]  [,3] [,4]
#> [1,] 0.385875 0.489125 0.125    0
#> [2,] 0.385875 0.489125 0.125    0
     
#Calculate LDR root system for trees (Schenck & Jackson 2002)
V2 <- root_ldrDistribution(Z50 = rep(200,ntree), 
                           Z95 = rep(1000,ntree),
                           Z100 = rep(NA, ntree), s$widths)
print(V2)     
#>           [,1]      [,2]       [,3]       [,4]
#> [1,] 0.6799879 0.2737911 0.03567632 0.01054468
#> [2,] 0.6799879 0.2737911 0.03567632 0.01054468
```
