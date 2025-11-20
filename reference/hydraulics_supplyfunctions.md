# Hydraulic supply functions

Set of functions used in the implementation of hydraulic supply
functions (Sperry and Love 2015).

## Usage

``` r
hydraulics_EXylem(
  psiPlant,
  psiUpstream,
  kxylemmax,
  c,
  d,
  allowNegativeFlux = TRUE,
  psiCav = 0
)

hydraulics_E2psiXylem(E, psiUpstream, kxylemmax, c, d, psiCav = 0)

hydraulics_E2psiXylemUp(E, psiDownstream, kxylemmax, c, d, psiCav = 0)

hydraulics_EVanGenuchten(psiRhizo, psiSoil, krhizomax, n, alpha, l = 0.5)

hydraulics_ECrit(psiUpstream, kxylemmax, c, d, pCrit = 0.001)

hydraulics_E2psiVanGenuchten(
  E,
  psiSoil,
  krhizomax,
  n,
  alpha,
  psiStep = -1e-04,
  psiMax = -10
)

hydraulics_E2psiTwoElements(
  E,
  psiSoil,
  krhizomax,
  kxylemmax,
  n,
  alpha,
  c,
  d,
  psiCav = 0,
  psiStep = -1e-04,
  psiMax = -10
)

hydraulics_E2psiBelowground(E, hydraulicNetwork, psiIni = as.numeric(c(0)))

hydraulics_E2psiAboveground(E, psiRootCrown, hydraulicNetwork)

hydraulics_E2psiNetwork(E, hydraulicNetwork, psiIni = as.numeric(c(0)))

hydraulics_supplyFunctionOneXylem(
  psiSoil,
  v,
  kstemmax,
  stemc,
  stemd,
  psiCav = 0,
  maxNsteps = 200L,
  dE = 0.01
)

hydraulics_supplyFunctionTwoElements(
  Emax,
  psiSoil,
  krhizomax,
  kxylemmax,
  n,
  alpha,
  c,
  d,
  psiCav = 0,
  dE = 0.1,
  psiMax = -10
)

hydraulics_supplyFunctionThreeElements(
  Emax,
  psiSoil,
  krhizomax,
  kxylemmax,
  kleafmax,
  n,
  alpha,
  stemc,
  stemd,
  leafc,
  leafd,
  psiCav = 0,
  dE = 0.1,
  psiMax = -10
)

hydraulics_supplyFunctionBelowground(
  hydraulicNetwork,
  minFlow = 0,
  pCrit = 0.001
)

hydraulics_supplyFunctionAboveground(
  Erootcrown,
  psiRootCrown,
  hydraulicNetwork
)

hydraulics_supplyFunctionNetwork(hydraulicNetwork, minFlow = 0, pCrit = 0.001)

hydraulics_regulatedPsiXylem(E, psiUpstream, kxylemmax, c, d, psiStep = -0.01)

hydraulics_regulatedPsiTwoElements(
  Emax,
  psiSoil,
  krhizomax,
  kxylemmax,
  n,
  alpha,
  c,
  d,
  dE = 0.1,
  psiMax = -10
)

hydraulics_initSperryNetworks(x)

hydraulics_supplyFunctionPlot(
  x,
  draw = TRUE,
  type = "E",
  speciesNames = FALSE,
  ylim = NULL
)
```

## Arguments

- psiPlant:

  Plant water potential (in MPa).

- psiUpstream:

  Water potential upstream (in MPa). In a one-component model
  corresponds to soil potential. In a two-component model corresponds to
  the potential inside the roots.

- kxylemmax:

  Maximum xylem hydraulic conductance (defined as flow per leaf surface
  unit and per pressure drop).

- c, d:

  Parameters of the Weibull function (generic xylem vulnerability
  curve).

- allowNegativeFlux:

  A boolean to indicate whether negative flux (i.e. from plant to soil)
  is allowed.

- psiCav:

  Minimum water potential (in MPa) experienced (for irreversible
  cavitation).

- E:

  Flow per surface unit.

- psiDownstream:

  Water potential upstream (in MPa).

- psiRhizo:

  Soil water potential (in MPa) in the rhizosphere (root surface).

- psiSoil:

  Soil water potential (in MPa). A scalar or a vector depending on the
  function.

- krhizomax:

  Maximum rhizosphere hydraulic conductance (defined as flow per leaf
  surface unit and per pressure drop).

- n, alpha, l:

  Parameters of the Van Genuchten function (rhizosphere vulnerability
  curve).

- pCrit:

  Critical water potential (in MPa).

- psiStep:

  Water potential precision (in MPa).

- psiMax:

  Minimum (maximum in absolute value) water potential to be considered
  (in MPa).

- hydraulicNetwork:

  List with the hydraulic characteristics of nodes in the hydraulic
  network.

- psiIni:

  Vector of initial water potential values (in MPa).

- psiRootCrown:

  Soil water potential (in MPa) at the root crown.

- v:

  Proportion of fine roots within each soil layer.

- kstemmax:

  Maximum stem xylem hydraulic conductance (defined as flow per leaf
  surface unit and per pressure drop).

- stemc, stemd:

  Parameters of the Weibull function for stems (stem xylem vulnerability
  curve).

- maxNsteps:

  Maximum number of steps in the construction of supply functions.

- dE:

  Increment of flow per surface unit.

- Emax:

  Maximum flow per surface unit.

- kleafmax:

  Maximum leaf hydraulic conductance (defined as flow per leaf surface
  unit and per pressure drop).

- leafc, leafd:

  Parameters of the Weibull function for leaves (leaf vulnerability
  curve).

- minFlow:

  Minimum flow in supply function.

- Erootcrown:

  Flow per surface unit at the root crown.

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- draw:

  A flag to indicate whether the supply function should be drawn or just
  returned.

- type:

  Plot type for `hydraulics_supplyFunctionPlot`, either `"E"`,
  `"ERhizo"`, `"StemPsi"`, `"RootPsi"` or `"dEdP"`).

- speciesNames:

  A flag to indicate the use of species names instead of cohort names in
  plots.

- ylim:

  Graphical parameter to override function defaults.

## Value

Values returned for each function are:

- `hydraulics_E2psiXylem`: The plant (leaf) water potential (in MPa)
  corresponding to the input flow, according to the xylem supply
  function and given an upstream (soil or root) water potential.

- `hydraulics_E2psiVanGenuchten`: The root water potential (in MPa)
  corresponding to the input flow, according to the rhizosphere supply
  function and given a soil water potential.

- `hydraulics_E2psiTwoElements`: The plant (leaf) water potential (in
  MPa) corresponding to the input flow, according to the rhizosphere and
  plant supply functions and given an input soil water potential.

- `hydraulics_E2psiNetwork`: The rhizosphere, root crown and plant (leaf
  water potential (in MPa) corresponding to the input flow, according to
  the vulnerability curves of rhizosphere, root and stem elements in a
  network.

- `hydraulics_Ecrit`: The critical flow according to the xylem supply
  function and given an input soil water potential.

- `hydraulics_EVanGenuchten`: The flow (integral of the vulnerability
  curve) according to the rhizosphere supply function and given an input
  drop in water potential (soil and rhizosphere).

- `hydraulics_EXylem`: The flow (integral of the vulnerability curve)
  according to the xylem supply function and given an input drop in
  water potential (rhizosphere and plant).

- `hydraulics_supplyFunctionOneXylem`,
  `hydraulics_supplyFunctionTwoElements` and
  `hydraulics_supplyFunctionNetwork`: A list with different numeric
  vectors with information of the two-element supply function:

  - `E`: Flow values (supply values).

  - `FittedE`: Fitted flow values (for
    `hydraulics_supplyFunctionTwoElements`).

  - `Elayers`: Flow values across the roots of each soil layer (only for
    `hydraulics_supplyFunctionNetwork`).

  - `PsiRhizo`: Water potential values at the root surface (only for
    `hydraulics_supplyFunctionNetwork`).

  - `PsiRoot`: Water potential values inside the root crown (not for
    `hydraulics_supplyFunctionOneXylem`).

  - `PsiPlant`: Water potential values at the canopy (leaf).

  - `dEdP`: Derivatives of the supply function.

- `hydraulics_supplyFunctionPlot`: If `draw = FALSE` a list with the
  result of calling `hydraulics_supplyFunctionNetwork` for each cohort.

- `hydraulics_regulatedPsiXylem`: Plant water potential after regulation
  (one-element loss function) given an input water potential.

- `hydraulics_regulatedPsiTwoElements`: Plant water potential after
  regulation (two-element loss function) given an input soil water
  potential.

## Details

Function `hydraulics_supplyFunctionPlot` draws a plot of the supply
function for the given `soil` object and network properties of each
plant cohort in `x`. Function `hydraulics_vulnerabilityCurvePlot` draws
a plot of the vulnerability curves for the given `soil` object and
network properties of each plant cohort in `x`.

## References

Sperry, J. S., F. R. Adler, G. S. Campbell, and J. P. Comstock. 1998.
Limitation of plant water use by rhizosphere and xylem conductance:
results from a model. Plant, Cell and Environment 21:347–359.

Sperry, J. S., and D. M. Love. 2015. What plant hydraulics can tell us
about responses to climate-change droughts. New Phytologist 207:14–27.

## See also

[`hydraulics_psi2K`](https://emf-creaf.github.io/medfate/reference/hydraulics_conductancefunctions.md),
[`hydraulics_maximumStemHydraulicConductance`](https://emf-creaf.github.io/medfate/reference/hydraulics_scalingconductance.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
kstemmax = 4 # in mmol·m-2·s-1·MPa-1
stemc = 3 
stemd = -4 # in MPa
psiVec = seq(-0.1, -7.0, by =-0.01)

#Vulnerability curve
kstem = unlist(lapply(psiVec, hydraulics_xylemConductance, kstemmax, stemc, stemd))
plot(-psiVec, kstem, type="l",ylab="Xylem conductance (mmol·m-2·s-1·MPa-1)", 
     xlab="Canopy pressure (-MPa)", lwd=1.5,ylim=c(0,kstemmax))

```
