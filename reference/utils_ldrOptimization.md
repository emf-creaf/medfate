# Optimization of root distribution

Functions `utils_ldrExploration` and `utils_ldrOptimization` are used to
find optimum the species root distribution within `spwb`, given the
arguments `x`, `meteo` and `psi_crit`.

## Usage

``` r
utils_ldrExploration(
  x,
  meteo,
  cohorts = NULL,
  RZmin = 301,
  RZmax = 4000,
  V1min = 0.01,
  V1max = 0.94,
  resolution = 10,
  heat_stop = 0,
  transformation = "identity",
  verbose = FALSE,
  ...
)

utils_ldrOptimization(y, psi_crit, opt_mode = 1)
```

## Arguments

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- meteo:

  A data frame with daily meteorological data series (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- cohorts:

  A character string with the names of cohorts to be explored. If `NULL`
  then all cohorts are explored.

- RZmin:

  The minimum value of RZ (the rooting depth) to be explored (in mm)

- RZmax:

  The maximum value of RZ (the rooting depth) to be explored (in mm)

- V1min:

  The minimum value of V1 (the root proportion in the first soil layer)
  to be explored

- V1max:

  The maximum value of V1 (the root proportion in the first soil layer)
  to be explored

- resolution:

  An integer defining the number of values to obtain by discretization
  of the root parameters RZ and V1. The number of parameter combinations
  and therefore the computation cost increases increase with the square
  of resolution

- heat_stop:

  An integer defining the number of days during to discard from the
  calculation of the optimal root distribution. Usefull if the soil
  water content initialization is not certain

- transformation:

  Function to modify the size of Z intervals to be explored (by default,
  bins are equal).

- verbose:

  A logical value. Print the internal messages of the function?

- ...:

  Additional parameters to function
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md).

- y:

  The result of calling `utils_ldrExploration`.

- psi_crit:

  A numerical vector of length iqual to the number of species in the
  plot containing the species values of water potential inducing
  hydraulic failure (in MPa). Use `NA` values to skip optimization for
  particular plant cohorts.

- opt_mode:

  Optimization mode:

  - `opt_mode = 1` maximizes transpiration along the line of stress
    equal to `psi_crit` (Cabon et al. 2018). The optimization is based
    on the eco-hydrological equilibrium hypothesis (Eagleson, 1982),
    which is formulated here as the root distribution for which plant
    transpiration is maximized while the plant water potential is close
    to the species-defined critical value `psi_crit` (Cabon et
    al.,2018).

  - `opt_mode = 2` maximizes transpiration among combinations with
    stress according to `psi_crit`).

  - `opt_mode = 3` maximizes photosynthesis among combinations with
    stress according to `psi_crit`).

  - `opt_mode = 4` maximizes transpiration, subject to root construction
    constrains, among combinations with stress according to `psi_crit`).

  - `opt_mode = 5` maximizes photosynthesis, subject to root
    construction constrains, among combinations with stress according to
    `psi_crit`).

## Value

Function `utils_ldrExploration` returns a list containing a list
containing the explored RZ and V1 combinations as well as arrays with
the values of average daily plant transpiration, average daily net
photosynthesis and the minimum plant water potential for each cohort and
parameter combination.

Function `utils_ldrOptimization` returns a data frame with containing
the species index used in medfate, `psi_crit` and the optimized values
of V1 and the LDR parameters Z50 and Z95 (see
[`root_ldrDistribution`](https://emf-creaf.github.io/medfate/reference/root.md))
and as many rows as the number of species.

## Details

For each combination of the parameters RZ and V1 the function
`utils_ldrExploration` runs `spwb`, setting the total soil depth equal
to RZ. The root proportion in each soil layer is derived from V1, the
depth of the first soil layer and RZ using the LDR root distribution
model (Schenk and Jackson, 2002) and assuming that the depth containing
95 percent of the roots is equal to RZ. Function `utils_ldrOptimization`
takes the result of the exploration and tries to find optimum root
distribution parameters. `psi_crit`, the species specific water
potential inducing hydraulic failure, can be approached by the water
potential inducing 50 percent of loss of conductance for the and
gymnosperms and 88 percent for the angiosperms (Urli et al., 2013,
Brodribb et al., 2010). Details of the hypothesis and limitations of the
optimization method are given in Cabon et al. (2019).

## References

Brodribb, T.J., Bowman, D.J.M.S., Nichols, S., Delzon, S., Burlett, R.,
2010. Xylem function and growth rate interact to determine recovery
rates after exposure to extreme water deficit. New Phytol. 188, 533–542.
doi:10.1111/j.1469-8137.2010.03393.x

Cabon, A., Martínez-Vilalta, J., Poyatos, R., Martínez de Aragón, J., De
Cáceres, M. (2018) Applying the eco-hydrological equilibrium hypothesis
to estimate root ditribution in water-limited forests. Ecohydrology 11:
e2015.

Eagleson, P.S., 1982. Ecological optimality in water-limited natural
soil-vegetation systems: 1. Theory and hypothesis. Water Resour. Res.
18, 325–340. doi:10.1029/WR018i002p00325

Schenk, H.J., Jackson, R.B., 2002. The Global Biogeography of Roots.
Ecol. Monogr. 72, 311. doi:10.2307/3100092

Urli, M., Porte, A.J., Cochard, H., Guengant, Y., Burlett, R., Delzon,
S., 2013. Xylem embolism threshold for catastrophic hydraulic failure in
angiosperm trees. Tree Physiol. 33, 672–683. doi:10.1093/treephys/tpt030

## See also

[`utils_rockOptimization`](https://emf-creaf.github.io/medfate/reference/utils_rockOptimization.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`soil`](https://emf-creaf.github.io/medfate/reference/soil.md),
[`root_ldrDistribution`](https://emf-creaf.github.io/medfate/reference/root.md)

## Author

Antoine Cabon, WSL

Miquel De Cáceres Ainsa, CREAF
