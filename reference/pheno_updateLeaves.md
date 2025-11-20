# Leaf phenology

Function `pheno_leafDevelopmentStatus` returns the expanded status (0
to 1) of leaves according to the growth degree days required to start
bud burst and leaf unfolding, as dictated by a simple ecodormancy
(one-phase) model (Chuine et al. 2013). Function
`pheno_leafSenescenceStatus` returns the 0/1 senescence status of leaves
according to the one-phase senescence model of Delpierre et al. (2009)
on the basis of photoperiod and temperature. Function
`pheno_updateLeaves` updates the status of expanded leaves and dead
leaves of object `x` given the photoperiod, temperature and wind of a
given day. It applies the development model for 1 \< doy \< 180 and the
senescence model for 181 \> doy \> 365.

## Usage

``` r
pheno_leafDevelopmentStatus(Sgdd, gdd, unfoldingDD = 300)

pheno_leafSenescenceStatus(Ssen, sen)

pheno_updatePhenology(x, doy, photoperiod, tmean)

pheno_updateLeaves(x, wind, fromGrowthModel)
```

## Arguments

- Sgdd:

  Degree days required for leaf budburst (in Celsius).

- gdd:

  Cumulative degree days (in Celsius)

- unfoldingDD:

  Degree-days for complete leaf unfolding after budburst has occurred.

- Ssen:

  Threshold to start leaf senescence.

- sen:

  Cumulative senescence variable.

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- doy:

  Day of the year.

- photoperiod:

  Day length (in hours).

- tmean:

  Average day temperature (in Celsius).

- wind:

  Average day wind speed (in m/s).

- fromGrowthModel:

  Boolean flag to indicate that routine is called from
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)
  simulation function.

## Value

Function `pheno_leafDevelopmentStatus` returns a vector of values
between 0 and 1, whereas function `pheno_leafSenescenceStatus` returns a
vector of 0 (senescent) and 1 (expanded) values. The other two functions
do not return any value (see note).

## Note

Functions `pheno_updatePhenology` and `pheno_updateLeaves` modify the
input object `x`. The first modifies phenological state and the second
modifies the leaf area accordingly.

## References

Chuine, I., De Cortazar-Atauri, I.G., Kramer, K., Hänninen, H., 2013.
Plant development models. Phenology: An Integrative Environmental
Science. Springer, pp. 275–293.

Delpierre N, Dufrêne E, Soudani K et al (2009) Modelling interannual and
spatial variability of leaf senescence for three deciduous tree species
in France. Agric For Meteorol 149:938–948.
doi:10.1016/j.agrformet.2008.11.014

## See also

[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)

## Author

Miquel De Cáceres Ainsa, CREAF
