medfate - Mediterranean Forest Simulation
================

## Introduction

Package **medfate** is designed to assist forest scientists to simulate
forest functioning and dynamics, using cohort-based description of
forest stands. The package provides functions to simulate the following
processes:

-   Soil water balance (De Cáceres et al. 2015)
-   Plant hydraulics, transpiration and photosynthesis (De Cáceres et
    al. 2021)
-   Plant growth (in preparation)
-   Forest dynamics (in preparation)

The models are parameterized for species of the Mediterranean region
(particularly for Catalonia, NE Spain), but forests with different
composition could be modelled with different parameter sets.

## Package installation

Package **medfate** can be found at
[CRAN](https://CRAN.R-project.org/package=medfate), where it is updated
every few months. Users can also download and install the latest stable
versions GitHub as follows (required package `devtools` should be
installed/updated first):

``` r
devtools::install_github("emf-creaf/medfate")
```

Additionally, users can have help to run package functions directly as
package vignettes, by forcing their inclusion in installation:

``` r
devtools::install_github("emf-creaf/medfate", 
                         build_opts = c("--no-resave-data", "--no-manual"),
                         build_vignettes = TRUE)
```

## Documentation

-   The package includes a number of *vignettes* that illustrate how to
    run simulation models in **medfate**.

-   Complete documentation of the models included in the package can be
    found at the *reference book*
    <https://emf-creaf.github.io/medfatebook/index.html>.

## References

-   De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P,
    Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance
    model with forest inventory data to predict drought stress: the role
    of forest structural changes vs. climate changes. Agricultural and
    Forest Meteorology 213: 77-90
    (<https://doi.org/10.1016/j.agrformet.2015.06.012>).

-   De Cáceres M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L,
    Poyatos R, Cabon A, Granda V, Forner A, Valladares F,
    Martínez-Vilalta J (2021) Unravelling the effect of species mixing
    on water use and drought stress in holm oak forests: a modelling
    approach. Agricultural and Forest Meteorology 296
    (<https://doi.org/10.1016/j.agrformet.2020.108233>).
