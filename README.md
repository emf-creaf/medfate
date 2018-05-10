medfate
================

## Introduction

Package `medfate` is designed to assist forest scientists to simulate
forest functioning and dynamics, using cohort-based description of
forest stands. The models are parameterized for species of the
Mediterranean region (particularly for Catalonia, NE Spain), but forests
with different composition could be modelled with different parameter
sets.

## Package installation and documentation

Package `medfate` can be found at [CRAN](https://cran.r-project.org/),
but the version in this repository is very old. We recommend users to
download the latest stable versions GitHub as follows:

``` r
devtools::install_github("miquelcaceres/medfate")
```

Documentation on the models included in `medfate` and how to run them
using the package functions can be found at
(<http://vegmod.ctfc.cat/medfateweb>). Additionally, users can have help
to run package functions directly as package vignettes, by forcing their
inclusion in installation:

``` r
devtools::install_github("miquelcaceres/medfate", build_vignettes=TRUE)
```

## References

  - De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P,
    Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance
    model with forest inventory data to predict drought stress: the role
    of forest structural changes vs. climate changes. Agricultural and
    Forest Meteorology 213: 77-90
    (<doi:10.1016/j.agrformet.2015.06.012>).
