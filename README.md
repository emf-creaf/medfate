medfate - Mediterranean forest functioning
================

## Introduction

Package **medfate** is designed to assist forest scientists to simulate
forest functioning and dynamics, using cohort-based description of
forest stands. The package provides functions to simulate the following
processes:

  - Soil water balance (De Cáceres et al. 2015)
  - Plant hydraulics, transpiration and photosynthesis
  - Plant growth

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
devtools::install_github("vegmod/medfate")
```

Additionally, users can have help to run package functions directly as
package vignettes, by forcing their inclusion in installation:

``` r
devtools::install_github("vegmod/medfate", 
                         build_opts = c("--no-resave-data", "--no-manual"))
```

## Documentation

The package installation includes a number of *vignettes* that
illustrate how to run simulation models in **medfate**. Vignettes and
documentation of package functions can also be read at the [GitHub pages
site](https://vegmod.github.io/medfate/) of the package.

The [medfate reference
book](https://vegmod.github.io/medfate/medfatebook/) provides full
documentation regarding the design and formulation of models included in
the package.

An interactive example application showing the behaviour of the water
balance model can be found at <http://vegmod.ctfc.cat/medfateweb>.

## References

  - De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P,
    Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance
    model with forest inventory data to predict drought stress: the role
    of forest structural changes vs. climate changes. Agricultural and
    Forest Meteorology 213: 77-90
    (<doi:10.1016/j.agrformet.2015.06.012>).
