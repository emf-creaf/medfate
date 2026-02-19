# Rcpp module: Running simulations

This Rcpp module allows running simulations using C++ background
objects. They are intended for non-standard users.

## Details

The module contains the following items:

classes:

- single_runner:

  An S4 class to run spwb_day() or growth_day() simulations on a single
  site.

- multiple_runner:

  An S4 class to run spwb_day() or growth_day() simulations on multiple
  sites, if required using parallelization.

- watershed_runner:

  An S4 class to run one day of watershed simulations.

## Examples

``` r
show( runners )
#> Uninitialized module named "runners" from package "medfate"
```
