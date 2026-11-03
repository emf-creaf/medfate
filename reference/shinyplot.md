# Shiny app with interactive plots

Creates a shiny app with interactive plots for simulation results and
evaluation

## Usage

``` r
shinyplot(x, ...)

# S3 method for class 'growth'
shinyplot(x, measuredData = NULL, ...)

# S3 method for class 'aspwb'
shinyplot(x, measuredData = NULL, ...)

# S3 method for class 'spwb'
shinyplot(x, measuredData = NULL, ...)

# S3 method for class 'pwb'
shinyplot(x, measuredData = NULL, ...)

# S3 method for class 'fordyn'
shinyplot(x, measuredData = NULL, ...)

# S3 method for class 'growth_day'
shinyplot(x, ...)

# S3 method for class 'spwb_day'
shinyplot(x, ...)

# S3 method for class 'pwb_day'
shinyplot(x, ...)
```

## Arguments

- x:

  An object of the right class.

- ...:

  Additional parameters.

- measuredData:

  A data frame with observed/measured values (see
  [`evaluation_plot`](https://emf-creaf.github.io/medfate/reference/evaluation.md)).

## Value

An object that represents the shiny app

## Details

Only run this function in interactive mode. When `measuredData` is not
`NULL`, an additional panel is shown for evaluation plots.

## See also

[`plot.spwb`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md),
[`evaluation_plot`](https://emf-creaf.github.io/medfate/reference/evaluation.md)

## Author

Miquel De Cáceres Ainsa, CREAF
