# Multiple model runs and function factories for optimization routines

Function factories to generate functions to be used in model
calibration, uncertainty or sensitivity analysis.

## Usage

``` r
multiple_runs(
  parMatrix,
  x,
  meteo,
  latitude,
  elevation = NA,
  slope = NA,
  aspect = NA,
  summary_function = NULL,
  args = NULL,
  verbose = TRUE
)

optimization_function(
  parNames,
  x,
  meteo,
  latitude,
  elevation = NA,
  slope = NA,
  aspect = NA,
  summary_function,
  args = NULL
)

optimization_evaluation_function(
  parNames,
  x,
  meteo,
  latitude,
  elevation = NA,
  slope = NA,
  aspect = NA,
  measuredData,
  type = "SWC",
  cohorts = NULL,
  temporalResolution = "day",
  metric = "loglikelihood"
)

optimization_multicohort_function(
  cohortParNames,
  cohortNames,
  x,
  meteo,
  latitude,
  otherParNames = NULL,
  elevation = NA,
  slope = NA,
  aspect = NA,
  summary_function,
  args = NULL
)

optimization_evaluation_multicohort_function(
  cohortParNames,
  cohortNames,
  x,
  meteo,
  latitude,
  otherParNames = NULL,
  elevation = NA,
  slope = NA,
  aspect = NA,
  measuredData,
  type = "SWC",
  cohorts = cohortNames,
  temporalResolution = "day",
  metric = "loglikelihood"
)
```

## Arguments

- parMatrix:

  A matrix of parameter values with runs in rows and parameters in
  columns. Column names should follow parameter modification naming
  rules (see examples and naming rules in
  [`modifyInputParams`](https://emf-creaf.github.io/medfate/reference/modifyParams.md)).

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  or
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- meteo, latitude, elevation, slope, aspect:

  Additional parameters to simulation functions
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md).

- summary_function:

  A function whose input is the result of
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) or
  [`growth`](https://emf-creaf.github.io/medfate/reference/growth.md).
  The function must return a numeric scalar in the case of
  `optimization_function`, but is not restricted in the case of
  `multiple_runs`.

- args:

  A list of additional arguments of `optimization_function`.

- verbose:

  A flag to indicate extra console output.

- parNames:

  A string vector of parameter names (see examples and naming rules in
  [`modifyInputParams`](https://emf-creaf.github.io/medfate/reference/modifyParams.md)).

- measuredData:

  A data frame with observed/measured values. Dates should be in row
  names, whereas columns should be named according to the type of output
  to be evaluated (see details).

- type:

  A string with the kind of model output to be evaluated. Accepted
  values are `"SWC"` (soil moisture content), `"REW"` relative
  extractable water, `"ETR"` (total evapotranspiration), `"E"`
  (transpiration per leaf area), `"LFMC"` (live fuel moisture content)
  and `"WP"` (plant water potentials).

- cohorts:

  A string or a vector of strings with the cohorts to be compared (e.g.
  "T1_68"). If several cohort names are provided, the function
  `optimization_cohorts_function` evaluates the performance for each one
  and provides the mean value. If `NULL` results for the first cohort
  will be evaluated.

- temporalResolution:

  A string to indicate the temporal resolution of the model evaluation,
  which can be "day", "week", "month" or "year". Observed and modelled
  values are aggregated temporally (using either means or sums) before
  comparison.

- metric:

  An evaluation metric (see
  [`evaluation_metric`](https://emf-creaf.github.io/medfate/reference/evaluation.md)).

- cohortParNames:

  A string vector of vegetation parameter names for cohorts (e.g. 'Z95'
  or 'psiExtract').

- cohortNames:

  A string vector of cohort names. All cohorts will be given the same
  parameter values for each parameter in 'cohortParNames'.

- otherParNames:

  A string vector of parameter names (see examples and naming rules in
  [`modifyInputParams`](https://emf-creaf.github.io/medfate/reference/modifyParams.md))
  for non-vegetation parameters (i.e. control parameters and soil
  parameters).

## Value

Function `multiple_runs` returns a list, whose elements are either the
result of calling simulation models or the result of calling
`summary_function` afterwards.

Function `optimization_function` returns a function whose parameters are
parameter values and whose return is a prediction scalar (e.g. total
transpiration).

Function `optimization_evaluation_function` returns a function whose
parameters are parameter values and whose return is an evaluation metric
(e.g. loglikelihood of the data observations given model predictions).
If evaluation data contains information for different cohorts (e.g.
plant water potentials or transpiration rates) then the evaluation is
performed for each cohort and the metrics are averaged.

Function `optimization_multicohorts_function` returns a function whose
parameters are parameter values and whose return is a prediction scalar
(e.g. total transpiration). The difference with `optimization_function`
is that multiple cohorts are set to the same parameter values.

Function `optimization_evaluation_multicohort_function` returns a
function whose parameters are parameter values and whose return is an
evaluation metric (e.g. loglikelihood of the data observations given
model predictions). If evaluation data contains information for
different cohorts (e.g. plant water potentials or transpiration rates)
then the evaluation is performed for each cohort and the metrics are
averaged. The difference with `optimization_evaluation_function` is that
multiple cohorts are set to the same parameter values.

## Details

See
[`evaluation`](https://emf-creaf.github.io/medfate/reference/evaluation.md)
for details regarding how to specify measured data.

Functions produced by these function factories should be useful for
sensitivity analyses using package 'sensitivity'.

Parameter naming (i.e. `parNames`) should follow the rules specified in
section details of
[`modifyInputParams`](https://emf-creaf.github.io/medfate/reference/modifyParams.md).
The exception to the naming rules applies when multiple cohorts are to
be modified to the same values with functions
`optimization_multicohort_function` and
`optimization_evaluation_multicohort_function`. Then, only a vector of
parameter names is supplied for `cohortParNames`.

## See also

[`evaluation_metric`](https://emf-creaf.github.io/medfate/reference/evaluation.md),
[`modifyInputParams`](https://emf-creaf.github.io/medfate/reference/modifyParams.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`growth`](https://emf-creaf.github.io/medfate/reference/growth.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
# \donttest{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

#Initialize control parameters
control <- defaultControl("Granier")

#Initialize input
x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

# Cohort name for Pinus halepensis
PH_coh <- paste0("T1_", SpParamsMED$SpIndex[SpParamsMED$Name=="Pinus halepensis"])
PH_coh 
#> [1] "T1_148"

#Parameter names of interest
parNames <- c(paste0(PH_coh,"/Z50"), paste0(PH_coh,"/Z95"))

#Specify parameter matrix
parMatrix <- cbind(c(200,300), c(500,1000))
colnames(parMatrix) <- parNames

#Define a summary function as the total transpiration over the simulated period
sf<-function(x) {sum(x$WaterBalance$Transpiration, na.rm=TRUE)}

#Perform two runs and evaluate the summary function
multiple_runs(parMatrix, 
              x1, examplemeteo, latitude = 42, elevation = 100,
              summary_function = sf)
#> 1. Parameter values = [200, 500] f = 247.055219438414
#> 2. Parameter values = [300, 1000] f = 248.279875459894
#> [[1]]
#> [1] 247.0552
#> 
#> [[2]]
#> [1] 248.2799
#> 

#Load observed data (in this case the same simulation results with some added error)  
# Generate a prediction function for total transpiration over the simulated period
# as a function of parameters "Z50" and "Z95" for Pinus halepensis cohort 
of<-optimization_function(parNames = parNames,
                          x = x1,
                          meteo = examplemeteo, 
                          latitude = 41.82592, elevation = 100,
                          summary_function = sf)

# Evaluate for the values of the parameter matrix
of(parMatrix[1, ])
#> [1] 247.1513
of(parMatrix)
#> [1] 247.1513 248.3765


# Generate a loglikelihood function for soil water content
# as a function of parameters "Z50" and "Z95" for Pinus halepensis cohort 
data(exampleobs)
oef<-optimization_evaluation_function(parNames = parNames,
                                      x = x1,
                                      meteo = examplemeteo, latitude = 41.82592, elevation = 100,
                                      measuredData = exampleobs, type = "SWC", 
                                      metric = "loglikelihood")

# Loglikelihood for the values of the parameter matrix
oef(parMatrix[1, ])
#> [1] 996.0599
oef(parMatrix)
#> [1] 996.0599 988.3677
# }
```
