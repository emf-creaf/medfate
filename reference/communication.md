# Internal communication

Functions for internal communication. Not to be called by users.

## Usage

``` r
aspwb_day_inner(
  internalCommunication,
  x,
  date,
  meteovec,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  runon = 0,
  lateralFlows = NULL,
  waterTableDepth = NA_real_,
  modifyInput = TRUE
)

copy_model_output(internalCommunication, x, model)

general_communication_structures(
  numCohorts,
  nlayers,
  ncanlayers,
  ntimesteps,
  model
)

instance_communication_structures(x, model)

growth_day_inner(
  internalCommunication,
  x,
  date,
  meteovec,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  runon = 0,
  lateralFlows = NULL,
  waterTableDepth = NA_real_,
  modifyInput = TRUE
)

spwb_day_inner(
  internalCommunication,
  x,
  date,
  meteovec,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  runon = 0,
  lateralFlows = NULL,
  waterTableDepth = NA_real_,
  modifyInput = TRUE
)
```

## Arguments

- internalCommunication:

  List for internal communication.

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  or
  [`growthInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- date:

  Date as string "yyyy-mm-dd".

- meteovec:

  A named numerical vector with weather data. See variable names in
  parameter `meteo` of
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md).

- latitude:

  Latitude (in degrees).

- elevation, slope, aspect:

  Elevation above sea level (in m), slope (in degrees) and aspect (in
  degrees from North).

- runon:

  Surface water amount running on the target area from upslope (in mm).

- lateralFlows:

  Lateral source/sink terms for each soil layer (interflow/to from
  adjacent locations) as mm/day.

- waterTableDepth:

  Water table depth (in mm). When not missing, capillarity rise will be
  allowed if lower than total soil depth.

- modifyInput:

  Boolean flag to indicate that the input `x` object is allowed to be
  modified during the simulation.

- model:

  String for model, either "spwb" or "growth".
