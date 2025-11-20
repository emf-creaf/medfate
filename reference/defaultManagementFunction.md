# Default forest management actions

Function `defaultManagementFunction` implements actions for 'regular'
and 'irregular' management models of monospecific or mixed stands,
whereas function `defaultManagementArguments` returns a list with
default values for the parameters regulating management. Both functions
are meant to be used in simulations with
[`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md).

## Usage

``` r
defaultManagementFunction(x, args, verbose = FALSE)

defaultManagementArguments()
```

## Arguments

- x:

  An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)

- args:

  A list of arguments regulating management actions, e.g. the list
  returned by `defaultManagementArguments`

- verbose:

  A logical flag enabling console printing

## Value

Function `defaultManagementFunction` returns a list with the following
items:

- `"action"`: A string identifying the action performed (e.g.
  "thinning").

- `"N_tree_cut"`: A vector with the density of trees removed.

- `"Cover_shrub_cut"`: A vector with the cover of shrubs removed.

- `"planted_forest"`: An object of class
  [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
  with the new plant cohorts resulting from tree/shrub planting.

- `"management_args"`: A list of management arguments to be used in the
  next call to the management function.

Function `defaultManagementArguments` returns a list with default
arguments:

- `"type"`: Management model, either 'regular' or 'irregular'.

- `"targetTreeSpecies"`: Either `"all"` for unspecific cuttings or a
  numeric vector of target tree species to be selected for cutting
  operations.

- `"thinning"`: Kind of thinning to be applied in irregular models or in
  regular models before the final cuts. Options are 'below', 'above',
  'systematic', 'below-systematic', 'above-systematic' or a string with
  the proportion of cuts to be applied to different diameter sizes (see
  details).

- `"thinningMetric"`: The stand-level metric used to decide whether
  thinning is applied, either 'BA' (basal area), 'N' (density) or 'HB'
  (Hart-Becking index).

- `"thinningThreshold"`: The threshold value of the stand-level metric
  causing the thinning decision.

- `"thinningPerc"`: Percentage of stand's basal area to be removed in
  thinning operations.

- `"minThinningInterval"`: Minimum number of years between thinning
  operations.

- `"yearsSinceThinning"`: State variable to count the years since the
  last thinning ocurred.

- `"finalMeanDBH"`: Mean DBH threshold to start final cuts.

- `"finalPerc"`: String with percentages of basal area to be removed in
  final cuts, separated by '-' (e.g. "40-60-100").

- `"finalPreviousStage"`: Integer state variable to store the stage of
  final cuts ('0' before starting final cuts).

- `"finalYearsBetweenCuts"`: Number of years separating final cuts.

- `"finalYearsToCut"`: State variable to count the years to be passed
  before new final cut is applied.

- `"plantingSpecies"`: Species code to be planted. If missing, planting
  does not occur and only natural regeneration is allowed.

- `"plantingDBH"`: Initial DBH (cm) of planted species.

- `"plantingHeight"`: Initial height (cm) of planted species.

- `"plantingDensity"`: Initial density (ind./ha) of the planted species.

- `"understoryMaximumCover"`: Percentage of overall shrub cover to be
  left after any silvicultural intervention. If missing, shrub cover
  will not be left unmodified.

## Details

This function implements silvicultural actions following either
'regular' or 'irregular' management models. Irregular models are
implemented by executing thinning operations only, whereas regular
models include both thinning and a set of final cuts. Thinning occurs
anytime a stand-level metric (e.g. basal area) crosses a given
threshold, and different kinds of thinning operations are allowed.
Unrealistic high frequency thinning can be avoided by setting a minimum
number of years to happen between thinning operations. Final cuts start
whenever mean DBH exceeds a given threshold, and may include different
cuts separated a number of years. The function can be applied to target
management of specific taxa (instead of assuming a monospecific stand),
but the thresholds that determine thinning operations apply to
stand-level metrics. Mean DBH will be calculated for the target species
only. Planting is only allowed under regular management models, and is
applied after the last final cut. Understory clearings are assumed to
occur anytime there is an intervention on trees, an only a residual
shrub cover is left.

*Thinning types*:

- `above`: Extract largest trees (according to DBH) until thinning
  objective is met.

- `below`: Extract smallest trees (according to DBH) until thinning
  objective is met.

- `systematic`: Extract equally from all size classes until thinning
  objective is met.

- `above-systematic`: Extract half the objective as systematic thinning
  and the other hald as above thinning.

- `below-systematic`: Extract half the objective as systematic thinning
  and the other hald as below thinning.

- `free string`: A string specifying the proportion of tree cuts from
  size classes, with size classes separated by "/" and each one composed
  of a number specifying the upper limit and a number indicating its
  proportion, separated by "-" (e.g. "10-50/40-30/60-20").

## See also

[`fordyn`](https://emf-creaf.github.io/medfate/reference/fordyn.md)

## Author

Miquel De Cáceres Ainsa, CREAF

Aitor Améztegui, UdL

Jose-Ramon Gonzalez Olabarria, CTFC

## Examples

``` r
# Load example forest object
data(exampleforest)
  
# Define arguments
args = defaultManagementArguments()
  
# Call management function
f = defaultManagementFunction(exampleforest, args)
  
#list names
names(f)
#> [1] "action"          "N_tree_cut"      "Cover_shrub_cut" "planted_forest" 
#> [5] "management_args"
  
# Action performed
f$action
#> [1] "thinning"
  
# Number of trees cut for each cohort
f$N_tree_cut
#> [1]   9.76362 384.00000
  
# Percent cover of shrubs removed
f$Cover_shrub_cut
#> [1] 0.75
  
```
