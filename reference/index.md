# Package index

## Forest definition

Forest definition and manipulation functions

- [`forest`](https://emf-creaf.github.io/medfate/reference/forest.md)
  [`exampleforest`](https://emf-creaf.github.io/medfate/reference/forest.md)
  [`exampleforest2`](https://emf-creaf.github.io/medfate/reference/forest.md)
  : Description of a forest stand.
- [`emptyforest()`](https://emf-creaf.github.io/medfate/reference/emptyforest.md)
  : Creation of an empty forest
- [`tree2forest()`](https://emf-creaf.github.io/medfate/reference/tree2forest.md)
  [`shrub2forest()`](https://emf-creaf.github.io/medfate/reference/tree2forest.md)
  : Single-cohort forests
- [`forest_mapTreeTable()`](https://emf-creaf.github.io/medfate/reference/forest_mapWoodyTables.md)
  [`forest_mapShrubTable()`](https://emf-creaf.github.io/medfate/reference/forest_mapWoodyTables.md)
  [`forest_mapWoodyTables()`](https://emf-creaf.github.io/medfate/reference/forest_mapWoodyTables.md)
  : Map forest plot data
- [`forest_mergeTrees()`](https://emf-creaf.github.io/medfate/reference/forest_simplification.md)
  [`forest_mergeShrubs()`](https://emf-creaf.github.io/medfate/reference/forest_simplification.md)
  [`forest_reduceToDominant()`](https://emf-creaf.github.io/medfate/reference/forest_simplification.md)
  : Forest complexity reduction
- [`poblet_trees`](https://emf-creaf.github.io/medfate/reference/poblet_trees.md)
  : Example forest inventory data

## Species parameters

Species parameter data and utility functions.

- [`SpParams`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  [`SpParamsDefinition`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  [`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)
  : Data tables with species parameter definitions and values
- [`modifySpParams()`](https://emf-creaf.github.io/medfate/reference/modifyParams.md)
  [`modifyCohortParams()`](https://emf-creaf.github.io/medfate/reference/modifyParams.md)
  [`modifyInputParams()`](https://emf-creaf.github.io/medfate/reference/modifyParams.md)
  : Modify parameters
- [`trait_family_means`](https://emf-creaf.github.io/medfate/reference/TaxonTraitMeans.md)
  : Parameter average values

## Forest structure

Summarizing and displaying forest attributes

- [`summary(`*`<forest>`*`)`](https://emf-creaf.github.io/medfate/reference/summary.forest.md)
  [`print(`*`<summary.forest>`*`)`](https://emf-creaf.github.io/medfate/reference/summary.forest.md)
  : Summary of forest structure
- [`plot(`*`<forest>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.forest.md)
  [`shinyplot(`*`<forest>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.forest.md)
  : Plot forest attributes

## Meteorology

Meteorological forcing

- [`examplemeteo`](https://emf-creaf.github.io/medfate/reference/examplemeteo.md)
  : Example daily meteorology data

## Soil

Soil initialization

- [`defaultSoilParams()`](https://emf-creaf.github.io/medfate/reference/defaultSoilParams.md)
  : Default soil parameters
- [`soil()`](https://emf-creaf.github.io/medfate/reference/soil.md)
  [`summary(`*`<soil>`*`)`](https://emf-creaf.github.io/medfate/reference/soil.md)
  : Soil initialization
- [`soil_redefineLayers()`](https://emf-creaf.github.io/medfate/reference/soil_redefineLayers.md)
  : Redefine soil layer widths
- [`soil_retentionCurvePlot()`](https://emf-creaf.github.io/medfate/reference/soil_retentionCurvePlot.md)
  [`soil_conductivityCurvePlot()`](https://emf-creaf.github.io/medfate/reference/soil_retentionCurvePlot.md)
  : Soil water retention and conductivity plots

## Simulation inputs

Simulation control and input objects

- [`spwbInput()`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  [`growthInput()`](https://emf-creaf.github.io/medfate/reference/modelInput.md)
  : Input for simulation models
- [`defaultControl()`](https://emf-creaf.github.io/medfate/reference/defaultControl.md)
  : Control parameters for simulation models
- [`defaultManagementFunction()`](https://emf-creaf.github.io/medfate/reference/defaultManagementFunction.md)
  [`defaultManagementArguments()`](https://emf-creaf.github.io/medfate/reference/defaultManagementFunction.md)
  : Default forest management actions
- [`resetInputs()`](https://emf-creaf.github.io/medfate/reference/resetInputs.md)
  : Reset simulation inputs

## Simulation functions

Simulation model functions

- [`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md) :
  Soil-plant water balance
- [`spwb_day()`](https://emf-creaf.github.io/medfate/reference/spwb_day.md)
  : Single-day soil-plant water balance
- [`growth()`](https://emf-creaf.github.io/medfate/reference/growth.md)
  : Forest growth
- [`growth_day()`](https://emf-creaf.github.io/medfate/reference/growth_day.md)
  : Single-day forest growth
- [`fordyn()`](https://emf-creaf.github.io/medfate/reference/fordyn.md)
  : Forest dynamics

## Plots and summaries

Summaries, extraction and plots of simulation results

- [`plot(`*`<spwb>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md)
  [`plot(`*`<aspwb>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md)
  [`plot(`*`<pwb>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md)
  [`plot(`*`<growth>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md)
  [`plot(`*`<fordyn>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md)
  : Plots simulation results
- [`plot(`*`<spwb_day>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md)
  [`plot(`*`<growth_day>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md)
  [`plot(`*`<pwb_day>`*`)`](https://emf-creaf.github.io/medfate/reference/plot.spwb_day.md)
  : Plots simulation results for one day
- [`summary(`*`<spwb>`*`)`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md)
  [`summary(`*`<aspwb>`*`)`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md)
  [`summary(`*`<pwb>`*`)`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md)
  [`summary(`*`<growth>`*`)`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md)
  [`summary(`*`<fordyn>`*`)`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md)
  : Summarize simulation results
- [`extract()`](https://emf-creaf.github.io/medfate/reference/extract.md)
  : Extracts model outputs
- [`shinyplot()`](https://emf-creaf.github.io/medfate/reference/shinyplot.md)
  : Shiny app with interactive plots

## Post-processing

Other post-processing functions

- [`resistances()`](https://emf-creaf.github.io/medfate/reference/resistances.md)
  : Soil-plant resistances
- [`droughtStress()`](https://emf-creaf.github.io/medfate/reference/droughtStress.md)
  : Drought stress indicators
- [`waterUseEfficiency()`](https://emf-creaf.github.io/medfate/reference/waterUseEfficiency.md)
  : Water use efficiency
- [`fireHazard()`](https://emf-creaf.github.io/medfate/reference/fireHazard.md)
  : Fire hazard

## Model analysis

Evaluation and optimization

- [`exampleobs`](https://emf-creaf.github.io/medfate/reference/exampleobs.md)
  : Example observed data
- [`evaluation_table()`](https://emf-creaf.github.io/medfate/reference/evaluation.md)
  [`evaluation_stats()`](https://emf-creaf.github.io/medfate/reference/evaluation.md)
  [`evaluation_plot()`](https://emf-creaf.github.io/medfate/reference/evaluation.md)
  [`evaluation_metric()`](https://emf-creaf.github.io/medfate/reference/evaluation.md)
  : Evaluation of simulations results
- [`multiple_runs()`](https://emf-creaf.github.io/medfate/reference/optimization.md)
  [`optimization_function()`](https://emf-creaf.github.io/medfate/reference/optimization.md)
  [`optimization_evaluation_function()`](https://emf-creaf.github.io/medfate/reference/optimization.md)
  [`optimization_multicohort_function()`](https://emf-creaf.github.io/medfate/reference/optimization.md)
  [`optimization_evaluation_multicohort_function()`](https://emf-creaf.github.io/medfate/reference/optimization.md)
  : Multiple model runs and function factories for optimization routines
- [`utils_rockOptimization()`](https://emf-creaf.github.io/medfate/reference/utils_rockOptimization.md)
  : Optimization of rock fragment content
- [`utils_ldrExploration()`](https://emf-creaf.github.io/medfate/reference/utils_ldrOptimization.md)
  [`utils_ldrOptimization()`](https://emf-creaf.github.io/medfate/reference/utils_ldrOptimization.md)
  : Optimization of root distribution

## Wildfire

Fuel characteristics and fire behaviour

- [`fire_FCCS()`](https://emf-creaf.github.io/medfate/reference/fire_behaviour.md)
  [`fire_Rothermel()`](https://emf-creaf.github.io/medfate/reference/fire_behaviour.md)
  : Fire behaviour functions
- [`fuel_stratification()`](https://emf-creaf.github.io/medfate/reference/fuel_properties.md)
  [`fuel_FCCS()`](https://emf-creaf.github.io/medfate/reference/fuel_properties.md)
  : Fuel stratification and fuel characteristics
- [`SFM_metric`](https://emf-creaf.github.io/medfate/reference/SFM_metric.md)
  : Standard fuel models (Albini 1976, Scott & Burgan 2005)
