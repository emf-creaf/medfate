# Plots simulation results for one day

Functions to plot the sub-daily simulation results of
[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
[`growth_day`](https://emf-creaf.github.io/medfate/reference/growth_day.md)
or the transpiration calculations of
[`transp_transpirationSperry`](https://emf-creaf.github.io/medfate/reference/transp_modes.md)
or
[`transp_transpirationSureau`](https://emf-creaf.github.io/medfate/reference/transp_modes.md).

## Usage

``` r
# S3 method for class 'spwb_day'
plot(
  x,
  type = "PlantTranspiration",
  bySpecies = FALSE,
  xlim = NULL,
  ylim = NULL,
  xlab = NULL,
  ylab = NULL,
  ...
)

# S3 method for class 'growth_day'
plot(
  x,
  type = "PlantTranspiration",
  bySpecies = FALSE,
  xlim = NULL,
  ylim = NULL,
  xlab = NULL,
  ylab = NULL,
  ...
)

# S3 method for class 'pwb_day'
plot(
  x,
  type = "PlantTranspiration",
  bySpecies = FALSE,
  xlim = NULL,
  ylim = NULL,
  xlab = NULL,
  ylab = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `spwb_day`, `growth_day` or `pwb_day`.

- type:

  The information to be plotted (see details).

- bySpecies:

  Allows aggregating output by species, before drawing plots.
  Aggregation can involve a sum (as for plant LAI or transpiration) or a
  LAI-weighted mean (as for plant stress or plant water potential).

- xlim:

  Range of values for x.

- ylim:

  Range of values for y.

- xlab:

  x-axis label.

- ylab:

  y-axis label.

- ...:

  Additional parameters for function `plot`.

## Value

An ggplot object

## Details

The following plots are available for `spwb_day` and `pwb_day`:

- `"LeafPsi"`: Leaf water potential (for shade and sunlit leaves).

- `"LeafPsiAverage"`: Average leaf water potential.

- `"RootPsi"`: Root crown water potential.

- `"StemPsi"`: Stem water potential.

- `"StemPLC"`: (Average) percentage of loss conductance in the stem
  conduits.

- `"StemRWC"`: (Average) relative water content in the stem.

- `"LeafRWC"`: Relative water content in the leaf.

- `"StemSympRWC"`: (Average) relative water content in the stem
  symplasm.

- `"LeafSympRWC"`: Relative water content in the leaf symplasm.

- `"SoilPlantConductance"`: Overall soil plant conductance (calculated
  as the derivative of the supply function).

- `"PlantExtraction"`: Water extracted from each soil layer.

- `"PlantTranspiration"`: Plant cohort transpiration per ground area.

- `"TranspirationPerLeaf"`: Plant cohort transpiration per leaf area.

- `"PlantGrossPhotosynthesis"`: Plant cohort gross photosynthesis per
  ground area.

- `"GrossPhotosynthesisPerLeaf"`: Plant cohort gross photosynthesis per
  leaf area.

- `"PlantNetPhotosynthesis"`: Plant cohort net photosynthesis per ground
  area.

- `"NetPhotosynthesisPerLeaf"`: Plant cohort net photosynthesis per leaf
  area.

- `"LeafTranspiration"`: Instantaneous transpiration per leaf area
  (differentiates sunlit and shade leaves).

- `"LeafGrossPhotosynthesis"`: Instantaneous gross photosynthesis per
  leaf area (differentiates sunlit and shade leaves).

- `"LeafNetPhotosynthesis"`: Instantaneous net photosynthesis per leaf
  area (differentiates sunlit and shade leaves).

- `"LeafAbsorbedSWR"`: Absorbed short wave radiation per leaf area
  (differentiates sunlit and shade leaves).

- `"LeafAbsorbedPAR"`: Absorbed photosynthetically-active radiation per
  leaf area (differentiates sunlit and shade leaves).

- `"LeafNetLWR"`: Net long wave radiation per leaf area (differentiates
  sunlit and shade leaves).

- `"LeafCi"`: Leaf intercellular CO2 concentration (differentiates
  sunlit and shade leaves).

- `"LeafIntrinsicWUE"`: Leaf intrinsic water use efficiency, i.e. the
  ratio between instantaneous photosynthesis and stomatal conductance
  (differentiates sunlit and shade leaves).

- `"LeafVPD"`: Leaf vapour pressure deficit (differentiates sunlit and
  shade leaves).

- `"LeafStomatalConductance"`: Leaf stomatal conductance to water vapour
  (differentiates sunlit and shade leaves).

- `"LeafTemperature"`: Leaf temperature (differentiates sunlit and shade
  leaves).

- `"Temperature"`: Above-canopy, inside-canopy and soil temperature.

- `"CanopyEnergyBalance"`: Canopy energy balance components.

- `"SoilEnergyBalance"`: Soil energy balance components.

- `"PlantWaterBalance"`: Difference between water extraction from the
  soil and transpired water per ground area.

- `"WaterBalancePerLeaf"`: Difference between water extraction from the
  soil and transpired water per leaf area.

And the following plots are additionally available for `growth_day`:

- `"GrossPhotosynthesis"`: Gross photosynthesis rate per dry weight.

- `"MaintenanceRespiration"`: Maintenance respiration cost per dry
  weight.

- `"RootExudation"`: Root exudation rate per dry weight.

- `"LabileCarbonBalance"`: Labile carbon balance per dry weight.

- `"SugarLeaf"`: Sugar concentration in leaves.

- `"StarchLeaf"`: Starch concentration in leaves.

- `"SugarSapwood"`: Sugar concentration in sapwood.

- `"StarchSapwood"`: Starch concentration in sapwood.

- `"SugarTransport"`: Phloem sugar transport rate.

## Note

Only for soil plant water balance simulations using
`transpirationMode = "Sperry"` or `transpirationMode = "Sureau"`. This
function can be used to display subdaily dynamics of corresponding to
single days on
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) runs, if
control option `subdailyResults` is set to `TRUE`. See also option
`subdaily` in
[`plot.spwb`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md).

## See also

[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
[`growth_day`](https://emf-creaf.github.io/medfate/reference/growth_day.md),
[`plot.spwb`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md)

## Author

Miquel De Cáceres Ainsa, CREAF

## Examples

``` r
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Define soil with default soil params (2 layers)
examplesoil <- defaultSoilParams(4)

#Switch to 'Sperry' transpiration mode
control <- defaultControl("Sperry")

#Enable subdaily results for plotting
control["subdailyResults"] <- TRUE

#Initialize model inputs
x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)

#Simulate one day only
d <- 100
date <- examplemeteo$dates[d]
meteovec <- unlist(examplemeteo[d,])
sd2 <- spwb_day(x2, date, meteovec,
              latitude = 41.82592, elevation = 100, slope= 0, aspect = 0)

#Display transpiration for subdaily steps
plot(sd2, "PlantTranspiration")

```
