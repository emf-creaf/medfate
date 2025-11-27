# Plant water balance

Function `pwb()` performs plant water balance only (i.e. soil moisture
dynamics is an input) at daily steps for a given forest stand during a
period specified in the input climatic data. It works much as
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md) but
imposing soil moisture dynamics. Plant transpiration and photosynthesis
processes are conducted with different level of detail depending on the
transpiration mode.

## Usage

``` r
pwb(
  x,
  meteo,
  W,
  latitude,
  elevation,
  slope = NA_real_,
  aspect = NA_real_,
  canopyEvaporation = numeric(0),
  snowMelt = numeric(0),
  soilEvaporation = numeric(0),
  herbTranspiration = numeric(0),
  CO2ByYear = numeric(0)
)
```

## Arguments

- x:

  An object of class
  [`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md).

- meteo:

  A data frame with daily meteorological data series (see
  [`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md)).

- W:

  A matrix with the same number of rows as `meteo` and as many columns
  as soil layers, containing the soil moisture of each layer as
  proportion of field capacity.

- latitude:

  Latitude (in degrees).

- elevation, slope, aspect:

  Elevation above sea level (in m), slope (in degrees) and aspect (in
  degrees from North).

- canopyEvaporation:

  A vector of daily canopy evaporation (from interception) values (mm).
  The length should match the number of rows in `meteo`.

- snowMelt:

  A vector of daily snow melt values (mm). The length should match the
  number of rows in `meteo`.

- soilEvaporation:

  A vector of daily bare soil evaporation values (mm). The length should
  match the number of rows in `meteo`.

- herbTranspiration:

  A vector of daily herbaceous transpiration values (mm). The length
  should match the number of rows in `meteo`.

- CO2ByYear:

  A named numeric vector with years as names and atmospheric CO2
  concentration (in ppm) as values. Used to specify annual changes in
  CO2 concentration along the simulation (as an alternative to
  specifying daily values in `meteo`).

## Value

A list of class 'pwb' with the same elements as explained in
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md).

## References

De Cáceres M, Martínez-Vilalta J, Coll L, Llorens P, Casals P, Poyatos
R, Pausas JG, Brotons L. (2015) Coupling a water balance model with
forest inventory data to predict drought stress: the role of forest
structural changes vs. climate changes. Agricultural and Forest
Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).

De Cáceres M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L,
Poyatos R, Cabon A, Granda V, Forner A, Valladares F, Martínez-Vilalta J
(2021) Unravelling the effect of species mixing on water use and drought
stress in holm oak forests: a modelling approach. Agricultural and
Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).

Granier A, Bréda N, Biron P, Villette S (1999) A lumped water balance
model to evaluate duration and intensity of drought constraints in
forest stands. Ecol Modell 116:269–283.
https://doi.org/10.1016/S0304-3800(98)00205-1.

Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022)
SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations
of plant water status and drought-induced mortality at the ecosystem
level. Geoscientific Model Development 15, 5593-5626
(doi:10.5194/gmd-15-5593-2022).

Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S.
Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to
the environment from the optimization of photosynthetic gain and
hydraulic cost. Plant Cell and Environment 40, 816-830 (doi:
10.1111/pce.12852).

## See also

[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md),
[`spwb_day`](https://emf-creaf.github.io/medfate/reference/spwb_day.md),
[`plot.spwb`](https://emf-creaf.github.io/medfate/reference/plot.spwb.md),
[`extract`](https://emf-creaf.github.io/medfate/reference/extract.md),
[`summary.spwb`](https://emf-creaf.github.io/medfate/reference/summary.spwb.md),
[`forest`](https://emf-creaf.github.io/medfate/reference/forest.md),
[`spwb`](https://emf-creaf.github.io/medfate/reference/spwb.md),
[`aspwb`](https://emf-creaf.github.io/medfate/reference/aspwb.md)

## Author

- Miquel De Cáceres Ainsa, CREAF

- Nicolas Martin-StPaul, URFM-INRAE
