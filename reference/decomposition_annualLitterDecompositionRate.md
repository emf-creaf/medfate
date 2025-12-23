# Low-level decomposition functions

Functions related to litter and soil carbon decomposition processes

## Usage

``` r
decomposition_annualLitterDecompositionRate(AET, lignin)

decomposition_snagFallProbability(DBH, decayClass, durabilityEffect = 0)

decomposition_litterMetabolicFraction(ligninPercent, Nmass)

decomposition_pHEffect(x, pool)

decomposition_moistureEffect(sand, clay, soilMoisture)

decomposition_temperatureEffect(soilTemperature)
```

## Arguments

- AET:

  Actual evapotranspiration (mm)

- lignin:

  Lignin percent

- DBH:

  Diameter at breast height

- decayClass:

  Decay class, from 1 to 5

- durabilityEffect:

  Effect of wood durability

- ligninPercent:

  lignin content (% of dry)

- Nmass:

  nitrogen content (mg N / g dry)

- x:

  Soil water pH (0-14)

- pool:

  String indicating the decomposition pool

- sand, clay:

  Soil texture values in percent volume.

- soilMoisture:

  Soil moisture content, relative to saturation.

- soilTemperature:

  Soil temperature (in Celsius).

## Value

Functions `decomposition_moistureEffect`, `decomposition_pHEffect` and
`decomposition_temperatureEffect` return a scalar value representing a
factor that should modify a decomposition rate. Function
`decomposition_annualLitterDecompositionRate` directly returns a scalar
value with the annual decomposition rate (yr-1). Function
`decomposition_litterMetabolicFraction` returns a scalar with the
fraction of litter that corresponds to metabolic carbon.

## Details

Function `decomposition_moistureEffect` follows Kelly et al. (2000)
Function `decomposition_snagFallProbability` follows Vanderwell et al.
(2006)

## References

Bonan, G. (2019). Climate change and terrestrial ecosystem modeling.
Cambridge University Press, Cambridge, UK.

Vanderwel et al. (2006) Snag dynamics in partially harvested and
unmanaged northern hardwood forests. Canadian Journal of Forest Research
36: 2769-2779.

Meentemeyer (1978)

Kelly et al (2000)

## See also

[`decomposition_DAYCENT`](https://emf-creaf.github.io/medfate/reference/decomposition_DAYCENT.md)

## Author

Miquel De Cáceres Ainsa, CREAF
