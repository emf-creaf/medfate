# Standard fuel models (Albini 1976, Scott & Burgan 2005)

Standard fuel models converted to metric system. Copied from package
'Rothermel' (Giorgio Vacchiano, Davide Ascoli).

## Format

A data frame including standard fuel models as in Albini (1976) and
Scott and Burgan (2005), to be used as input of
[`fire_Rothermel`](https://emf-creaf.github.io/medfate/reference/fire_behaviour.md)
function. All values converted to metric format.

- `Fuel_Model_Type`:

  A factor with levels `D` (for dynamic) or `S` (for static).

- `Load_1h`:

  Loading of 1h fuel class \[t/ha\].

- `Load_10h`:

  Loading of 10h fuel class \[t/ha\].

- `Load_100h`:

  Loading of 100h fuel class \[t/ha\]

- `Load_Live_Herb`:

  Loading of herbaceous fuels \[t/ha\]

- `Load_Live_Woody`:

  Loading of woody fuels \[t/ha\]

- `SA/V_1h`:

  Surface area to volume ratio of 1h fuel class \[m2/m3\]

- `SA/V_10h`:

  Surface area to volume ratio of 10h fuel class \[m2/m3\]

- `SA/V_100h`:

  Surface area to volume ratio of 100h fuel class \[m2/m3\]

- `SA/V_Live_Herb`:

  Surface area to volume ratio of herbaceous fuels \[m2/m3\]

- `SA/V_Live_Woody`:

  Surface area to volume ratio of woody fuels \[m2/m3\]

- `Fuel_Bed_Depth`:

  Fuel bed depth \[cm\]

- `Mx_dead`:

  Dead fuel moisture of extinction \[percent\]

- `Heat_1h`:

  Heat content of 1h fuel class \[kJ/kg\]

- `Heat_10h`:

  Heat content of 10h fuel class \[kJ/kg\]

- `Heat_100h`:

  Heat content of 100h fuel class \[kJ/kg\]

- `Heat_Live_Herb`:

  Heat content of herbaceous fuels \[kJ/kg\]

- `Heat_Live_Woody`:

  Heat content of woody fuels \[kJ/kg\]

## Source

Albini, F. A. (1976). Computer-based models of wildland fire behavior: A
users' manual. Ogden, UT: US Department of Agriculture, Forest Service,
Intermountain Forest and Range Experiment Station.

Scott, J., and Burgan, R. E. (2005). A new set of standard fire behavior
fuel models for use with Rothermel's surface fire spread model. Gen.
Tech. Rep. RMRSGTR-153. Fort Collins, CO: US Department of Agriculture,
Forest Service, Rocky Mountain Research Station.

## See also

[`fire_Rothermel`](https://emf-creaf.github.io/medfate/reference/fire_behaviour.md)

## Examples

``` r
data(SFM_metric)
```
