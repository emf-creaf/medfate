# Model evaluation in experimental plots

## Introduction

This document presents **medfate** (**ver. 4.9.0**) model evaluation
results at stand-level, using data from a set of **20 experimental
forest plots**. The main source of observed data are SAPFLUXNET database
([Poyatos et
al. 2021](https://essd.copernicus.org/articles/13/2607/2021/)) and
FLUXNET 2015 dataset ([Pastorello et
al. 2020](https://doi.org/10.1038/s41597-020-0534-3)).

### List of sites

The table below lists the experimental forest plots used in the report
and the data sources available.

| Country | Plot | Stand | SAPFLUXNET | FLUXNET/ICOS |
|:---|:---|:---|:---|:---|
| Australia | Wombat | Mixed eucalyptus forest | AUS_WOM | AU-Wom |
| Australia | Euc-FACE | Eucalyptus trees in ambient (control) plots of a CO2 enrichment experiment | AUS_RIC_EUC_ELE | AU-Cum |
| Denmark | Soroe | European beech forest |  | DK-Sor |
| France | Puéchabon | Dense evergreen forest dominated by Q. ilex | FRA_PUE | FR-Pue |
| France | Hesse | Naturally regenerated, managed beech forest | FRA_HES_HE2_NON | FR-Hes |
| France | Fontainebleau-Barbeau | Mixed deciduous forest | FRA_FON | FR-Fon |
| France | Font-Blanche | Mixed forest with P. halepensis and Q. ilex |  | FR-Fbn |
| Italy | Collelongo | European beech forest |  | IT-Col |
| Portugal | Mitra II | Evergreen forest dominated by Quercus ilex subsp. rotundifolia | PRT_MIT | PT-Mi1 |
| Spain | Rinconada | Young, homogeneous, Quercus pyrenaica regrowth forest | ESP_RIN |  |
| Spain | Vallcebre (Cal Barrol) | Semi-deciduous sub-Mediterranean oak forest | ESP_VAL_BAR |  |
| Spain | Vallcebre (Cal Sort) | Pinus sylvestris forest in a terraced area | ESP_VAL_SOR |  |
| Spain | Prades (Tillar valley) | Mixed forest with P. sylvestris (overstory) Q. ilex (midstory) | ESP_TIL_MIX |  |
| Spain | Can Balasc | Mixed forest dominated by Q. ilex | ESP_CAN |  |
| Spain | Alto Tajo (Armallones) | Sparse mixed forest dominated by three species | ESP_ALT_ARM |  |
| Spain | Ronda (Pilones) | Mixed gimnosperm forest dominated by Abies pinsapo | ESP_RON_PIL |  |
| Switzerland | Davos Seehornwald | Subalpine coniferous (spruce) forest | CHE_DAV_SEE | CH-Dav |
| Switzerland | Lötschental | Mixed evergreen Norway spruce and deciduous European larch forest | CHE_LOT_NOR |  |
| USA | Morgan-Mornoe | Mixed temperate forest | USA_MOR_SF | US-MMS |
| USA | Sevilleta | Mixed pine-juniper forest | USA_PJS_P04_AMB |  |

### Parametrization and simulations

Forest water balance simulations (i.e. function
[`spwb()`](https://emf-creaf.github.io/medfate/reference/spwb.md)) have
been conducted using the three transpiration modes (i.e. `Granier`,
`Sperry` or `Sureau`).

The set of control parameters modified from defaults in simulations are
the following:

| transpirationMode | soilDomains | stemCavitationRecovery | leafCavitationRecovery | segmentedXylemVulnerability | subdailyResults |
|:---|:---|:---|:---|:---|:---|
| Granier | dual | rate | total | NA | NA |
| Sperry | dual | rate | total | TRUE | FALSE |
| Sureau | dual | rate | rate | FALSE | FALSE |

Soil characteristics have been tuned to modulate total available water
and fit observed saturation and residual moisture values, but
calibration exercises have not been conducted. When available, however,
local leaf area to sapwood area ratios have been used. Thus, the
evaluation exercise is meant to be more or less representative of
simulations with default species-level trait data.

### Evaluation variables

The table below lists the set of predicted variables that are evaluated
and the data sources used:

| Variable | Level | Observation source | Units |
|----|----|----|----|
| Sensible heat turbulent flux | Stand | FLUXNET / ICOS | MJ/m2 |
| Latent heat turbulent flux | Stand | FLUXNET / ICOS | MJ/m2 |
| Gross primary productivity | Stand | FLUXNET / ICOS | gC/m2 |
| Soil moisture content (topsoil) | Stand | SAPFLUXNET / FLUXNET / ICOS | % vol. |
| Transpiration per leaf area | Plant | SAPFLUXNET | l/m2 |
| Predawn/midday leaf water potential | Plant | SAPFLUXNET (addition) | MPa |

### Structure of site reports

The following contains as many sections as forest stands included in the
evaluation. The following sub-sections are reported for each stand:

1.  **General information**: General information about the site,
    topography, soil and climate, as well as data sources used.
2.  **Model inputs**: Description of model inputs (vegetation, soil,
    custom species parameters and parameterization remarks).
3.  **Climate**: Graphical description of climate inputs and predicted
    soil/canopy temperatures (under Sperry).
4.  **Evaluation results**: Evaluation results are presented for
    variables with available measurements.

## Wombat

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Wombat |
| Country | Australia |
| SAPFLUXNET code | AUS_WOM |
| SAPFLUXNET contributor (affiliation) | Anne Griebel (U. of Melbourne) |
| FLUXNET/ICOS code | AU-Wom |
| FLUXNET/ICOS contributor (affiliation) | Stefan Arndt (U. of Melbourne) |
| Latitude (º) | -37.4222 |
| Longitude (º) | 144.0944 |
| Elevation (m) | 705 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material |  |
| Soil texture | Loam |
| MAT (ºC) | 10.9 |
| MAP (mm) | 1024 |
| Forest stand | Mixed eucalyptus forest |
| Stand LAI | 2.2 |
| Stand description DOI | 10.1016/j.foreco.2016.12.017 |
| Species simulated | Eucalyptus obliqua, E. radiata, E. rubida |
| Species parameter table | SpParamsAU |
| Simulation period | 2013-2015 |
| Evaluation period | 2013-2015 |

### Model inputs

#### Vegetation

| Species            | DBH | Height |      N | Z50 |  Z95 |   LAI |
|:-------------------|----:|-------:|-------:|----:|-----:|------:|
| Eucalyptus obliqua |  23 |   2200 | 712.60 | 300 | 1000 | 1.540 |
| Eucalyptus rubida  |  23 |   2200 | 213.78 | 300 | 1000 | 0.462 |
| Eucalyptus radiata |  23 |   2200 |  91.62 | 300 | 1000 | 0.198 |

#### Soil

| widths |  clay |     sand |    om |        bd |   rfc |
|-------:|------:|---------:|------:|----------:|------:|
|    100 | 24.30 | 55.13333 | 8.340 | 0.9866667 | 13.60 |
|    200 | 25.00 | 50.00000 | 6.000 | 1.1000000 | 14.00 |
|    700 | 35.25 | 44.60000 | 4.465 | 1.2100000 | 14.35 |
|   1000 | 33.10 | 45.90000 | 4.160 | 1.2800000 | 50.00 |
|   2000 | 33.10 | 45.90000 | 0.000 | 1.3500000 | 90.00 |

#### Custom traits

| Species            |   Al2As |
|:-------------------|--------:|
| Eucalyptus obliqua | 4590.63 |
| Eucalyptus rubida  | 4590.63 |
| Eucalyptus radiata | 4590.63 |

#### Custom control

||
||
||

#### Remarks

| Title           | Remark                                         |
|:----------------|:-----------------------------------------------|
| Soil            | Taken from SoilGrids                           |
| Vegetation      | No understory or secondary species considered. |
| Weather         |                                                |
| Sapflow         | Species-level Huber value used for scaling     |
| Eddy covariance | Variables taken: LE_CORR and GPP_NT_VUT_REF    |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-16-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-17-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-18-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| WOMBAT | sperry | 730 | 2.600475 | 75.25273 | 5.016844 | 145.1778 | 0.5621591 | -1.0500799 | -0.2253586 |
| WOMBAT | sureau | 730 | 2.586075 | 74.83601 | 4.768960 | 138.0045 | 0.5667960 | -0.7745536 | -0.1648134 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-38-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-39-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| WOMBAT | sperry | 730 | -3.157123 | -54.76089 | 3.374487 | 58.53110 | 0.5492935 | -0.7500186 | -0.3466214 |
| WOMBAT | sureau | 730 | -3.263691 | -56.60932 | 3.385140 | 58.71588 | 0.6969171 | -0.6589523 | -0.3508727 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-45-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-46-1.png)

#### Gross primary productivity

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| WOMBAT | granier | 730 | -2.6944993 | -48.86060 | 2.752716 | 49.91627 | 0.7225719 | -0.6685083 | -0.3158651 |
| WOMBAT | sperry | 730 | -0.8669382 | -15.72059 | 1.594796 | 28.91918 | 0.6586017 | 0.3147340 | 0.2376486 |
| WOMBAT | sureau | 730 | -1.4026171 | -25.43431 | 1.753773 | 31.80198 | 0.7285550 | 0.1788762 | 0.1616539 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-52-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-53-1.png)

#### Soil water content (SWC.1)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| WOMBAT | granier | 1035 | -0.0817882 | -38.93786 | 0.0883318 | 42.05314 | 0.3699423 | -1.499753 | -0.7253290 |
| WOMBAT | sperry | 1035 | -0.0797165 | -37.95156 | 0.0872436 | 41.53508 | 0.3575807 | -1.428781 | -0.7040743 |
| WOMBAT | sureau | 1035 | -0.0806969 | -38.41834 | 0.0876563 | 41.73156 | 0.3692498 | -1.457076 | -0.7121353 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-59-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-60-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| WOMBAT | T1_7431 | granier | 1035 | -0.1687791 | -41.41273 | 0.1757917 | 43.13338 | 0.8809295 | 0.2753780 | 0.1776524 |
| WOMBAT | T1_7431 | sperry | 1035 | -0.1394164 | -34.20811 | 0.1455198 | 35.70568 | 0.9290043 | 0.5589456 | 0.3192633 |
| WOMBAT | T1_7431 | sureau | 1035 | -0.1385786 | -34.00254 | 0.1580950 | 38.79122 | 0.9145502 | 0.4969340 | 0.2604368 |
| WOMBAT | T2_7526 | granier | 1035 | -0.1655800 | -36.58994 | 0.1837911 | 40.61424 | 0.8687882 | 0.4248388 | 0.3363398 |
| WOMBAT | T2_7526 | sperry | 1035 | -0.1707542 | -37.73334 | 0.1730738 | 38.24591 | 0.9442202 | 0.5848070 | 0.3750396 |
| WOMBAT | T2_7526 | sureau | 1035 | -0.1361686 | -30.09060 | 0.1509832 | 33.36434 | 0.9374027 | 0.6947895 | 0.4548073 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-66-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-66-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-66-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-66-4.png)

## Euc-FACE

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Euc-FACE |
| Country | Australia |
| SAPFLUXNET code | AUS_RIC_EUC_ELE |
| SAPFLUXNET contributor (affiliation) | Teresa Gimeno (CREAF) |
| FLUXNET/ICOS code | AU-Cum |
| FLUXNET/ICOS contributor (affiliation) | Elise Pendall (U. Western Sidney) |
| Latitude (º) | -33.61778 |
| Longitude (º) | 150.74028 |
| Elevation (m) | 23 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material |  |
| Soil texture | Sandy loam |
| MAT (ºC) | 17.6 |
| MAP (mm) | 899 |
| Forest stand | Eucalyptus trees in ambient (control) plots of a CO2 enrichment experiment |
| Stand LAI | 2 |
| Stand description DOI | 10.1111/1365-2435.12532 |
| Species simulated | Eucalyptus tereticornis |
| Species parameter table | SpParamsAU |
| Simulation period | 2012-2014 |
| Evaluation period | 2012-2014 |

### Model inputs

#### Vegetation

| Species                 | DBH | Height |   N | Z50 |  Z95 | LAI | Cover |
|:------------------------|----:|-------:|----:|----:|-----:|----:|------:|
| Eucalyptus tereticornis |  21 |   2200 | 800 | 200 | 3000 |   2 |    NA |
| Herbaceous layer        |  NA |     10 |  NA |  NA |   NA |  NA |    10 |

#### Soil

| widths |     clay |     sand |        om |       bd |       rfc | VG_theta_res | VG_theta_sat |
|-------:|---------:|---------:|----------:|---------:|----------:|-------------:|-------------:|
|    300 | 18.16667 | 61.86667 | 1.8700000 | 1.246667 |  6.566667 |         0.03 |          0.4 |
|    500 | 31.00000 | 52.46667 | 0.6566667 | 1.313333 | 20.000000 |         0.03 |          0.4 |
|    500 | 31.10000 | 52.73750 | 0.6912500 | 1.331250 | 20.000000 |         0.03 |          0.4 |
|   2700 | 30.60000 | 53.30000 | 0.5600000 | 1.350000 | 90.000000 |         0.03 |          0.4 |

#### Custom traits

| Species                 | Vmax298 | Jmax298 |    Al2As |
|:------------------------|--------:|--------:|---------:|
| Eucalyptus tereticornis |      91 |     159 | 6896.552 |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Taken from SoilGrids |
| Vegetation | No understory or secondary species considered. 10% Herbaceous cover |
| Weather | CO2 set to 390 ppm |
| Sapflow | Species-level Huber value used for scaling |
| Eddy covariance | Variables taken: LE_CORR and GPP_NT_VUT_REF |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-77-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-78-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-79-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| EUCFACE | sperry | 884 | -2.110592 | -43.68872 | 4.305725 | 89.12742 | 0.3106459 | -1.0354596 | -0.3672240 |
| EUCFACE | sureau | 884 | -1.883877 | -38.99578 | 4.124443 | 85.37494 | 0.2865083 | -0.8663871 | -0.3096606 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-99-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-100-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| EUCFACE | sperry | 884 | -2.816318 | -50.52309 | 2.885080 | 51.75664 | 0.6093248 | -1.700928 | -0.9232982 |
| EUCFACE | sureau | 884 | -3.162415 | -56.73186 | 3.204562 | 57.48795 | 0.5183815 | -2.366945 | -1.1362760 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-106-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-107-1.png)

#### Gross primary productivity

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| EUCFACE | granier | 884 | -0.6786855 | -22.34654 | 1.041018 | 34.27678 | 0.3164765 | -0.4680075 | -0.2074953 |
| EUCFACE | sperry | 884 | 2.4767438 | 81.54979 | 2.562458 | 84.37204 | 0.3481429 | -6.0762564 | -1.9722406 |
| EUCFACE | sureau | 884 | 1.5056873 | 49.57658 | 1.655174 | 54.49860 | 0.2956801 | -2.2537359 | -0.9198653 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-113-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-114-1.png)

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| EUCFACE | granier | 710 | -0.0057752 | -6.057183 | 0.0473376 | 49.64857 | 0.4868848 | 0.1701432 | 0.1335998 |
| EUCFACE | sperry | 710 | -0.0096740 | -10.146286 | 0.0454911 | 47.71190 | 0.4798745 | 0.1709765 | 0.1673959 |
| EUCFACE | sureau | 710 | -0.0078969 | -8.282440 | 0.0463297 | 48.59143 | 0.4635703 | 0.1651038 | 0.1520477 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-120-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-121-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| EUCFACE | T1_7590 | granier | 708 | -0.3155961 | -48.56817 | 0.3260133 | 50.17131 | 0.7267545 | -0.7770994 | -0.4622488 |
| EUCFACE | T1_7590 | sperry | 708 | -0.1699725 | -26.15765 | 0.2011028 | 30.94840 | 0.7715523 | 0.1811446 | 0.0980053 |
| EUCFACE | T1_7590 | sureau | 708 | -0.2231318 | -34.33852 | 0.2581643 | 39.72979 | 0.7113742 | -0.2364998 | -0.1579294 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-127-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-127-2.png)

#### Leaf water potential

| Site | Cohort | WP | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| EUCFACE | T1_7590 | Midday | sperry | 15 | 0.6268778 | 31.47533 | 0.7376916 | 37.03925 | 0.2071056 | -1.5067114 | -0.5923862 |
| EUCFACE | T1_7590 | Midday | sureau | 15 | 1.0943407 | 54.94649 | 1.1532759 | 57.90561 | 0.1509508 | -4.6859314 | -1.4894695 |
| EUCFACE | T1_7590 | Predawn | sperry | 8 | 0.1375796 | 21.37955 | 0.2778156 | 43.17188 | 0.3955593 | -1.3943445 | -0.8812367 |
| EUCFACE | T1_7590 | Predawn | sureau | 8 | 0.0764864 | 11.88580 | 0.2165934 | 33.65810 | 0.5136820 | -0.5813806 | -0.4666687 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-132-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-132-2.png)

## Soroe

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Soroe |
| Country | Denmark |
| SAPFLUXNET code |  |
| SAPFLUXNET contributor (affiliation) |  |
| FLUXNET/ICOS code | DK-Sor |
| FLUXNET/ICOS contributor (affiliation) | Andreas Ibrom (Technical University of Denmark) |
| Latitude (º) | 11.6446 |
| Longitude (º) | 55.4859 |
| Elevation (m) | 40 |
| Slope (º) | 0 |
| Aspect (º) | NA |
| Parent material |  |
| Soil texture |  |
| MAT (ºC) | 8.2 |
| MAP (mm) | 660 |
| Forest stand | European beech forest |
| Stand LAI | 4.5 |
| Stand description DOI | 10.1016/j.agrformet.2011.02.013 |
| Species simulated | Fagus sylvatica |
| Species parameter table | SpParamsFR |
| Simulation period | 2003-2006 |
| Evaluation period | 2003-2006 |

### Model inputs

#### Vegetation

| Species         |      DBH | Height |       N | Z50 | Z95 | LAI |
|:----------------|---------:|-------:|--------:|:----|:----|----:|
| Fagus sylvatica | 24.28232 |   2099 | 352.987 | NA  | NA  | 4.5 |

#### Soil

| widths |  clay |  sand |   om |   bd |  rfc |
|-------:|------:|------:|-----:|-----:|-----:|
|    400 | 15.27 | 58.73 | 5.42 | 1.28 | 7.48 |
|    200 | 19.20 | 55.70 | 2.18 | 1.58 | 6.10 |
|   3000 | 19.60 | 55.45 | 1.78 | 1.65 | 7.55 |

#### Custom traits

| Species         | Vmax298 | Jmax298 |
|:----------------|--------:|--------:|
| Fagus sylvatica |    94.5 |   159.9 |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Taken from SoilGrids |
| Vegetation | No understory |
| Weather | From V. saponaro |
| Sapflow | Not available |
| Soil moisture | Taken from FLUXNET |
| Eddy covariance | Variables H_F_MDS and LE_F_MDS for sensible and latent heat |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-142-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-143-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-144-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| SOROE | sperry | 1461 | -2.821307 | -1810.811 | 6.351118 | 4076.365 | 0.3466944 | -1.846792 | -0.7945913 |
| SOROE | sureau | 1461 | -1.802605 | -1156.974 | 6.127848 | 3933.063 | 0.3904302 | -1.735913 | -0.7315035 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-164-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-165-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| SOROE | sperry | 1461 | 1.029663 | 35.618677 | 1.804545 | 62.42379 | 0.6811450 | 0.3244145 | 0.3215193 |
| SOROE | sureau | 1461 | -0.100055 | -3.461158 | 2.407022 | 83.26501 | 0.3427894 | -0.0285535 | 0.0949973 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-171-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-172-1.png)

#### Gross primary productivity

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| SOROE | granier | 1461 | 0.6686326 | 12.20259 | 2.772429 | 50.59701 | 0.8347640 | 0.6514127 | 0.4559557 |
| SOROE | sperry | 1461 | 2.3387976 | 42.68320 | 2.897244 | 52.87489 | 0.8810292 | 0.5980873 | 0.4314628 |
| SOROE | sureau | 1461 | 1.2300464 | 22.44842 | 2.528592 | 46.14696 | 0.8677881 | 0.6954012 | 0.5038048 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-178-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-179-1.png)

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| SOROE | granier | 1461 | -0.0717646 | -31.52245 | 0.0772213 | 33.91927 | 0.5606022 | -1.3741448 | -0.6940781 |
| SOROE | sperry | 1461 | -0.0649644 | -28.53545 | 0.0710859 | 31.22433 | 0.5277637 | -1.0820531 | -0.5594810 |
| SOROE | sureau | 1461 | -0.0552180 | -24.25439 | 0.0661443 | 29.05374 | 0.4436313 | -0.7859524 | -0.4510720 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-185-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-186-1.png)

## Puéchabon

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Puéchabon |
| Country | France |
| SAPFLUXNET code | FRA_PUE |
| SAPFLUXNET contributor (affiliation) | Jean-Marc Limousin (CEFE-CNRS) |
| FLUXNET/ICOS code | FR-Pue |
| FLUXNET/ICOS contributor (affiliation) | Jean-Marc Limousin (CEFE-CNRS) |
| Latitude (º) | 43.74 |
| Longitude (º) | 3.6 |
| Elevation (m) | 270 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material | Limestone |
| Soil texture | Silty clay loam |
| MAT (ºC) | 13.4 |
| MAP (mm) | 720 |
| Forest stand | Dense evergreen forest dominated by Q. ilex |
| Stand LAI | 2 |
| Stand description DOI | 10.1111/j.1365-2486.2009.01852.x |
| Species simulated | Quercus ilex, Buxus sempervirens |
| Species parameter table | SpParamsFR |
| Simulation period | 2004-2006 |
| Evaluation period | 2004-2006 |

### Model inputs

#### Vegetation

| Species            |    DBH |   Height |    N | Z50 |  Z95 | LAI | Cover |
|:-------------------|-------:|---------:|-----:|----:|-----:|----:|------:|
| Quercus ilex       | 9.1156 | 530.2222 | 1750 | 529 | 2287 | 2.0 |    NA |
| Buxus sempervirens |     NA | 200.0000 |   NA | 390 | 1470 | 0.2 |    13 |
| Herbaceous layer   |     NA |  20.0000 |   NA |  NA |   NA |  NA |    10 |

#### Soil

| widths | clay | sand |  om |   bd | rfc | VG_theta_sat | VG_theta_res |
|-------:|-----:|-----:|----:|-----:|----:|-------------:|-------------:|
|    100 |   39 |   26 |   6 | 1.45 |  75 |         0.27 |        0.015 |
|    200 |   39 |   26 |   4 | 1.45 |  75 |         0.27 |        0.015 |
|    200 |   39 |   26 |   3 | 1.45 |  75 |         0.27 |        0.015 |
|   1500 |   39 |   26 |   1 | 1.45 |  80 |         0.27 |        0.015 |
|   2500 |   39 |   26 |   1 | 1.45 |  90 |         0.27 |        0.015 |

#### Custom traits

| Species | SLA | VCleaf_P12 | VCleaf_P50 | VCleaf_P88 | VCleaf_slope | VCstem_P12 | VCstem_P50 | VCstem_P88 | VCstem_slope | VCroot_P12 | VCroot_P50 | VCroot_P88 | VCroot_slope | VCleaf_kmax | Kmax_stemxylem | LeafEPS | LeafPI0 | LeafAF | StemEPS | StemPI0 | StemAF | Gswmin | Gswmax | Gs_P50 | Gs_slope | Al2As |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Quercus ilex | 4.55 | -4.004731 | -5.25 | -6.495269 | 40 | -4.739642 | -6.4 | -8.060358 | 30 | -2.004731 | -3.25 | -4.495269 | 40 | 2.63 | 0.20 | 15 | -2.5 | 0.4 | 15 | -2.5 | 0.4 | 0.002 | 0.20 | -2.114188 | 44.70588 | 1540.671 |
| Buxus sempervirens | 5.19 | -5.004731 | -6.25 | -7.495269 | 40 | NA | NA | NA | NA | NA | NA | NA | NA | 2.00 | 0.15 | NA | NA | NA | NA | NA | NA | 0.002 | 0.18 | NA | NA | NA |

#### Custom control

||
||
||

#### Remarks

| Title      | Remark                              |
|:-----------|:------------------------------------|
| Soil       | Adjusted theta_res and theta_sat    |
| Vegetation | Using B. sempervirens as understory |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-199-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-200-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-201-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAPUE | sperry | 1096 | 0.7091794 | 20.03411 | 4.126571 | 116.5744 | 0.7012380 | -0.0710422 | 0.0910885 |
| FRAPUE | sureau | 1096 | 0.8773274 | 24.78424 | 3.792415 | 107.1346 | 0.6812086 | 0.0966535 | 0.1646892 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-221-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-222-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAPUE | sperry | 1096 | -0.6213170 | -19.52938 | 2.123809 | 66.75606 | 0.1528356 | -0.5932759 | -0.1132033 |
| FRAPUE | sureau | 1096 | -0.9235667 | -29.02976 | 1.751178 | 55.04342 | 0.5261058 | 0.0789480 | 0.0821131 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-228-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-229-1.png)

#### Gross primary productivity

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAPUE | granier | 1096 | -0.0564101 | -1.709898 | 1.199534 | 36.36019 | 0.6572415 | 0.3420145 | 0.2373708 |
| FRAPUE | sperry | 1096 | -0.7486825 | -22.694012 | 1.522117 | 46.13830 | 0.3724089 | -0.4864733 | 0.0322818 |
| FRAPUE | sureau | 1096 | -0.3031693 | -9.189647 | 0.889379 | 26.95880 | 0.7978724 | 0.6035443 | 0.4345584 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-235-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-236-1.png)

#### Soil water content (SWC.2)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAPUE | granier | 1021 | -0.0187856 | -11.79326 | 0.0293413 | 18.41991 | 0.7824730 | 0.4044871 | 0.2350578 |
| FRAPUE | sperry | 1021 | -0.0447104 | -28.06832 | 0.0480113 | 30.14061 | 0.8133482 | -0.4637887 | -0.2516798 |
| FRAPUE | sureau | 1021 | -0.0186275 | -11.69401 | 0.0269961 | 16.94766 | 0.8403371 | 0.5047975 | 0.2961974 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-242-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-243-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAPUE | T1_2854 | granier | 1096 | -0.0875490 | -25.52221 | 0.1297054 | 37.81162 | 0.7584299 | 0.4113520 | 0.2834507 |
| FRAPUE | T1_2854 | sperry | 1096 | -0.0330665 | -9.63952 | 0.2732000 | 79.64306 | 0.2378637 | -1.8985760 | -0.5092762 |
| FRAPUE | T1_2854 | sureau | 1096 | -0.0857442 | -24.99608 | 0.1350925 | 39.38208 | 0.8306098 | 0.3821038 | 0.2536897 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-249-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-249-2.png)

#### Leaf water potential

| Site | Cohort | WP | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAPUE | T1_2854 | Midday | sperry | 16 | -0.8667751 | -28.60570 | 1.1129919 | 36.73146 | 0.2206411 | -4.2696050 | -1.5256593 |
| FRAPUE | T1_2854 | Midday | sureau | 28 | 0.3260463 | 10.44945 | 0.7664665 | 24.56448 | 0.7572531 | -0.7769527 | -0.2838443 |
| FRAPUE | T1_2854 | Predawn | sperry | 16 | -1.7673364 | -147.18221 | 1.7809574 | 148.31656 | 0.4836091 | -5.0731345 | -1.3541743 |
| FRAPUE | T1_2854 | Predawn | sureau | 28 | -0.9489590 | -64.87618 | 1.0380478 | 70.96680 | 0.8682520 | -0.1800865 | -0.0916681 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-254-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-254-2.png)

## Hesse

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Hesse |
| Country | France |
| SAPFLUXNET code | FRA_HES_HE2_NON |
| SAPFLUXNET contributor (affiliation) | André Granier (INRAE) |
| FLUXNET/ICOS code | FR-Hes |
| FLUXNET/ICOS contributor (affiliation) | Matthias Cuntz (INRAE) |
| Latitude (º) | 48.6742 |
| Longitude (º) | 7.0647 |
| Elevation (m) | 300 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material |  |
| Soil texture | Silt loam |
| MAT (ºC) | 10 |
| MAP (mm) | 1003 |
| Forest stand | Naturally regenerated, managed beech forest |
| Stand LAI | 7 |
| Stand description DOI | 10.1051/forest:2008052 |
| Species simulated | Fagus sylvatica |
| Species parameter table | SpParamsFR |
| Simulation period | 2001-2003 |
| Evaluation period | 2001-2003 |

### Model inputs

#### Vegetation

| Species          |   DBH | Height |    N | Z50 |  Z95 | LAI | Cover |
|:-----------------|------:|-------:|-----:|----:|-----:|----:|------:|
| Fagus sylvatica  | 12.91 |   1300 | 3203 | 300 | 1200 |   7 |    NA |
| Herbaceous layer |    NA |     20 |   NA |  NA |   NA |  NA |     5 |

#### Soil

| widths | sand | clay |  om | rfc |   bd | VG_theta_sat |
|-------:|-----:|-----:|----:|----:|-----:|-------------:|
|    200 |    8 |   25 |   6 |   9 | 1.16 |         0.46 |
|    300 |    8 |   35 |   3 |  13 | 1.37 |         0.43 |
|    300 |    8 |   45 |   1 |  15 | 1.58 |         0.38 |
|    400 |    8 |   45 |   0 |  40 | 1.58 |         0.35 |
|   2100 |    8 |   45 |   0 |  90 | 1.58 |         0.30 |

#### Custom traits

| Species         |   Al2As |
|:----------------|--------:|
| Fagus sylvatica | 2076.12 |

#### Custom control

||
||
||

#### Remarks

| Title      | Remark                                      |
|:-----------|:--------------------------------------------|
| Soil       | VG_theta_sat modified                       |
| Vegetation | No woody understory but 5% herbaceous layer |
| Sapflow    | Scaling using species-level Huber value     |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-264-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-265-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-266-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAHES | sperry | 1070 | -6.424645 | -437.5994 | 8.715119 | 593.6096 | 0.3387258 | -6.405685 | -2.244986 |
| FRAHES | sureau | 1070 | -5.942324 | -404.7473 | 8.405016 | 572.4877 | 0.3513127 | -6.122615 | -2.129522 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-286-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-287-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAHES | sperry | 1064 | 4.924529 | 208.2519 | 4.995029 | 211.2332 | 0.5282436 | -3.427270 | -1.206155 |
| FRAHES | sureau | 1064 | 4.411126 | 186.5407 | 4.607555 | 194.8475 | 0.3372953 | -3.368079 | -1.035019 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-293-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-294-1.png)

#### Soil water content (SWC.1)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAHES | granier | 43 | -0.1167927 | -35.96654 | 0.1195657 | 36.82048 | -0.0447996 | -3.205517 | -1.0635067 |
| FRAHES | sperry | 43 | -0.1015041 | -31.25840 | 0.1044324 | 32.16018 | 0.3739676 | -2.208586 | -0.8023324 |
| FRAHES | sureau | 43 | -0.0961584 | -29.61217 | 0.1009050 | 31.07389 | 0.1508324 | -2.102857 | -0.7414540 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-301-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-302-1.png)

#### Soil water content (SWC.2)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAHES | granier | 43 | -0.0554048 | -16.16960 | 0.0626661 | 18.28876 | 0.6292911 | -0.4740954 | -0.1568336 |
| FRAHES | sperry | 43 | -0.0697226 | -20.34816 | 0.0697226 | 20.34816 | 0.9252582 | -0.4952064 | -0.2870983 |
| FRAHES | sureau | 43 | -0.0433818 | -12.66075 | 0.0507298 | 14.80521 | 0.8114872 | -0.0018787 | 0.0635142 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-307-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-308-1.png)

#### Soil water content (SWC.3)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAHES | granier | 43 | -0.0131000 | -3.908371 | 0.0206488 | 6.160557 | 0.8145394 | 0.3367829 | 0.2724901 |
| FRAHES | sperry | 43 | -0.0537498 | -16.036236 | 0.0574243 | 17.132536 | 0.8655994 | -3.8996764 | -1.0232081 |
| FRAHES | sureau | 43 | -0.0261403 | -7.798945 | 0.0333691 | 9.955679 | 0.9139930 | -0.5801261 | -0.1756818 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-313-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-314-1.png)

#### Soil water content (SWC.4)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAHES | granier | 43 | -0.0253122 | -7.811079 | 0.0305826 | 9.437468 | 0.8157366 | -8.176814 | -2.084371 |
| FRAHES | sperry | 43 | -0.0732815 | -22.613911 | 0.0766628 | 23.657348 | 0.7468941 | -56.591412 | -6.731739 |
| FRAHES | sureau | 43 | -0.0538238 | -16.609474 | 0.0581363 | 17.940258 | 0.7958251 | -33.469512 | -4.863269 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-319-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-320-1.png)

#### Soil water content (SWC.5)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAHES | granier | 43 | -0.0021158 | -0.7151926 | 0.0041871 | 1.415296 | 0.3674432 | -1.845976 | -0.7531172 |
| FRAHES | sperry | 43 | -0.0041774 | -1.4120265 | 0.0049052 | 1.658027 | 0.4283351 | -2.948199 | -1.0537856 |
| FRAHES | sureau | 43 | -0.0033983 | -1.1486742 | 0.0043343 | 1.465055 | 0.4656803 | -2.149878 | -0.8147531 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-325-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-326-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FRAHES | T1_1396 | granier | 559 | -0.2070140 | -40.91862 | 0.2843263 | 56.20025 | 0.7589582 | 0.0346228 | 0.1242706 |
| FRAHES | T1_1396 | sperry | 559 | -0.1590868 | -31.44527 | 0.2246880 | 44.41207 | 0.7743905 | 0.3955688 | 0.3079576 |
| FRAHES | T1_1396 | sureau | 559 | -0.2002308 | -39.57783 | 0.2586176 | 51.11863 | 0.7276322 | 0.1796239 | 0.2034540 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-332-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-332-2.png)

## Fontainebleau-Barbeau

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Fontainebleau-Barbeau |
| Country | France |
| SAPFLUXNET code | FRA_FON |
| SAPFLUXNET contributor (affiliation) | Nicolas Delpierre (Univ. Paris-Sud) |
| FLUXNET/ICOS code | FR-Fon |
| FLUXNET/ICOS contributor (affiliation) | Nicolas Delpierre (Univ. Paris-Sud) |
| Latitude (º) | 48.47634 |
| Longitude (º) | 2.78014 |
| Elevation (m) | 105 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material | Millstone |
| Soil texture | Loam |
| MAT (ºC) | 11.2 |
| MAP (mm) | 697 |
| Forest stand | Mixed deciduous forest |
| Stand LAI | 6 |
| Stand description DOI | 10.1111/nph.13771 |
| Species simulated | Quercus petraea, Carpinus betulus |
| Species parameter table | SpParamsFR |
| Simulation period | 2006-2008 |
| Evaluation period | 2006-2008 |

### Model inputs

#### Vegetation

| Species          | DBH | Height |     N | Z50 | Z95 |  LAI |
|:-----------------|----:|-------:|------:|:----|:----|-----:|
| Quercus petraea  |  33 |   2800 | 220.8 | NA  | NA  | 4.74 |
| Carpinus betulus |  10 |    500 | 883.2 | NA  | NA  | 1.26 |

#### Soil

| widths |     clay |     sand |        om |       bd |      rfc |
|-------:|---------:|---------:|----------:|---------:|---------:|
|    300 | 24.13333 | 26.56667 | 3.6733330 | 1.286667 | 14.43333 |
|    600 | 30.48571 | 25.22857 | 0.6228571 | 1.504286 | 15.54286 |
|   1100 | 30.20000 | 25.64000 | 0.4720000 | 1.518000 | 80.00000 |
|   2500 | 29.90000 | 25.80000 | 0.4500000 | 1.520000 | 90.00000 |

#### Custom traits

| Species          |    Al2As |
|:-----------------|---------:|
| Quercus petraea  | 1075.267 |
| Carpinus betulus | 1075.267 |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Taken from SoilGrids with theta_sat and theta_res modified |
| Vegetation | Plantation |
| Weather | From fluxnet data |
| Sapflow | Sapwood area estimated from dbh for trees within missing data |
| Eddy covariance | Variables H_F_MDS and LE_F_MDS for sensible and latent heat |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-343-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-344-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-345-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONTAINEBLEAU | sperry | 842 | -3.606711 | -272.5310 | 6.359181 | 480.5136 | 0.3470153 | -8.020759 | -1.991221 |
| FONTAINEBLEAU | sureau | 842 | -3.189159 | -240.9798 | 6.100736 | 460.9849 | 0.3540096 | -7.459807 | -1.869654 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-365-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-366-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONTAINEBLEAU | sperry | 842 | 0.4050401 | 8.350385 | 1.723501 | 35.53203 | 0.7197648 | 0.4641103 | 0.3767237 |
| FONTAINEBLEAU | sureau | 842 | 0.0763086 | 1.573193 | 1.930411 | 39.79773 | 0.6693312 | 0.3689898 | 0.3018981 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-372-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-373-1.png)

#### Gross primary productivity

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONTAINEBLEAU | granier | 842 | -1.600863 | -25.33391 | 3.223897 | 51.01869 | 0.7971266 | 0.3927816 | 0.2646893 |
| FONTAINEBLEAU | sperry | 842 | 3.762197 | 59.53738 | 3.924695 | 62.10894 | 0.8491455 | 0.0980251 | 0.1048504 |
| FONTAINEBLEAU | sureau | 842 | 3.398290 | 53.77850 | 3.524887 | 55.78191 | 0.8642702 | 0.2725630 | 0.1960391 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-379-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-380-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONTAINEBLEAU | T1_2856 | granier | 740 | -0.2345005 | -48.51650 | 0.2792133 | 57.76725 | 0.7503480 | -0.0574703 | 0.0125488 |
| FONTAINEBLEAU | T1_2856 | sperry | 740 | -0.2065900 | -42.74201 | 0.2542840 | 52.60957 | 0.6687282 | 0.0373690 | 0.1007122 |
| FONTAINEBLEAU | T1_2856 | sureau | 740 | -0.2259329 | -46.74392 | 0.2712834 | 56.12663 | 0.6330826 | -0.0920662 | 0.0405930 |
| FONTAINEBLEAU | T2_730 | granier | 732 | -0.0400328 | -25.03275 | 0.0559700 | 34.99831 | 0.8127820 | 0.5267715 | 0.4033449 |
| FONTAINEBLEAU | T2_730 | sperry | 732 | -0.0698322 | -43.66644 | 0.0871546 | 54.49823 | 0.7351887 | 0.0718052 | 0.0709082 |
| FONTAINEBLEAU | T2_730 | sureau | 732 | -0.0961035 | -60.09402 | 0.0989634 | 61.88232 | 0.7107221 | -0.2278494 | -0.0549766 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-387-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-387-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-387-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-387-4.png)

## Font-Blanche

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Font-Blanche |
| Country | France |
| SAPFLUXNET code |  |
| SAPFLUXNET contributor (affiliation) |  |
| FLUXNET/ICOS code | FR-Fbn |
| FLUXNET/ICOS contributor (affiliation) | Nicolas Martin-StPaul (INRAE) |
| Latitude (º) | 43.24 |
| Longitude (º) | 5.68 |
| Elevation (m) | 420 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material | Cretaceous limestone |
| Soil texture | Clay loam |
| MAT (ºC) | 13.5 |
| MAP (mm) | 722 |
| Forest stand | Mixed forest with P. halepensis and Q. ilex |
| Stand LAI | 2 |
| Stand description DOI | 10.1016/j.agrformet.2021.108472 |
| Species simulated | Quercus ilex, Pinus halepensis, Phillyrea latifolia |
| Species parameter table | SpParamsFR |
| Simulation period | 2014-2018 |
| Evaluation period | 2014-2018 |

### Model inputs

#### Vegetation

| Species             |       DBH |    Height |    N | Z50 |  Z95 |       LAI | Cover |
|:--------------------|----------:|----------:|-----:|----:|-----:|----------:|------:|
| Phillyrea latifolia |  2.587859 |  323.0000 | 1248 | 390 | 1470 | 0.0000000 |    NA |
| Pinus halepensis    | 26.759914 | 1195.7667 |  256 | 300 | 1200 | 0.9843761 |    NA |
| Quercus ilex        |  6.220031 |  495.5532 | 3104 | 500 | 2287 | 1.7156239 |    NA |
| Herbaceous layer    |        NA |   10.0000 |   NA |  NA |   NA |        NA |     5 |

#### Soil

| widths | clay | sand |  om |   bd | rfc |
|-------:|-----:|-----:|----:|-----:|----:|
|    300 |   39 |   26 |   6 | 1.45 |  50 |
|    700 |   39 |   26 |   3 | 1.45 |  65 |
|   1000 |   39 |   26 |   1 | 1.45 |  90 |
|   2500 |   39 |   26 |   1 | 1.45 |  95 |

#### Custom traits

| Species | VCstem_P12 | VCstem_P50 | VCstem_P88 | VCstem_slope | VCroot_P12 | VCroot_P50 | VCroot_P88 | VCleaf_kmax | LeafEPS | LeafPI0 | LeafAF | StemEPS | StemPI0 | StemAF | Gswmin | Gswmax | Gs_P50 | Gs_slope | Al2As |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Phillyrea latifolia | -1.971750 | -6.50 | -11.028250 | 11 | NA | NA | NA | 3.00 | 12.38 | -2.13 | 0.5 | 12.38 | -2.13 | 0.4 | 0.002 | 0.2200 | -2.207094 | 89.41176 | NA |
| Pinus halepensis | -3.707158 | -4.79 | -5.872842 | 46 | -1 | -1.741565 | -2.301482 | 4.00 | 5.31 | -1.50 | 0.6 | 5.00 | -1.65 | 0.4 | 0.001 | 0.2175 | -1.871216 | 97.43590 | 631.000 |
| Quercus ilex | -4.739642 | -6.40 | -8.060358 | 30 | NA | NA | NA | 2.63 | 15.00 | -2.50 | 0.4 | 15.00 | -2.50 | 0.4 | 0.002 | 0.2200 | -2.114188 | 44.70588 | 1540.671 |

#### Custom control

||
||
||

#### Remarks

| Title      | Remark                        |
|:-----------|:------------------------------|
| Soil       | Equal to Puechabon            |
| Vegetation |                               |
| Weather    | Missing values for some dates |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-398-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-399-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-400-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONBLA | sperry | 1004 | -1.0312666 | -21.03789 | 4.148988 | 84.63955 | 0.6794009 | -0.0111644 | 0.1032998 |
| FONBLA | sureau | 1004 | -0.5974221 | -12.18744 | 3.923233 | 80.03413 | 0.6714751 | 0.0901733 | 0.1520913 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-420-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-421-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONBLA | sperry | 1026 | -0.3652071 | -12.45583 | 1.470337 | 50.14761 | 0.4707739 | -0.0835329 | 0.0279308 |
| FONBLA | sureau | 1026 | -0.5430908 | -18.52277 | 1.572532 | 53.63310 | 0.5235614 | -0.0819579 | -0.0396324 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-427-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-428-1.png)

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONBLA | granier | 1760 | -0.1699596 | -38.13552 | 0.1700243 | 38.15005 | 0.8828316 | -0.7937427 | -0.3302394 |
| FONBLA | sperry | 1760 | -0.1751688 | -39.30438 | 0.1752191 | 39.31565 | 0.9086845 | -0.8189105 | -0.3708824 |
| FONBLA | sureau | 1760 | -0.1779463 | -39.92757 | 0.1779938 | 39.93825 | 0.8990984 | -0.8467524 | -0.3925914 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-435-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-436-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONBLA | T2_2631 | granier | 300 | 0.1401927 | 68.163335 | 0.1584582 | 77.04425 | 0.5902966 | -2.4193670 | -0.6013013 |
| FONBLA | T2_2631 | sperry | 300 | 0.0275698 | 13.404748 | 0.0872013 | 42.39833 | 0.6916213 | -0.0917146 | 0.1187855 |
| FONBLA | T2_2631 | sureau | 300 | 0.1513536 | 73.589901 | 0.2121985 | 103.17343 | 0.4361402 | -7.6028654 | -1.1443750 |
| FONBLA | T3_2854 | granier | 309 | -0.0093674 | -3.236196 | 0.0624460 | 21.57338 | 0.8856051 | 0.7812994 | 0.5745685 |
| FONBLA | T3_2854 | sperry | 309 | 0.1010601 | 34.913486 | 0.1267154 | 43.77669 | 0.9037622 | 0.1600268 | 0.1367146 |
| FONBLA | T3_2854 | sureau | 309 | 0.0036870 | 1.273770 | 0.0789663 | 27.28068 | 0.8932891 | 0.6079389 | 0.4620194 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-442-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-442-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-442-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-442-4.png)

#### Leaf water potential

| Site | Cohort | WP | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| FONBLA | T2_2631 | Midday | sperry | 3 | 0.8040548 | 30.139498 | 0.8040548 | 30.13950 | 0.9821158 | -12.5520966 | -3.2989070 |
| FONBLA | T2_2631 | Midday | sureau | 3 | 0.3309003 | 12.403593 | 0.5346218 | 20.03997 | 0.9850975 | -5.3541430 | -1.8583742 |
| FONBLA | T2_2631 | Predawn | sperry | 3 | 0.6687932 | 34.385254 | 0.6687932 | 34.38525 | 0.9789705 | -1.7476112 | -0.8809808 |
| FONBLA | T2_2631 | Predawn | sureau | 3 | -0.0566764 | -2.913956 | 0.4424417 | 22.74765 | 0.9902574 | -0.4115934 | -0.2443673 |
| FONBLA | T3_2854 | Midday | sperry | 3 | 0.6338100 | 23.225938 | 0.6338100 | 23.22594 | 0.9419772 | -2.0661325 | -0.9579944 |
| FONBLA | T3_2854 | Midday | sureau | 3 | 0.4927896 | 18.058251 | 0.6171812 | 22.61657 | 0.9879538 | -2.3052748 | -0.9066238 |
| FONBLA | T3_2854 | Predawn | sperry | 3 | 0.2788256 | 18.370645 | 0.2788256 | 18.37064 | 0.9917150 | 0.6655061 | 0.3700176 |
| FONBLA | T3_2854 | Predawn | sureau | 3 | -0.4311340 | -28.405610 | 0.4762455 | 31.37781 | 0.9977193 | -0.5152720 | -0.0760358 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-447-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-447-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-447-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-447-4.png)

## Collelongo

### General information

| Attribute                              | Value                              |
|:---------------------------------------|:-----------------------------------|
| Plot name                              | Collelongo                         |
| Country                                | Italy                              |
| SAPFLUXNET code                        |                                    |
| SAPFLUXNET contributor (affiliation)   |                                    |
| FLUXNET/ICOS code                      | IT-Col                             |
| FLUXNET/ICOS contributor (affiliation) | Giorgio Matteucci (IEIF CNR)       |
| Latitude (º)                           | 13.5881                            |
| Longitude (º)                          | 41.8494                            |
| Elevation (m)                          | 1560                               |
| Slope (º)                              | 19.29                              |
| Aspect (º)                             | 252                                |
| Parent material                        | Calcareous                         |
| Soil texture                           | Silt loam                          |
| MAT (ºC)                               | 6.3                                |
| MAP (mm)                               | 1180                               |
| Forest stand                           | European beech forest              |
| Stand LAI                              | 5.5                                |
| Stand description DOI                  | 10.1111/j.1365-2486.1996.tb00072.x |
| Species simulated                      | Fagus sylvatica                    |
| Species parameter table                | SpParamsFR                         |
| Simulation period                      | 2011-2013                          |
| Evaluation period                      | 2011-2013                          |

### Model inputs

#### Vegetation

| Species         |  DBH | Height |        N | Z50 | Z95 | LAI |  CR |
|:----------------|-----:|-------:|---------:|:----|:----|----:|----:|
| Fagus sylvatica | 20.2 | 1898.9 | 899.9668 | NA  | NA  | 5.5 | 0.5 |

#### Soil

| widths |     clay |  sand |    om |   bd |  rfc |
|-------:|---------:|------:|------:|-----:|-----:|
|    300 | 27.23333 | 32.50 | 3.000 | 1.37 | 17.4 |
|    700 | 30.90000 | 32.15 | 1.955 | 1.37 | 20.9 |
|   2000 | 31.20000 | 33.70 | 1.430 | 1.44 | 22.8 |

#### Custom traits

| Species         | Vmax298 | Jmax298 |
|:----------------|--------:|--------:|
| Fagus sylvatica |    94.5 |   159.9 |

#### Custom control

| freeDrainage |
|:-------------|
| FALSE        |

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Taken from SoilGrids |
| Vegetation | No understory |
| Weather | From V. saponaro |
| Sapflow | Not available |
| Eddy covariance | Variables H_F_MDS and LE_F_MDS for sensible and latent heat |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-457-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-458-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-459-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| COLLELONGO | sperry | 1096 | -1.2353106 | -41.81980 | 4.522891 | 153.1165 | 0.3852017 | -2.950934 | -0.8492370 |
| COLLELONGO | sureau | 1096 | -0.4686333 | -15.86496 | 4.505644 | 152.5326 | 0.4063043 | -2.870977 | -0.8421852 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-479-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-480-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| COLLELONGO | sperry | 1096 | 2.469258 | 91.75938 | 2.700000 | 100.33390 | 0.5709577 | -1.0376325 | -0.3203585 |
| COLLELONGO | sureau | 1096 | 1.556553 | 57.84263 | 2.236573 | 83.11262 | 0.4595437 | -0.6997558 | -0.0937326 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-486-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-487-1.png)

#### Gross primary productivity

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| COLLELONGO | granier | 1096 | 3.235170 | 80.99168 | 3.980176 | 99.64272 | 0.6740509 | -0.0570138 | 0.0527023 |
| COLLELONGO | sperry | 1096 | 4.145140 | 103.77253 | 4.656980 | 116.58633 | 0.6502941 | -0.4528247 | -0.1083797 |
| COLLELONGO | sureau | 1096 | 3.050072 | 76.35779 | 3.794291 | 94.98913 | 0.6881477 | 0.0345015 | 0.0969437 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-493-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-494-1.png)

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| COLLELONGO | granier | 1096 | -0.0658845 | -21.19950 | 0.0679496 | 21.86398 | 0.7151824 | -0.1954632 | -0.1093978 |
| COLLELONGO | sperry | 1096 | -0.0631716 | -20.32657 | 0.0664443 | 21.37961 | 0.7174353 | -0.1398915 | -0.0848206 |
| COLLELONGO | sureau | 1096 | -0.0504427 | -16.23083 | 0.0566589 | 18.23102 | 0.7105910 | 0.0893095 | 0.0749419 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-500-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-501-1.png)

## Mitra

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Mitra II |
| Country | Portugal |
| SAPFLUXNET code | PRT_MIT |
| SAPFLUXNET contributor (affiliation) | Teresa David (INIAV IP) |
| FLUXNET/ICOS code | PT-Mi1 |
| FLUXNET/ICOS contributor (affiliation) | Joao Santos Pereira |
| Latitude (º) | 38.54056 |
| Longitude (º) | -8.00028 |
| Elevation (m) | 235 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material | Granite |
| Soil texture | Sand |
| MAT (ºC) | 16.5 |
| MAP (mm) | 584 |
| Forest stand | Evergreen forest dominated by Quercus ilex subsp. rotundifolia |
| Stand LAI | 0.55 (trees) |
| Stand description DOI | 10.1093/treephys/27.6.793 |
| Species simulated | Quercus ilex |
| Species parameter table | SpParamsES |
| Simulation period | 2001-2003 |
| Evaluation period | 2001-2003 |

### Model inputs

#### Vegetation

| Species          |  DBH | Height |   N | Z50 |  Z95 |  LAI | Cover |
|:-----------------|-----:|-------:|----:|----:|-----:|-----:|------:|
| Quercus ilex     | 38.9 |    750 |  30 | 529 | 2287 | 0.55 |    NA |
| Herbaceous layer |   NA |     15 |  NA |  NA |   NA |   NA |   100 |

#### Soil

| widths |     clay | sand |   om |       bd |      rfc |
|-------:|---------:|-----:|-----:|---------:|---------:|
|    300 | 14.86667 | 58.4 | 2.92 | 1.463333 | 20.46667 |
|    700 | 15.05000 | 63.8 | 1.10 | 1.535000 | 21.65000 |
|   1000 | 14.70000 | 64.1 | 0.87 | 1.550000 | 80.00000 |
|   2000 | 14.70000 | 64.1 | 0.00 | 1.550000 | 90.00000 |

#### Custom traits

| Species | SLA | VCleaf_P12 | VCleaf_P50 | VCleaf_P88 | VCleaf_slope | VCstem_P12 | VCstem_P50 | VCstem_P88 | VCstem_slope | VCroot_P12 | VCroot_P50 | VCroot_P88 | VCroot_slope | VCleaf_kmax | Kmax_stemxylem | LeafEPS | LeafPI0 | LeafAF | StemEPS | StemPI0 | StemAF | Gswmin | Gswmax | Gs_P50 | Gs_slope | Al2As |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Quercus ilex | 4.55 | -4.004731 | -5.25 | -6.495269 | 40 | -4.739642 | -6.4 | -8.060358 | 30 | -2.004731 | -3.25 | -4.495269 | 40 | 2.63 | 0.2 | 15 | -2.5 | 0.4 | 15 | -2.5 | 0.4 | 0.002 | 0.2 | -2.114188 | 44.70588 | 1540.671 |

#### Custom control

||
||
||

#### Remarks

| Title                | Remark                                  |
|:---------------------|:----------------------------------------|
| Soil                 | Taken from SoilGrids                    |
| Vegetation           | No understory but 100% herbaceous cover |
| Soil moisture data   | Not available                           |
| Eddy covariance data | Not enough quality for evaluation       |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-514-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-515-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-516-1.png)

### Evaluation results

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| MITRA | T1_394 | granier | 1060 | -0.6229574 | -60.16683 | 0.6353503 | 61.36377 | 0.8356874 | -1.1755837 | -0.6869876 |
| MITRA | T1_394 | sperry | 1060 | -0.1281138 | -12.37356 | 0.4068453 | 39.29417 | 0.8428604 | -0.0564714 | -0.0802592 |
| MITRA | T1_394 | sureau | 1060 | -0.4413545 | -42.62715 | 0.4576423 | 44.20027 | 0.8848686 | -0.1708583 | -0.2151357 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-540-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-540-2.png)

## Rinconada

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Rinconada |
| Country | Spain |
| SAPFLUXNET code | ESP_RIN |
| SAPFLUXNET contributor (affiliation) | Virginia Hernandez-Santana (IRNAS-CSIC) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 40.600278 |
| Longitude (º) | -6.016667 |
| Elevation (m) | 1200 |
| Slope (º) | 10 |
| Aspect (º) | 0 |
| Parent material |  |
| Soil texture | Silty loam |
| MAT (ºC) | 10 |
| MAP (mm) | 1000 |
| Forest stand | Young, homogeneous, Quercus pyrenaica regrowth forest |
| Stand LAI | 3.4 |
| Stand description DOI | 10.1016/j.foreco.2008.03.004 |
| Species simulated | Quercus pyrenaica |
| Species parameter table | SpParamsES |
| Simulation period | 2006-2007 |
| Evaluation period | 2006-2007 |

### Model inputs

#### Vegetation

| Species           |  DBH | Height |    N | Z50 |  Z95 | LAI | Cover |
|:------------------|-----:|-------:|-----:|----:|-----:|----:|------:|
| Quercus pyrenaica | 11.7 |    740 | 1975 | 300 | 1500 | 3.4 |    NA |
| Herbaceous layer  |   NA |     10 |   NA |  NA |   NA |  NA |     5 |

#### Soil

| widths |  clay |     sand |    om |   bd |   rfc | VG_theta_sat | VG_theta_res |
|-------:|------:|---------:|------:|-----:|------:|-------------:|-------------:|
|    250 | 19.10 | 45.33333 | 4.000 | 1.48 | 10.00 |         0.35 |         0.03 |
|    250 | 23.95 | 41.60000 | 2.000 | 1.48 | 21.00 |         0.35 |         0.03 |
|    500 | 23.95 | 41.60000 | 1.315 | 1.48 | 24.85 |         0.35 |         0.03 |
|   1000 | 24.50 | 42.30000 | 0.820 | 1.51 | 60.00 |         0.35 |         0.03 |
|   2500 | 24.50 | 42.30000 | 0.000 | 1.56 | 85.00 |         0.35 |         0.03 |

#### Custom traits

| Species           | Kmax_stemxylem | VCleaf_kmax | Gswmin | Gswmax |    Al2As |
|:------------------|---------------:|------------:|-------:|-------:|---------:|
| Quercus pyrenaica |              1 |           4 |  0.003 |    0.3 | 4189.325 |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Taken from SoilGrids, with modification of theta_sat and theta_res |
| Vegetation | Understory not considered |
| Weather | Available weather complemented with interpolation |
| Sapflow | Sapflow scaling needs to be revised |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-551-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-552-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-553-1.png)

### Evaluation results

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| ESPRIN | granier | 28 | -0.0195095 | -11.23664 | 0.0480846 | 27.69469 | 0.7098820 | 0.3592953 | 0.2698148 |
| ESPRIN | sperry | 28 | -0.0281330 | -16.20344 | 0.0387631 | 22.32592 | 0.8707302 | 0.4977594 | 0.4113653 |
| ESPRIN | sureau | 28 | -0.0212111 | -12.21672 | 0.0459664 | 26.47471 | 0.6964435 | 0.3840233 | 0.3019803 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-576-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-577-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| ESPRIN | T1_402 | granier | 103 | 0.0433034 | 11.96871 | 0.0728124 | 20.12477 | 0.6405361 | -0.1256516 | -0.229276 |
| ESPRIN | T1_402 | sperry | 103 | 0.1983072 | 54.81056 | 0.2280273 | 63.02496 | 0.5927588 | -9.5199068 | -2.849737 |
| ESPRIN | T1_402 | sureau | 103 | 0.1329784 | 36.75419 | 0.2084606 | 57.61688 | 0.4911223 | -9.1444055 | -2.519397 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-583-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-583-2.png)

#### Leaf water potential

| Site | Cohort | WP | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| ESPRIN | T1_402 | Midday | sperry | 15 | 0.1044321 | 6.603999 | 0.3103223 | 19.62393 | 0.7753966 | 0.5280601 | 0.2771297 |
| ESPRIN | T1_402 | Midday | sureau | 15 | -0.5315740 | -33.615273 | 0.6189078 | 39.13802 | 0.7824509 | -1.8029680 | -0.4416945 |
| ESPRIN | T1_402 | Predawn | sperry | 15 | -0.1655488 | -41.345477 | 0.2690001 | 67.18224 | 0.6183757 | -2.2962836 | -0.6655821 |
| ESPRIN | T1_402 | Predawn | sureau | 15 | -0.6366121 | -158.992625 | 0.6834435 | 170.68866 | 0.6111314 | -24.4482717 | -3.2317129 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-588-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-588-2.png)

## Vallcebre (Barrol)

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Vallcebre (Cal Barrol) |
| Country | Spain |
| SAPFLUXNET code | ESP_VAL_BAR |
| SAPFLUXNET contributor (affiliation) | Rafael Poyatos (CREAF) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 42.202933 |
| Longitude (º) | 1.820486 |
| Elevation (m) | 1102 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material | Limestone |
| Soil texture | Silty clay loam |
| MAT (ºC) | 9.3 |
| MAP (mm) | 603 |
| Forest stand | Semi-deciduous sub-Mediterranean oak forest |
| Stand LAI | 2.1 |
| Stand description DOI | 10.1093/treephys/27.4.537 |
| Species simulated | Quercus pubescens |
| Species parameter table | SpParamsES |
| Simulation period | 2004-2005 |
| Evaluation period | 2004-2005 |

### Model inputs

#### Vegetation

| Species            |      DBH | Height |   N | Z50 | Z95 | LAI | Cover |
|:-------------------|---------:|-------:|----:|:----|:----|----:|------:|
| Quercus pubescens  | 21.82917 | 1162.5 | 828 | NA  | NA  | 2.1 |    NA |
| Buxus sempervirens |       NA |  100.0 |  NA | NA  | NA  |  NA |    20 |
| Herbaceous layer   |       NA |   20.0 |  NA | NA  | NA  |  NA |     5 |

#### Soil

| widths | clay | sand | om  |   bd | rfc |
|-------:|-----:|-----:|:----|-----:|----:|
|    100 | 31.0 |  9.7 | NA  | 1.23 |   5 |
|    100 | 31.0 |  9.7 | NA  | 1.30 |  10 |
|    100 | 31.0 |  9.7 | NA  | 1.30 |  15 |
|    200 | 28.6 |  8.8 | NA  | 1.50 |  20 |
|   1000 | 28.6 |  8.8 | NA  | 1.50 |  70 |
|   2000 | 28.6 |  8.8 | NA  | 1.50 |  90 |

#### Custom traits

| Species            |    Al2As |
|:-------------------|---------:|
| Quercus pubescens  | 2342.065 |
| Buxus sempervirens |       NA |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Soil depth is 50 cm but additional rocky layers were added |
| Vegetation | Understory modelled using B. sempervirens |
| Weather | Missing values have been complemented with interpolated data |
| Sapflow | Sapflux density has been scaled to cohort level using measured plant Huber values |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-598-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-599-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-600-1.png)

### Evaluation results

#### Soil water content (SWC.2)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| QVALLCEBRE | granier | 445 | 0.0063995 | 2.153539 | 0.0350804 | 11.80519 | 0.3762598 | 0.0938464 | 0.0276431 |
| QVALLCEBRE | sperry | 445 | 0.0038097 | 1.282038 | 0.0351415 | 11.82577 | 0.3700599 | 0.1043881 | 0.0259479 |
| QVALLCEBRE | sureau | 445 | 0.0056629 | 1.905677 | 0.0353343 | 11.89063 | 0.3545169 | 0.0762063 | 0.0206054 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-623-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-624-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| QVALLCEBRE | T1_400 | granier | 279 | -0.6776365 | -87.78636 | 0.6830502 | 88.48770 | 0.4276319 | -3.955211 | -1.586766 |
| QVALLCEBRE | T1_400 | sperry | 279 | -0.7719149 | -99.99992 | 0.7719149 | 99.99992 | 0.0113448 | -5.251250 | -1.923304 |
| QVALLCEBRE | T1_400 | sureau | 279 | -0.7719155 | -100.00000 | 0.7719155 | 100.00000 | NA | -5.251259 | -1.923306 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-630-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-630-2.png)

## Vallcebre (Sort)

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Vallcebre (Cal Sort) |
| Country | Spain |
| SAPFLUXNET code | ESP_VAL_SOR |
| SAPFLUXNET contributor (affiliation) | Rafael Poyatos (CREAF) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 42.196053 |
| Longitude (º) | 1.813561 |
| Elevation (m) | 1257 |
| Slope (º) | 10 |
| Aspect (º) | 0 |
| Parent material | Limestone |
| Soil texture | Sandy clay loam |
| MAT (ºC) | 8.5 |
| MAP (mm) | 623 |
| Forest stand | Pinus sylvestris forest in a terraced area |
| Stand LAI | 2.4 |
| Stand description DOI | 10.5194/hess-9-493-2005 |
| Species simulated | Pinus sylvestris |
| Species parameter table | SpParamsES |
| Simulation period | 2003-2005 |
| Evaluation period | 2003-2005 |

### Model inputs

#### Vegetation

| Species            |  DBH |   Height |    N | Z50 | Z95 | LAI | Cover |
|:-------------------|-----:|---------:|-----:|:----|:----|----:|------:|
| Pinus sylvestris   | 16.2 | 1076.923 | 2165 | NA  | NA  | 2.4 |    NA |
| Buxus sempervirens |   NA |  100.000 |   NA | NA  | NA  |  NA |     5 |
| Herbaceous layer   |   NA |   20.000 |   NA | NA  | NA  |  NA |     5 |

#### Soil

| widths | clay | sand |  om |   bd | rfc |
|-------:|-----:|-----:|----:|-----:|----:|
|    100 |   22 |   59 |   4 | 1.18 |  10 |
|    100 |   21 |   60 |   3 | 1.28 |  15 |
|    100 |   20 |   61 |   2 | 1.38 |  19 |
|    350 |   18 |   62 |   1 | 1.48 |  20 |
|    350 |   18 |   62 |   0 | 1.50 |  50 |

#### Custom traits

| Species            |    Al2As |
|:-------------------|---------:|
| Pinus sylvestris   | 681.6332 |
| Buxus sempervirens |       NA |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Soil depth is 65 cm (30 + 35) but an additional layer of 35 cm is considered with 50% rocks |
| Vegetation | Understory modelled using B. sempervirens |
| Weather | Missing values have been complemented with interpolated data |
| Sapflow | Sapflux density has been scaled to cohort level using measured plant Huber values |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-641-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-642-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-643-1.png)

### Evaluation results

#### Soil water content (SWC.2)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| PVALLCEBRE | granier | 743 | -0.0509432 | -23.73060 | 0.0536858 | 25.00814 | 0.6243016 | -0.5221338 | -0.1332389 |
| PVALLCEBRE | sperry | 743 | -0.0543417 | -25.31367 | 0.0559257 | 26.05154 | 0.7073752 | -0.5384204 | -0.1805204 |
| PVALLCEBRE | sureau | 743 | -0.0518721 | -24.16328 | 0.0539346 | 25.12405 | 0.6808283 | -0.4838932 | -0.1384911 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-666-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-667-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| PVALLCEBRE | T1_361 | granier | 733 | -0.4791828 | -64.67002 | 0.4933712 | 66.58487 | 0.7180062 | -0.6175715 | -0.2284064 |
| PVALLCEBRE | T1_361 | sperry | 733 | -0.3783439 | -51.06091 | 0.4034255 | 54.44589 | 0.7134381 | -0.1279842 | -0.0044575 |
| PVALLCEBRE | T1_361 | sureau | 733 | -0.4580281 | -61.81501 | 0.4807692 | 64.88412 | 0.6066207 | -0.5637218 | -0.1970296 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-673-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-673-2.png)

## Prades

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Prades (Tillar valley) |
| Country | Spain |
| SAPFLUXNET code | ESP_TIL_MIX |
| SAPFLUXNET contributor (affiliation) | Rafael Poyatos (CREAF) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 41.33263 |
| Longitude (º) | 1.014429 |
| Elevation (m) | 1018 |
| Slope (º) | 35 |
| Aspect (º) | 8.53 |
| Parent material | Fractured schist |
| Soil texture | Clay loam |
| MAT (ºC) | 10.5 |
| MAP (mm) | 651.274727491089 |
| Forest stand | Mixed forest with P. sylvestris (overstory) Q. ilex (midstory) |
| Stand LAI | 3.27 |
| Stand description DOI | 10.1111/nph.12278 |
| Species simulated | Quercus ilex, Pinus sylvestris |
| Species parameter table | SpParamsES |
| Simulation period | 2010-2013 |
| Evaluation period | 2010-2013 |

### Model inputs

#### Vegetation

| Species          |  DBH | Height |    N | Z50 |  Z95 |  LAI | Cover |
|:-----------------|-----:|-------:|-----:|----:|-----:|-----:|------:|
| Pinus sylvestris | 27.7 |   1424 |  257 | 300 | 1200 | 0.58 |    NA |
| Quercus ilex     |  8.4 |    500 | 2913 | 529 | 2287 | 2.69 |    NA |
| Herbaceous layer |   NA |     20 |   NA |  NA |   NA |   NA |    10 |

#### Soil

| widths | clay | sand |  om |  bd | rfc |
|-------:|-----:|-----:|----:|----:|----:|
|    300 |   21 |   47 |   4 | 1.5 |  45 |
|    700 |   19 |   48 |   4 | 1.5 |  70 |
|   1000 |   19 |   48 |   4 | 1.5 |  85 |
|   2500 |   19 |   48 |   4 | 1.5 |  90 |

#### Custom traits

| Species | VCleaf_P12 | VCleaf_P50 | VCleaf_P88 | VCleaf_slope | VCstem_P12 | VCstem_P50 | VCstem_P88 | VCstem_slope | VCroot_P12 | VCroot_P50 | VCroot_P88 | VCroot_slope | VCleaf_kmax | LeafEPS | LeafPI0 | LeafAF | StemEPS | StemPI0 | StemAF | Gswmin | Gswmax | Gs_P50 | Gs_slope | Al2As |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pinus sylvestris | NA | NA | NA | NA | -0.9930548 | -3.2 | -5.406945 | 22.57 | -0.2474341 | -1.65 | -3.052566 | 35.51402 | 4.00 | 5.31 | -1.5 | 0.6 | 5 | -1.65 | 0.4 | 0.001 | 0.18 | -1.871216 | 97.43590 | 594.5372 |
| Quercus ilex | -4.004731 | -5.25 | -6.495269 | 40 | -4.7396415 | -6.4 | -8.060358 | 30.00 | -2.0047311 | -3.25 | -4.495269 | 40.00000 | 2.63 | 15.00 | -2.5 | 0.4 | 15 | -2.50 | 0.4 | 0.002 | 0.20 | -2.114188 | 44.70588 | 1387.0312 |

#### Custom control

||
||
||

#### Remarks

| Title      | Remark                                      |
|:-----------|:--------------------------------------------|
| Soil       | Additional rocky layer considered           |
| Vegetation | Understory not considered                   |
| Sapflow    | No scaling required (already per leaf area) |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-684-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-685-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-686-1.png)

### Evaluation results

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| PRADES | granier | 1219 | 0.0126859 | 6.658040 | 0.0277414 | 14.55971 | 0.9037732 | 0.7603946 | 0.5356481 |
| PRADES | sperry | 1219 | -0.0341963 | -17.947483 | 0.0392720 | 20.61141 | 0.8562335 | 0.4985917 | 0.3426418 |
| PRADES | sureau | 1219 | 0.0007093 | 0.372253 | 0.0225653 | 11.84312 | 0.9122336 | 0.8319628 | 0.6222882 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-709-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-710-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| PRADES | T1_361 | granier | 1157 | 0.0240317 | 9.151837 | 0.1383128 | 52.67284 | 0.7185727 | 0.3949991 | 0.2946511 |
| PRADES | T1_361 | sperry | 1157 | -0.1675409 | -63.803595 | 0.1765273 | 67.22586 | 0.6477598 | -0.0182773 | 0.0997696 |
| PRADES | T1_361 | sureau | 1157 | -0.1551231 | -59.074616 | 0.1917968 | 73.04084 | 0.6909798 | 0.0464280 | 0.0219005 |
| PRADES | T2_394 | granier | 908 | 0.0645376 | 39.003423 | 0.0866384 | 52.36007 | 0.8293234 | 0.0794978 | 0.1450859 |
| PRADES | T2_394 | sperry | 908 | 0.2009877 | 121.467344 | 0.3005758 | 181.65364 | 0.4345681 | -14.5269977 | -1.9659672 |
| PRADES | T2_394 | sureau | 908 | 0.1496196 | 90.422936 | 0.1759120 | 106.31279 | 0.8200349 | -3.6734347 | -0.7358323 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-716-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-716-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-716-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-716-4.png)

#### Leaf water potential

| Site | Cohort | WP | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| PRADES | T1_361 | Midday | sperry | 13 | -2.9939678 | -174.690735 | 3.1125045 | 181.60706 | 0.7892417 | -137.0450274 | -12.3921068 |
| PRADES | T1_361 | Midday | sureau | 13 | -0.5858098 | -34.180575 | 1.0191117 | 59.46269 | 0.8358349 | -17.9219085 | -3.3849103 |
| PRADES | T1_361 | Predawn | sperry | 13 | -2.9303878 | -244.970124 | 3.1070671 | 259.73989 | 0.8887888 | -68.7034130 | -7.7761247 |
| PRADES | T1_361 | Predawn | sureau | 13 | -0.7148896 | -59.762255 | 1.1098711 | 92.78133 | 0.9493691 | -10.6162747 | -2.1349073 |
| PRADES | T2_394 | Midday | sperry | 9 | -1.9387935 | -67.343343 | 2.1691204 | 75.34367 | 0.5430751 | -4.6989131 | -1.5499884 |
| PRADES | T2_394 | Midday | sureau | 9 | 0.6122760 | 21.267202 | 0.8465349 | 29.40411 | 0.9600492 | 0.1342736 | 0.0048252 |
| PRADES | T2_394 | Predawn | sperry | 9 | -1.9412741 | -91.893768 | 2.1773976 | 103.07111 | 0.6976748 | -3.6843692 | -1.3949843 |
| PRADES | T2_394 | Predawn | sureau | 9 | 0.0807108 | 3.820595 | 0.7749422 | 36.68331 | 0.9693239 | 0.5026515 | 0.1476181 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-721-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-721-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-721-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-721-4.png)

## Can Balasc

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Can Balasc |
| Country | Spain |
| SAPFLUXNET code | ESP_CAN |
| SAPFLUXNET contributor (affiliation) | Elisenda Sánchez-Costa (IDAEA-CSIC) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 41.43099 |
| Longitude (º) | 2.0736 |
| Elevation (m) | 270 |
| Slope (º) | 0.86 |
| Aspect (º) | 90 |
| Parent material | Shales and granite |
| Soil texture | Sandy loam |
| MAT (ºC) | 17 |
| MAP (mm) | 585 |
| Forest stand | Mixed forest dominated by Q. ilex |
| Stand LAI | 3.2 |
| Stand description DOI | 10.1016/j.agrformet.2015.03.012 |
| Species simulated | Quercus ilex, Quercus pubescens, Pinus halepensis, Arbutus unedo |
| Species parameter table | SpParamsES |
| Simulation period | 2011-2012 |
| Evaluation period | 2011-2012 |

### Model inputs

#### Vegetation

| Species                |  DBH |  Height |    N | Z50 |  Z95 |       LAI | Cover |
|:-----------------------|-----:|--------:|-----:|----:|-----:|----------:|------:|
| Arbutus unedo          |  9.6 |  810.00 |   76 | 390 | 1470 | 0.1402566 |    NA |
| Pinus halepensis       | 33.7 | 1710.00 |   53 | 300 | 1200 | 0.3854073 |    NA |
| Quercus pubescens      | 12.0 |  960.00 |  150 | 529 | 2287 | 0.3641975 |    NA |
| Quercus ilex           | 11.9 | 1020.00 | 1150 | 529 | 2287 | 2.3101386 |    NA |
| Arbutus unedo          |   NA |  174.00 |   NA |  NA |   NA |        NA |  4.83 |
| Phillyrea angustifolia |   NA |  153.33 |   NA |  NA |   NA |        NA |  7.25 |
| Pistacia lentiscus     |   NA |  118.33 |   NA |  NA |   NA |        NA | 13.50 |
| Quercus ilex           |   NA |   78.00 |   NA |  NA |   NA |        NA |  9.67 |
| Viburnum spp.          |   NA |  138.33 |   NA |  NA |   NA |        NA |  9.67 |
| Herbaceous layer       |   NA |   20.00 |   NA |  NA |   NA |        NA |  5.00 |

#### Soil

| widths |  clay | sand |   om |  bd | rfc |
|-------:|------:|-----:|-----:|----:|----:|
|    300 | 20.23 | 48.9 | 2.70 | 1.5 |  20 |
|    700 | 24.58 | 52.4 | 1.00 | 1.5 |  30 |
|   1000 | 27.66 | 45.6 | 0.61 | 1.5 |  85 |
|   2500 | 27.66 | 45.6 | 0.61 | 1.5 |  90 |

#### Custom traits

| Species | VCstem_P12 | VCstem_P50 | VCstem_P88 | VCstem_slope | VCroot_P12 | VCroot_P50 | VCroot_P88 | VCleaf_kmax | LeafEPS | LeafPI0 | LeafAF | StemEPS | StemPI0 | StemAF | Gswmin | Gswmax | Gs_P50 | Gs_slope | Al2As |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Arbutus unedo | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 1297.3893 |
| Pinus halepensis | -3.707158 | -4.79 | -5.872842 | 46 | -1 | -1.741565 | -2.301482 | 4.00 | 5.31 | -1.5 | 0.6 | 5 | -1.65 | 0.4 | 0.001 | 0.2175 | -1.871216 | 97.43590 | 631.1784 |
| Quercus pubescens | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | 1487.5004 |
| Quercus ilex | -4.739642 | -6.40 | -8.060358 | 30 | NA | NA | NA | 2.63 | 15.00 | -2.5 | 0.4 | 15 | -2.50 | 0.4 | 0.002 | 0.2200 | -2.114188 | 44.70588 | 1009.0124 |
| Phillyrea angustifolia | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| Pistacia lentiscus | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| Viburnum spp. | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Soil description from local samples |
| Vegetation | Understory composed of multiple species |
| Sapflow | No scaling required (already per leaf area) but tree selection may be required |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-731-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-732-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-733-1.png)

### Evaluation results

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| CANBALASC | granier | 676 | 0.0473977 | 30.16485 | 0.0486203 | 30.94294 | 0.7837562 | -0.5697617 | -0.4046595 |
| CANBALASC | sperry | 676 | 0.0227535 | 14.48081 | 0.0374101 | 23.80855 | 0.7848474 | -0.0671201 | -0.0807927 |
| CANBALASC | sureau | 676 | 0.0386241 | 24.58114 | 0.0450748 | 28.68653 | 0.7509600 | -0.4852941 | -0.3022297 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-756-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-757-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| CANBALASC | T1_35 | granier | 644 | 0.0494477 | 16.643602 | 0.1357419 | 45.68939 | 0.6400983 | -0.1034996 | 0.0395725 |
| CANBALASC | T1_35 | sperry | 644 | 0.0279878 | 9.420433 | 0.1986003 | 66.84690 | 0.8221526 | -1.3167017 | -0.4051753 |
| CANBALASC | T1_35 | sureau | 644 | -0.0290858 | -9.790005 | 0.1849189 | 62.24191 | 0.4793934 | -0.7324706 | -0.3083746 |
| CANBALASC | T2_356 | granier | 611 | 0.2042640 | 78.764532 | 0.2715383 | 104.70559 | 0.2072645 | -4.1388446 | -1.0356453 |
| CANBALASC | T2_356 | sperry | 611 | 0.0312130 | 12.035781 | 0.1953536 | 75.32867 | 0.6019596 | -1.4808097 | -0.4645106 |
| CANBALASC | T2_356 | sureau | 611 | 0.1343729 | 51.814398 | 0.3096059 | 119.38452 | 0.1766265 | -6.0043058 | -1.3210273 |
| CANBALASC | T3_400 | granier | 679 | -0.3222220 | -81.260809 | 0.3264160 | 82.31849 | 0.0855845 | -2.2387665 | -0.8925127 |
| CANBALASC | T3_400 | sperry | 679 | -0.3964901 | -99.990398 | 0.3964901 | 99.99040 | 0.0670838 | -3.3697095 | -1.2987921 |
| CANBALASC | T3_400 | sureau | 679 | -0.3791250 | -95.611103 | 0.3791250 | 95.61110 | 0.1764338 | -3.0672721 | -1.1981115 |
| CANBALASC | T4_394 | granier | 644 | -0.0473004 | -16.098599 | 0.1184825 | 40.32531 | 0.7451627 | 0.5040159 | 0.3535929 |
| CANBALASC | T4_394 | sperry | 644 | 0.0961102 | 32.710923 | 0.2020600 | 68.77073 | 0.7741900 | -0.3487563 | -0.1023820 |
| CANBALASC | T4_394 | sureau | 644 | 0.0538849 | 18.339641 | 0.1685543 | 57.36714 | 0.6858141 | -0.1182911 | 0.0804155 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-763-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-763-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-763-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-763-4.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-763-5.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-763-6.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-763-7.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-763-8.png)

## Alto-Tajo Armallones

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Alto Tajo (Armallones) |
| Country | Spain |
| SAPFLUXNET code | ESP_ALT_ARM |
| SAPFLUXNET contributor (affiliation) | Alicia Forner (MNCN-CSIC) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 40.7769 |
| Longitude (º) | -2.3283 |
| Elevation (m) | 1079 |
| Slope (º) | 25.64 |
| Aspect (º) | 270 |
| Parent material | Cretaceous and Jurassic limestone |
| Soil texture | Clay |
| MAT (ºC) | 10.1 |
| MAP (mm) | 495 |
| Forest stand | Sparse mixed forest dominated by three species |
| Stand LAI | 10.1007/S11258-014-0351-x |
| Stand description DOI | 1.09 |
| Species simulated | Pinus nigra, Quercus faginea, Quercus ilex |
| Species parameter table | SpParamsES |
| Simulation period | 2012-2013 |
| Evaluation period | 2012-2013 |

### Model inputs

#### Vegetation

| Species          |     DBH |    Height |         N | Z50 |  Z95 |       LAI | Cover |
|:-----------------|--------:|----------:|----------:|----:|-----:|----------:|------:|
| Pinus nigra      | 25.4720 | 1208.1085 |  94.36517 | 300 | 1200 | 0.2318096 |    NA |
| Quercus faginea  | 16.7450 |  752.6195 | 240.75280 | 529 | 2287 | 0.5673799 |    NA |
| Quercus ilex     | 22.3075 |  902.0182 |  90.43751 | 529 | 2287 | 0.2908105 |    NA |
| Herbaceous layer |      NA |   20.0000 |        NA |  NA |   NA |        NA |    10 |

#### Soil

| widths |    clay | sand |   om |       bd | rfc |
|-------:|--------:|-----:|-----:|---------:|----:|
|    300 | 21.8667 | 41.4 | 4.26 | 1.243333 |  45 |
|    700 | 23.8000 | 42.2 | 0.87 | 1.510000 |  65 |
|   1000 | 24.1000 | 41.7 | 0.55 | 1.560000 |  90 |
|   2500 | 24.1000 | 41.7 | 0.55 | 1.560000 |  95 |

#### Custom traits

| Species | VCleaf_P12 | VCleaf_P50 | VCleaf_P88 | VCleaf_slope | VCstem_P12 | VCstem_P50 | VCstem_P88 | VCstem_slope | VCroot_P12 | VCroot_P50 | VCroot_P88 | VCroot_slope | VCleaf_kmax | LeafEPS | LeafPI0 | LeafAF | StemEPS | StemPI0 | StemAF | Gswmin | Gswmax | Gs_P50 | Gs_slope | Al2As |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pinus nigra | NA | NA | NA | NA | -0.9930548 | -3.2 | -5.406945 | 22.57 | -0.2474341 | -1.65 | -3.052566 | 35.51402 | NA | 5.31 | -1.5 | 0.6 | 5 | -1.65 | 0.4 | 0.001 | 0.18 | -1.871216 | 97.43590 | 1272 |
| Quercus faginea | -4.004731 | -5.25 | -6.495269 | 40 | -4.7396415 | -6.4 | -8.060358 | 30.00 | -2.0047311 | -3.25 | -4.495269 | 40.00000 | NA | 15.00 | -2.5 | 0.4 | 15 | -2.50 | 0.4 | 0.002 | 0.22 | -2.114188 | 44.70588 | 1488 |
| Quercus ilex | -4.004731 | -5.25 | -6.495269 | 40 | -4.7396415 | -6.4 | -8.060358 | 30.00 | -2.0047311 | -3.25 | -4.495269 | 40.00000 | 2.63 | 15.00 | -2.5 | 0.4 | 15 | -2.50 | 0.4 | 0.002 | 0.22 | -2.114188 | 44.70588 | 1541 |

#### Custom control

||
||
||

#### Remarks

| Title      | Remark                                  |
|:-----------|:----------------------------------------|
| Soil       | Taken from Soilgrids                    |
| Vegetation | Understory not considered               |
| Sapflow    | Scaling done using tree density and LAI |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-774-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-775-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-776-1.png)

### Evaluation results

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| ESPALTARM | granier | 715 | 0.1433899 | 259.1233 | 0.1433899 | 259.1233 | 0.7870933 | -14.63805 | -3.623193 |
| ESPALTARM | sperry | 715 | 0.1247126 | 225.3710 | 0.1247126 | 225.3710 | 0.8700708 | -11.15875 | -3.020994 |
| ESPALTARM | sureau | 715 | 0.1454450 | 262.8371 | 0.1454450 | 262.8371 | 0.7472993 | -15.22292 | -3.689452 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-799-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-800-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| ESPALTARM | T1_357 | granier | 627 | -1.1582563 | -77.28377 | 1.1918684 | 79.52651 | 0.5283444 | -0.5309958 | -0.0059475 |
| ESPALTARM | T1_357 | sperry | 627 | -1.2524345 | -83.56774 | 1.2525842 | 83.57772 | 0.8292409 | -0.5009595 | -0.0571921 |
| ESPALTARM | T1_357 | sureau | 627 | -1.4392528 | -96.03304 | 1.4424187 | 96.24429 | 0.0356331 | -1.0305610 | -0.2174142 |
| ESPALTARM | T2_392 | granier | 537 | -1.2373645 | -71.96913 | 1.2804475 | 74.47498 | 0.8003660 | -0.6686507 | -0.1476499 |
| ESPALTARM | T2_392 | sperry | 537 | -0.6744082 | -39.22577 | 0.9357470 | 54.42608 | 0.6625246 | 0.1133615 | 0.1613011 |
| ESPALTARM | T2_392 | sureau | 537 | -1.3775515 | -80.12286 | 1.4091252 | 81.95929 | 0.3292837 | -1.1445880 | -0.2629822 |
| ESPALTARM | T3_394 | granier | 628 | -1.2903450 | -78.91057 | 1.2918604 | 79.00325 | 0.8114469 | -1.0250536 | -0.3376533 |
| ESPALTARM | T3_394 | sperry | 628 | -0.7492384 | -45.81940 | 0.7655518 | 46.81704 | 0.8850627 | 0.3447780 | 0.2073115 |
| ESPALTARM | T3_394 | sureau | 628 | -1.3730589 | -83.96891 | 1.4195981 | 86.81500 | 0.3015672 | -1.4120144 | -0.4699189 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-806-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-806-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-806-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-806-4.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-806-5.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-806-6.png)

#### Leaf water potential

| Site | Cohort | WP | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| ESPALTARM | T1_357 | Midday | sperry | 4 | -1.900146 | -107.58083 | 1.900146 | 107.58083 | 0.9662285 | -67.0809648 | -5.0804685 |
| ESPALTARM | T1_357 | Midday | sureau | 4 | -13.553574 | -767.36444 | 13.571665 | 768.38870 | 0.5914709 | -2027.8222032 | -42.4293293 |
| ESPALTARM | T1_357 | Predawn | sperry | 4 | -1.118683 | -99.88239 | 1.274139 | 113.76239 | 0.9917442 | -9.5546907 | -1.9459855 |
| ESPALTARM | T1_357 | Predawn | sureau | 4 | -13.956630 | -1246.12772 | 14.171109 | 1265.27762 | 0.3265365 | -1003.1520667 | -31.7655708 |
| ESPALTARM | T2_392 | Midday | sperry | 4 | -1.059349 | -41.55327 | 1.059349 | 41.55327 | 0.8768157 | -16.2007805 | -3.5198872 |
| ESPALTARM | T2_392 | Midday | sureau | 4 | -11.299509 | -443.22664 | 12.169158 | 477.33887 | -0.4405119 | -3493.6472894 | -50.9217400 |
| ESPALTARM | T2_392 | Predawn | sperry | 4 | -1.278712 | -172.94498 | 1.299780 | 175.79447 | 0.8664038 | -9.8864451 | -2.0560597 |
| ESPALTARM | T2_392 | Predawn | sureau | 4 | -12.815870 | -1733.33833 | 12.997189 | 1757.86158 | -0.6728001 | -1238.8587688 | -29.5591513 |
| ESPALTARM | T3_394 | Midday | sperry | 4 | -1.427386 | -60.56266 | 1.427386 | 60.56266 | 0.9929924 | -3.8180473 | -1.8177888 |
| ESPALTARM | T3_394 | Midday | sureau | 4 | -5.560706 | -235.93555 | 6.117674 | 259.56719 | -0.3687785 | -128.1877828 | -11.0768398 |
| ESPALTARM | T3_394 | Predawn | sperry | 4 | -0.882941 | -77.40853 | 1.003353 | 87.96516 | 0.8401473 | -0.7370877 | -0.2259367 |
| ESPALTARM | T3_394 | Predawn | sureau | 4 | -5.530471 | -484.86318 | 6.074577 | 532.56567 | -0.5194238 | -73.6648544 | -6.4221638 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-811-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-811-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-811-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-811-4.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-811-5.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-811-6.png)

## Ronda

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Ronda (Pilones) |
| Country | Spain |
| SAPFLUXNET code | ESP_RON_PIL |
| SAPFLUXNET contributor (affiliation) | Víctor Lechuga (UJaén) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 36.692881 |
| Longitude (º) | -5.019568 |
| Elevation (m) | 1734 |
| Slope (º) | 15 |
| Aspect (º) | 315 |
| Parent material |  |
| Soil texture | Silty loam |
| MAT (ºC) | 8.1 |
| MAP (mm) | 925 |
| Forest stand | Mixed gimnosperm forest dominated by Abies pinsapo |
| Stand LAI | 10.3390/f10121132 |
| Stand description DOI | NA |
| Species simulated | Abies pinsapo, Taxus baccata |
| Species parameter table | SpParamsES |
| Simulation period | 2011-2013 |
| Evaluation period | 2011-2013 |

### Model inputs

#### Vegetation

| Species          |  DBH | Height |   N | Z50 | Z95 | Cover |
|:-----------------|-----:|-------:|----:|:----|:----|------:|
| Abies pinsapo    | 15.3 |   1256 | 486 | NA  | NA  |    NA |
| Taxus baccata    | 15.0 |    630 |  15 | NA  | NA  |    NA |
| Herbaceous layer |   NA |     10 |  NA | NA  | NA  |    10 |

#### Soil

| widths |     clay |     sand |   om |       bd |      rfc | VG_theta_sat | VG_theta_res |
|-------:|---------:|---------:|-----:|---------:|---------:|-------------:|-------------:|
|    300 | 18.96667 | 44.46667 | 3.97 | 1.276667 | 19.46667 |         0.55 |          0.1 |
|    700 | 19.65000 | 45.20000 | 1.18 | 1.420000 | 40.00000 |         0.55 |          0.1 |
|   1000 | 20.10000 | 45.70000 | 0.65 | 1.480000 | 80.00000 |         0.55 |          0.1 |
|   2500 | 20.10000 | 45.70000 | 0.00 | 1.480000 | 90.00000 |         0.55 |          0.1 |

#### Custom traits

| Species       |    Al2As | LeafAngle |
|:--------------|---------:|----------:|
| Abies pinsapo | 2587.510 |        NA |
| Taxus baccata | 6790.546 |        30 |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Taken from SoilGrids with theta_sat and theta_res modified |
| Vegetation | LAI not available. Understory not considered except herbaceous layer |
| Weather | Complemented with interpolated weather |
| Sapflow | Species-level Huber value used for scaling. Revise scaling. |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-821-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-822-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-823-1.png)

### Evaluation results

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| RONDA | granier | 731 | -0.0259815 | -9.803531 | 0.0470927 | 17.76938 | 0.7778536 | 0.4064516 | 0.3414886 |
| RONDA | sperry | 731 | -0.0315778 | -11.915174 | 0.0452281 | 17.06582 | 0.8168955 | 0.4185891 | 0.3675616 |
| RONDA | sureau | 731 | -0.0352805 | -13.312309 | 0.0445116 | 16.79545 | 0.8403874 | 0.4313748 | 0.3775813 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-846-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-847-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| RONDA | T1_2 | granier | 730 | -0.2194396 | -37.7685724 | 0.2581159 | 44.42529 | 0.7506175 | 0.0836503 | 0.0572780 |
| RONDA | T1_2 | sperry | 730 | -0.0193042 | -3.3225120 | 0.1691480 | 29.11270 | 0.7882288 | 0.5419261 | 0.3822173 |
| RONDA | T1_2 | sureau | 730 | -0.0029819 | -0.5132313 | 0.2350946 | 40.46303 | 0.7579859 | 0.1313242 | 0.1413589 |
| RONDA | T2_490 | granier | 712 | 0.4074759 | 136.2486973 | 0.4189946 | 140.10022 | 0.7591982 | -3.5353514 | -1.2325057 |
| RONDA | T2_490 | sperry | 712 | 0.3249809 | 108.6646355 | 0.3437419 | 114.93779 | 0.6871920 | -2.1953304 | -0.8315407 |
| RONDA | T2_490 | sureau | 712 | 0.4124187 | 137.9014375 | 0.4466308 | 149.34102 | 0.6614605 | -5.6726249 | -1.3797583 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-853-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-853-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-853-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-853-4.png)

## Davos Seehornwald

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Davos Seehornwald |
| Country | Switzerland |
| SAPFLUXNET code | CHE_DAV_SEE |
| SAPFLUXNET contributor (affiliation) | Roman Zweifel (WSL) |
| FLUXNET/ICOS code | CH-Dav |
| FLUXNET/ICOS contributor (affiliation) | Nina Buchmann (ETH) |
| Latitude (º) | 46.81668 |
| Longitude (º) | 9.856198 |
| Elevation (m) | 1650 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material |  |
| Soil texture | Loamy sand |
| MAT (ºC) | 3.8 |
| MAP (mm) | 840 |
| Forest stand | Subalpine coniferous (spruce) forest |
| Stand LAI | 3.9 |
| Stand description DOI | 10.1007/s10021-011-9481-3 |
| Species simulated | Picea abies |
| Species parameter table | SpParamsFR |
| Simulation period | 2009-2011 |
| Evaluation period | 2009-2011 |

### Model inputs

#### Vegetation

| Species     | DBH | Height |   N | Z50 | Z95 | LAI |
|:------------|----:|-------:|----:|:----|:----|----:|
| Picea abies |  20 |   2800 | 830 | NA  | NA  | 3.9 |

#### Soil

| widths |     clay |     sand |       om |       bd |      rfc |
|-------:|---------:|---------:|---------:|---------:|---------:|
|    300 | 13.53333 | 50.96667 | 9.176667 | 1.016667 | 12.73333 |
|    300 | 14.20000 | 51.80000 | 3.140000 | 1.300000 | 22.30000 |
|    400 | 15.30000 | 50.00000 | 3.530000 | 1.390000 | 23.90000 |
|   1000 | 15.20000 | 52.10000 | 3.550000 | 1.410000 | 80.00000 |
|   2500 | 15.20000 | 52.10000 | 3.550000 | 1.410000 | 90.00000 |

#### Custom traits

| Species     | Al2As |
|:------------|------:|
| Picea abies |  1975 |

#### Custom control

||
||
||

#### Remarks

| Title      | Remark                                        |
|:-----------|:----------------------------------------------|
| Soil       | Taken from SoilGrids                          |
| Vegetation | No understory or secondary species considered |
| Weather    | From fluxnet data                             |
| Sapflow    | Species-level Huber value used for scaling    |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-864-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-865-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-866-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| DAVOS | sperry | 1065 | -0.7657610 | -39.40939 | 5.585792 | 287.4692 | -0.0032606 | -3.275767 | -1.0616317 |
| DAVOS | sureau | 1065 | 0.7908215 | 40.69912 | 5.215662 | 268.4207 | 0.0741696 | -3.031016 | -0.9250223 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-886-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-887-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| DAVOS | sperry | 1065 | -0.6836765 | -12.31346 | 2.336597 | 42.08365 | 0.6566201 | 0.3910246 | 0.3287846 |
| DAVOS | sureau | 1065 | -2.6964103 | -48.56412 | 3.488684 | 62.83350 | 0.4947136 | -0.2591815 | -0.0021661 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-893-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-894-1.png)

#### Gross primary productivity

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| DAVOS | granier | 1065 | -0.4496603 | -12.69573 | 1.515589 | 42.79120 | 0.7523455 | 0.5225207 | 0.3683756 |
| DAVOS | sperry | 1065 | 3.5287276 | 99.63024 | 3.590429 | 101.37231 | 0.8295993 | -1.1537128 | -0.4963177 |
| DAVOS | sureau | 1065 | 2.5218941 | 71.20326 | 2.627119 | 74.17419 | 0.8419396 | -0.1377968 | -0.0948567 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-900-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-901-1.png)

#### Soil water content (SWC.1)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| DAVOS | granier | 365 | -0.0998496 | -42.45747 | 0.1000014 | 42.52202 | -0.1428519 | -7.452678 | -2.378283 |
| DAVOS | sperry | 365 | -0.0994479 | -42.28666 | 0.0999469 | 42.49886 | 0.3099641 | -7.226880 | -2.376444 |
| DAVOS | sureau | 365 | -0.0720454 | -30.63471 | 0.0907165 | 38.57395 | -0.5589951 | -5.807497 | -2.064618 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-907-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-908-1.png)

#### Soil water content (SWC.3)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| DAVOS | granier | 335 | -0.0721579 | -25.19281 | 0.0721579 | 25.19281 | -0.2659240 | -24.17779 | -4.715613 |
| DAVOS | sperry | 335 | -0.1463298 | -51.08881 | 0.1463298 | 51.08881 | 0.4148438 | -97.77566 | -10.590763 |
| DAVOS | sureau | 335 | -0.0710973 | -24.82254 | 0.0711698 | 24.84784 | 0.5853325 | -22.09837 | -4.637349 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-913-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-914-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| DAVOS | T1_2602 | granier | 364 | 0.0396852 | 31.21337 | 0.0637066 | 50.10682 | 0.8624273 | 0.6652376 | 0.4820753 |
| DAVOS | T1_2602 | sperry | 364 | 0.1450038 | 114.04912 | 0.1776036 | 139.68971 | 0.7222512 | -1.5776436 | -0.4438905 |
| DAVOS | T1_2602 | sureau | 364 | -0.0147712 | -11.61792 | 0.0722213 | 56.80391 | 0.7555197 | 0.5530534 | 0.4128514 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-920-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-920-2.png)

## Lötschental

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Lötschental |
| Country | Switzerland |
| SAPFLUXNET code | CHE_LOT_NOR |
| SAPFLUXNET contributor (affiliation) | Patrick Fonti (WSL) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 46.3918 |
| Longitude (º) | 7.7613 |
| Elevation (m) | 1300 |
| Slope (º) | 36.87 |
| Aspect (º) | 0 |
| Parent material | Calcareous |
| Soil texture | Loam |
| MAT (ºC) | 5 |
| MAP (mm) | 716 |
| Forest stand | Mixed evergreen Norway spruce and deciduous European larch forest |
| Stand LAI | 3 |
| Stand description DOI | 10.1016/j.agrformet.2012.08.002 |
| Species simulated | Picea abies, Larix decidua subsp. decidua |
| Species parameter table | SpParamsFR |
| Simulation period | 2013-2016 |
| Evaluation period | 2013-2016 |

### Model inputs

#### Vegetation

| Species                      |      DBH |   Height |        N | Z50 | Z95 |      LAI |
|:-----------------------------|---------:|---------:|---------:|:----|:----|---------:|
| Larix decidua subsp. decidua | 30.76667 | 1910.000 | 345.8887 | NA  | NA  | 1.417577 |
| Picea abies                  | 38.60000 | 1783.333 | 385.5873 | NA  | NA  | 1.580276 |

#### Soil

| widths | clay |     sand |       om |       bd | rfc | VG_theta_sat | VG_theta_res |
|-------:|-----:|---------:|---------:|---------:|----:|-------------:|-------------:|
|    100 | 23.2 | 42.13333 | 9.266667 | 1.083333 |  40 |         0.35 |         0.07 |
|    100 | 24.0 | 42.52000 | 5.220000 | 1.178000 |  40 |         0.40 |         0.10 |
|    200 | 23.3 | 43.70000 | 3.290000 | 1.273333 |  40 |         0.60 |         0.20 |
|    200 | 21.9 | 45.10000 | 2.550000 | 1.380000 |  45 |         0.32 |         0.10 |
|    400 | 22.2 | 47.40000 | 2.570000 | 1.450000 |  80 |         0.32 |         0.10 |
|   1000 | 22.8 | 46.10000 | 2.350000 | 1.490000 |  90 |         0.32 |         0.10 |
|   2000 | 22.8 | 46.10000 | 0.000000 | 1.490000 |  97 |         0.32 |         0.10 |

#### Custom traits

| Species                      | Kmax_stemxylem | Al2As | LeafAngle | LeafAngleSD |
|:-----------------------------|---------------:|------:|----------:|------------:|
| Larix decidua subsp. decidua |             20 |  3305 |        30 |          21 |
| Picea abies                  |             NA |  1975 |        42 |           2 |

#### Custom control

||
||
||

#### Remarks

| Title | Remark |
|:---|:---|
| Soil | Taken from SoilGrids with theta_sat and theta_res modified |
| Vegetation | Understory not considered. LAI not available |
| Weather | Complemented with weather station in the same valley for Radiation, Precipitation and WindSpeed (2013-2016); and all variables for 2016 |
| Soil moisture | Provided by A. Cabon |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-931-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-932-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-933-1.png)

### Evaluation results

#### Soil water content (SWC.1)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| LOTSCHENTAL | granier | 1262 | -0.0947313 | -57.23762 | 0.0947703 | 57.26123 | 0.3830117 | -2.522469 | -0.8178025 |
| LOTSCHENTAL | sperry | 1262 | -0.0944304 | -57.05585 | 0.0944748 | 57.08265 | 0.0203554 | -2.510329 | -0.8121335 |
| LOTSCHENTAL | sureau | 1262 | -0.0947084 | -57.22380 | 0.0947485 | 57.24804 | 0.2418909 | -2.522678 | -0.8173840 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-956-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-957-1.png)

#### Soil water content (SWC.3)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| LOTSCHENTAL | granier | 1253 | -0.0772424 | -18.05830 | 0.1169106 | 27.33224 | 0.4913111 | -0.2708428 | -0.1171372 |
| LOTSCHENTAL | sperry | 1253 | -0.1048706 | -24.51743 | 0.1490001 | 34.83435 | -0.0942243 | -1.1772032 | -0.4237675 |
| LOTSCHENTAL | sureau | 1253 | -0.0918774 | -21.47978 | 0.1250625 | 29.23804 | 0.2878240 | -0.5303322 | -0.1950321 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-962-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-963-1.png)

#### Soil water content (SWC.4)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| LOTSCHENTAL | granier | 1262 | -0.0302744 | -16.24065 | 0.0563511 | 30.22949 | 0.4193756 | -0.1488663 | -0.0504372 |
| LOTSCHENTAL | sperry | 1262 | -0.0406996 | -21.83327 | 0.0641166 | 34.39531 | 0.0775572 | -0.4873555 | -0.1951941 |
| LOTSCHENTAL | sureau | 1262 | -0.0340843 | -18.28449 | 0.0573594 | 30.77041 | 0.3154501 | -0.2296359 | -0.0692332 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-968-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-969-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| LOTSCHENTAL | T1_1956 | granier | 572 | -0.6139748 | -88.204672 | 0.6212300 | 89.24697 | 0.3345069 | -1.6696735 | -0.5442750 |
| LOTSCHENTAL | T1_1956 | sperry | 572 | -0.6656370 | -95.626543 | 0.6665345 | 95.75548 | 0.3016867 | -1.9589186 | -0.6568943 |
| LOTSCHENTAL | T1_1956 | sureau | 572 | -0.4274111 | -61.402610 | 0.4941607 | 70.99197 | 0.3804645 | -0.7123096 | -0.2284016 |
| LOTSCHENTAL | T2_2602 | granier | 695 | -0.1072958 | -39.502620 | 0.1677483 | 61.75913 | 0.5015690 | -0.0042430 | 0.0918991 |
| LOTSCHENTAL | T2_2602 | sperry | 695 | -0.1818994 | -66.969078 | 0.2787827 | 102.63818 | 0.3451237 | -1.9540071 | -0.5091831 |
| LOTSCHENTAL | T2_2602 | sureau | 695 | 0.0110793 | 4.079001 | 0.1465618 | 53.95899 | 0.5906479 | 0.1556494 | 0.2065916 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-975-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-975-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-975-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-975-4.png)

## Morgan-Monroe

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Morgan-Mornoe |
| Country | USA |
| SAPFLUXNET code | USA_MOR_SF |
| SAPFLUXNET contributor (affiliation) | Koong Yi (Indiana University Bloomington) |
| FLUXNET/ICOS code | US-MMS |
| FLUXNET/ICOS contributor (affiliation) | Kim Novick (Indiana University) |
| Latitude (º) | 39.323011 |
| Longitude (º) | -86.413321 |
| Elevation (m) | 275 |
| Slope (º) | 0 |
| Aspect (º) | 0 |
| Parent material |  |
| Soil texture | Silty clay loam |
| MAT (ºC) | 12 |
| MAP (mm) | 1159 |
| Forest stand | Mixed temperate forest |
| Stand LAI | 5 |
| Stand description DOI |  |
| Species simulated | Acer saccharum, Liriodendron tulipifera, Quercus rubra, Quercus alba |
| Species parameter table | SpParamsUS |
| Simulation period | 2011-2013 |
| Evaluation period | 2011-2013 |

### Model inputs

#### Vegetation

| Species                 |   DBH | Height | N   | Z50 | Z95 |  LAI | Cover |
|:------------------------|------:|-------:|:----|:----|:----|-----:|------:|
| Acer saccharum          | 42.50 |   2700 | NA  | NA  | NA  | 1.25 |    NA |
| Liriodendron tulipifera | 62.75 |   2700 | NA  | NA  | NA  | 1.25 |    NA |
| Quercus rubra           | 44.60 |   2700 | NA  | NA  | NA  | 1.25 |    NA |
| Quercus alba            | 37.50 |   2700 | NA  | NA  | NA  | 1.25 |    NA |
| Herbaceous layer        |    NA |     20 | NA  | NA  | NA  |   NA |    10 |

#### Soil

| widths |  clay |      sand |    om |    bd |  rfc |
|-------:|------:|----------:|------:|------:|-----:|
|    300 | 24.70 |  8.166667 | 1.940 | 1.430 |  3.0 |
|    700 | 27.95 | 11.000000 | 0.265 | 1.645 |  4.3 |
|   1000 | 35.30 | 20.100000 | 0.180 | 1.740 | 75.0 |
|   2500 | 35.30 | 20.100000 | 0.000 | 1.740 | 95.0 |

#### Custom traits

| Species                 |    Al2As |
|:------------------------|---------:|
| Acer saccharum          | 7365.000 |
| Liriodendron tulipifera | 6658.333 |
| Quercus rubra           | 5730.000 |
| Quercus alba            | 6880.000 |

#### Custom control

||
||
||

#### Remarks

| Title         | Remark                                   |
|:--------------|:-----------------------------------------|
| Soil          | Taken from SoilGrids                     |
| Vegetation    | Understory not considered                |
| Weather       | Using FLUXNET weather                    |
| Soil moisture | Taken from FLUXNET data                  |
| Sapflow       | Scaling and tree selection to be revised |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-986-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-987-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-988-1.png)

### Evaluation results

#### Sensible heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| USAMORSF | sperry | 1096 | -9.381381 | -285.2483 | 12.69905 | 386.1247 | 0.1719677 | -21.70736 | -3.997474 |
| USAMORSF | sureau | 1096 | -9.001336 | -273.6927 | 12.21847 | 371.5121 | 0.1605327 | -20.16373 | -3.808348 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1008-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1009-1.png)

#### Latent heat turbulent flux

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| USAMORSF | sperry | 1096 | -0.3354250 | -8.327449 | 3.007687 | 74.6705 | 0.3986226 | 0.1162605 | 0.1930587 |
| USAMORSF | sureau | 1096 | -0.3669866 | -9.111014 | 3.041175 | 75.5019 | 0.3951234 | 0.1129860 | 0.1840740 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1015-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1016-1.png)

#### Gross primary productivity

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| USAMORSF | granier | 1096 | -0.6868668 | -15.895021 | 3.636825 | 84.16102 | 0.6299199 | 0.2229482 | 0.1634666 |
| USAMORSF | sperry | 1096 | -0.1856127 | -4.295328 | 2.385368 | 55.20062 | 0.8625133 | 0.6393227 | 0.4513237 |
| USAMORSF | sureau | 1096 | 1.0964125 | 25.372458 | 2.305396 | 53.34996 | 0.8945023 | 0.6813222 | 0.4697186 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1022-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1023-1.png)

#### Soil water content (SWC)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| USAMORSF | granier | 1096 | 0.0093041 | 2.845939 | 0.0573540 | 17.54353 | 0.6525159 | 0.3437470 | 0.1678458 |
| USAMORSF | sperry | 1096 | 0.0161682 | 4.945548 | 0.0561446 | 17.17358 | 0.6818479 | 0.3733632 | 0.1853937 |
| USAMORSF | sureau | 1096 | 0.0195814 | 5.989600 | 0.0579825 | 17.73578 | 0.6462042 | 0.3242294 | 0.1587265 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1029-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1030-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| USAMORSF | T1_57 | granier | 159 | 0.1752095 | 99.54112 | 0.1752095 | 99.54112 | 0.4680211 | -9.452734 | -3.1649754 |
| USAMORSF | T1_57 | sperry | 159 | 0.3883679 | 220.64195 | 0.4076379 | 231.58978 | 0.4614672 | -61.608630 | -8.6901236 |
| USAMORSF | T1_57 | sureau | 159 | 0.1736351 | 98.64667 | 0.2004242 | 113.86625 | 0.3758731 | -18.086265 | -3.7643641 |
| USAMORSF | T2_1644 | granier | 154 | 0.1152704 | 51.21445 | 0.1421627 | 63.16261 | 0.1670790 | -1.340378 | -0.7192184 |
| USAMORSF | T2_1644 | sperry | 154 | 0.1493060 | 66.33641 | 0.1972770 | 87.64983 | 0.2169393 | -3.434234 | -1.3857340 |
| USAMORSF | T2_1644 | sureau | 154 | 0.1529975 | 67.97652 | 0.2261393 | 100.47327 | -0.0361855 | -6.134140 | -1.7347742 |
| USAMORSF | T3_2591 | granier | 119 | -0.0843033 | -64.13740 | 0.0887489 | 67.51956 | 0.1772907 | -6.044433 | -3.0697432 |
| USAMORSF | T3_2591 | sperry | 119 | -0.1310793 | -99.72425 | 0.1310793 | 99.72425 | 0.1095720 | -14.317584 | -5.0108818 |
| USAMORSF | T3_2591 | sureau | 119 | 0.0297873 | 22.66198 | 0.1080429 | 82.19830 | 0.1119496 | -13.251304 | -3.9545047 |
| USAMORSF | T4_2536 | granier | 119 | 0.2669789 | 269.42267 | 0.2669789 | 269.42267 | 0.1556662 | -115.127852 | -15.7938254 |
| USAMORSF | T4_2536 | sperry | 119 | 0.3527855 | 356.01468 | 0.3599122 | 363.20658 | 0.1519395 | -223.277532 | -21.6396237 |
| USAMORSF | T4_2536 | sureau | 119 | 0.3452678 | 348.42819 | 0.3487056 | 351.89739 | 0.1324050 | -243.505257 | -20.9346921 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1036-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1036-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1036-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1036-4.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1036-5.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1036-6.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1036-7.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1036-8.png)

## Sevilleta

### General information

| Attribute | Value |
|:---|:---|
| Plot name | Sevilleta |
| Country | USA |
| SAPFLUXNET code | USA_PJS_P04_AMB |
| SAPFLUXNET contributor (affiliation) | William Pockman (University of New Mexico, USA) |
| FLUXNET/ICOS code |  |
| FLUXNET/ICOS contributor (affiliation) |  |
| Latitude (º) | 34.386389 |
| Longitude (º) | -106.529444 |
| Elevation (m) | 1911 |
| Slope (º) | 1 |
| Aspect (º) | 0 |
| Parent material |  |
| Soil texture | Sandy loam |
| MAT (ºC) | 12.7 |
| MAP (mm) | 311 |
| Forest stand | Mixed pine-juniper forest |
| Stand LAI | 0.71 |
| Stand description DOI | 10.1890/ES11-00369.1 |
| Species simulated | Pinus edulis, Juniperus monosperma |
| Species parameter table | SpParamsUS |
| Simulation period | 2010-2016 |
| Evaluation period | 2010-2016 |

### Model inputs

#### Vegetation

| Species              |   DBH | Height |   N | Z50 |  Z95 |       LAI | Cover |
|:---------------------|------:|-------:|----:|----:|-----:|----------:|------:|
| Pinus edulis         | 19.18 |  484.2 |  77 | 150 | 2500 | 0.0844801 |    NA |
| Juniperus monosperma | 32.66 |  379.2 | 273 | 150 | 2500 | 0.6244180 |    NA |
| Herbaceous layer     |    NA |   20.0 |  NA |  NA |   NA |        NA |     5 |

#### Soil

| widths |     clay |     sand |        om |    bd |  rfc | VG_theta_sat | VG_theta_res |
|-------:|---------:|---------:|----------:|------:|-----:|-------------:|-------------:|
|    100 | 18.66667 | 51.96667 | 1.3033333 | 1.450 | 18.6 |          0.2 |        0.001 |
|    100 | 22.68000 | 48.82000 | 1.1380000 | 1.502 | 40.0 |          0.2 |        0.001 |
|    200 | 24.50000 | 47.43333 | 0.8200000 | 1.530 | 50.0 |          0.2 |        0.001 |
|    300 | 23.18182 | 48.30909 | 0.4518182 | 1.550 | 55.0 |          0.2 |        0.001 |
|    300 | 21.80000 | 49.40000 | 0.2700000 | 1.550 | 65.0 |          0.2 |        0.001 |
|   1000 | 21.00000 | 51.90000 | 0.1600000 | 1.540 | 90.0 |          0.2 |        0.001 |
|   2500 | 21.00000 | 51.90000 | 0.0000000 | 1.540 | 95.0 |          0.2 |        0.001 |

#### Custom traits

| Species              | Kmax_stemxylem |    Al2As |
|:---------------------|---------------:|---------:|
| Pinus edulis         |           0.25 | 1872.270 |
| Juniperus monosperma |           0.25 | 2577.026 |

#### Custom control

||
||
||

#### Remarks

| Title      | Remark                                                     |
|:-----------|:-----------------------------------------------------------|
| Soil       | Taken from SoilGrids with theta_sat and theta_res modified |
| Vegetation | Understory not considered                                  |

### Macroclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1047-1.png)

### Microclimate

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1048-1.png)

### Runoff & deep drainage

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1049-1.png)

### Evaluation results

#### Soil water content (SWC.2)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| SEVILLETA | granier | 2178 | 0.0002639 | 0.4230835 | 0.0201238 | 32.26474 | 0.6037450 | 0.3538163 | 0.2405065 |
| SEVILLETA | sperry | 2178 | -0.0075679 | -12.1336950 | 0.0185520 | 29.74462 | 0.6512673 | 0.3713968 | 0.2998286 |
| SEVILLETA | sureau | 2178 | -0.0024045 | -3.8551443 | 0.0202728 | 32.50352 | 0.5851166 | 0.3368382 | 0.2348858 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1072-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1073-1.png)

#### Soil water content (SWC.4)

| Site | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| SEVILLETA | granier | 2178 | 0.0240140 | 39.37384 | 0.0297670 | 48.80667 | 0.0764376 | -2.627320 | -1.189548 |
| SEVILLETA | sperry | 2178 | 0.0074503 | 12.21571 | 0.0282489 | 46.31760 | 0.0942935 | -2.262767 | -1.077884 |
| SEVILLETA | sureau | 2178 | 0.0206725 | 33.89514 | 0.0317713 | 52.09288 | 0.0417967 | -2.979154 | -1.336973 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1078-1.png)

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1079-1.png)

#### Transpiration per leaf area

| Site | Cohort | Mode | n | Bias | Bias.rel | MAE | MAE.rel | r | NSE | NSE.abs |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| SEVILLETA | T1_2222 | granier | 1787 | 0.4344575 | 275.94520 | 0.4349358 | 276.2490 | 0.3198403 | -16.950405 | -4.026980 |
| SEVILLETA | T1_2222 | sperry | 1787 | -0.0489987 | -31.12148 | 0.2049202 | 130.1548 | -0.0212776 | -4.039247 | -1.368463 |
| SEVILLETA | T1_2222 | sureau | 1787 | -0.0312583 | -19.85366 | 0.2250501 | 142.9403 | -0.0130848 | -7.124662 | -1.601125 |
| SEVILLETA | T2_1548 | granier | 1776 | 0.2266323 | 172.41553 | 0.2274706 | 173.0533 | 0.4968503 | -7.370381 | -2.126675 |
| SEVILLETA | T2_1548 | sperry | 1776 | 0.5906476 | 449.34821 | 0.5908405 | 449.4949 | 0.6100420 | -56.594584 | -7.121339 |
| SEVILLETA | T2_1548 | sureau | 1776 | 0.2818470 | 214.42130 | 0.2912637 | 221.5853 | 0.3536526 | -22.181278 | -3.003536 |

![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1085-1.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1085-2.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1085-3.png)![](StandLevelEvaluation_files/figure-html/unnamed-chunk-1085-4.png)
