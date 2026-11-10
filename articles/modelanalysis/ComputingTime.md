# Computing time estimates

## About this vignette

The aim of this vignette is to provide users with a rough estimation of
computing times for simulation models included in package **medfate**.

The results presented here were obtained using **1 year of simulation**
with the **example data sets** on a laptop (32 GiB memory and Intel®
Core™ Ultra 9 185H × 22) with Ubuntu Linux OS 26.04.1 LTS.

## Table of computational times

Computational times were estimated using
[`system.time()`](https://rdrr.io/r/base/system.time.html), are in
**seconds** and are shown by package version.

| function | transpirationMode | soilDomains | rhizosphereOverlap | 5.2.0 |
|:---------|:------------------|:------------|:-------------------|------:|
| spwb     | Granier           | buckets     | total              | 0.054 |
| spwb     | Granier           | single      | total              | 0.113 |
| spwb     | Granier           | dual        | total              | 0.125 |
| spwb     | Sperry            | buckets     | total              | 3.109 |
| spwb     | Sperry            | single      | total              | 3.096 |
| spwb     | Sperry            | dual        | total              | 3.084 |
| spwb     | Sureau            | buckets     | total              | 1.638 |
| spwb     | Sureau            | single      | total              | 1.679 |
| spwb     | Sureau            | dual        | total              | 1.694 |
| growth   | Granier           | buckets     | total              | 0.060 |
| growth   | Granier           | single      | total              | 0.102 |
| growth   | Granier           | dual        | total              | 0.124 |
| growth   | Sperry            | buckets     | total              | 5.167 |
| growth   | Sperry            | single      | total              | 5.057 |
| growth   | Sperry            | dual        | total              | 5.027 |
| growth   | Sureau            | buckets     | total              | 3.491 |
| growth   | Sureau            | single      | total              | 3.505 |
| growth   | Sureau            | dual        | total              | 3.534 |
| fordyn   | Granier           | buckets     | total              | 0.116 |
| fordyn   | Granier           | single      | total              | 0.145 |
| fordyn   | Granier           | dual        | total              | 0.162 |
| fordyn   | Sperry            | buckets     | total              | 5.075 |
| fordyn   | Sperry            | single      | total              | 5.102 |
| fordyn   | Sperry            | dual        | total              | 5.089 |
| fordyn   | Sureau            | buckets     | total              | 3.507 |
| fordyn   | Sureau            | single      | total              | 3.559 |
| fordyn   | Sureau            | dual        | total              | 3.584 |
