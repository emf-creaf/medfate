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
| spwb     | Granier           | buckets     | total              | 0.027 |
| spwb     | Granier           | single      | total              | 0.048 |
| spwb     | Granier           | dual        | total              | 2.706 |
| spwb     | Sperry            | buckets     | total              | 1.095 |
| spwb     | Sperry            | single      | total              | 1.509 |
| spwb     | Sperry            | dual        | total              | 4.306 |
| spwb     | Sureau            | buckets     | total              | 0.395 |
| spwb     | Sureau            | single      | total              | 0.759 |
| spwb     | Sureau            | dual        | total              | 6.742 |
| growth   | Granier           | buckets     | total              | 0.026 |
| growth   | Granier           | single      | total              | 0.036 |
| growth   | Granier           | dual        | total              | 2.553 |
| growth   | Sperry            | buckets     | total              | 2.186 |
| growth   | Sperry            | single      | total              | 2.487 |
| growth   | Sperry            | dual        | total              | 5.609 |
| growth   | Sureau            | buckets     | total              | 1.517 |
| growth   | Sureau            | single      | total              | 1.882 |
| growth   | Sureau            | dual        | total              | 7.925 |
| fordyn   | Granier           | buckets     | total              | 0.067 |
| fordyn   | Granier           | single      | total              | 0.072 |
| fordyn   | Granier           | dual        | total              | 2.688 |
| fordyn   | Sperry            | buckets     | total              | 2.381 |
| fordyn   | Sperry            | single      | total              | 2.610 |
| fordyn   | Sperry            | dual        | total              | 5.857 |
| fordyn   | Sureau            | buckets     | total              | 1.718 |
| fordyn   | Sureau            | single      | total              | 2.038 |
| fordyn   | Sureau            | dual        | total              | 8.143 |
