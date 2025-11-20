# Computing time estimates

## About this vignette

The aim of this vignette is to provide users with a rough estimation of
computing times for simulation models included in package **medfate**.

The results presented here were obtained using **1 year of simulation**
with the **example data sets** on a laptop (16 GiB memory and 11th Gen
Inter Core I5 processor @ 2.40 GHz x 8) with Ubuntu Linux OS.

## Table of computational times

Computational times were estimated using
[`system.time()`](https://rdrr.io/r/base/system.time.html), are in
**seconds** and are shown by package version.

| function | transpirationMode | soilDomains |  4.4.0 |  4.7.0 |  4.8.0 |  4.8.1 |  4.8.2 |  4.8.3 |  4.8.4 |  4.8.5 |
|:---------|:------------------|:------------|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|
| spwb     | Granier           | buckets     |  0.382 |  0.117 |  0.136 |  0.117 |  0.105 |  0.097 |  0.120 |  0.149 |
| spwb     | Granier           | single      |  0.426 |  0.143 |  0.181 |  0.155 |  0.148 |  0.137 |  0.168 |  0.245 |
| spwb     | Granier           | dual        |  6.850 |  5.090 |  5.954 |  5.509 |  6.181 |  5.158 |  5.250 | 10.119 |
| spwb     | Sperry            | buckets     |  9.531 |  6.681 |  7.139 |  7.256 |  6.921 |  6.620 |  7.268 | 15.059 |
| spwb     | Sperry            | single      | 10.208 |  6.726 |  6.689 |  7.210 |  6.820 |  6.789 |  7.621 | 15.583 |
| spwb     | Sperry            | dual        | 18.375 | 12.623 | 12.412 | 13.090 | 12.724 | 12.790 | 15.282 | 27.505 |
| spwb     | Sureau            | buckets     |  9.730 |  6.679 |  5.806 |  6.077 |  5.327 |  2.680 |  2.905 |  6.078 |
| spwb     | Sureau            | single      | 10.565 |  6.555 |  5.947 |  6.481 |  5.408 |  2.718 |  3.018 |  6.582 |
| spwb     | Sureau            | dual        | 19.013 | 15.222 | 13.566 | 14.561 | 13.509 | 10.639 | 11.573 | 23.596 |
| growth   | Granier           | buckets     |  1.576 |  0.209 |  0.179 |  0.177 |  0.169 |  0.174 |  0.256 |  0.308 |
| growth   | Granier           | single      |  1.762 |  0.254 |  0.219 |  0.246 |  0.225 |  0.236 |  0.314 |  0.464 |
| growth   | Granier           | dual        |  8.468 |  5.805 |  5.305 |  5.500 |  5.687 |  5.735 |  6.289 | 12.407 |
| growth   | Sperry            | buckets     | 17.586 | 15.088 | 14.140 | 15.051 | 15.084 | 14.942 | 17.610 | 25.363 |
| growth   | Sperry            | single      | 16.947 | 15.025 | 14.281 | 15.679 | 15.318 | 14.657 | 20.996 | 25.133 |
| growth   | Sperry            | dual        | 25.108 | 21.341 | 20.170 | 22.366 | 22.867 | 21.458 | 25.284 | 34.520 |
| growth   | Sureau            | buckets     | 17.009 | 14.063 | 13.986 | 14.535 | 13.653 | 10.783 | 11.825 | 16.444 |
| growth   | Sureau            | single      | 18.023 | 14.200 | 14.426 | 13.553 | 12.887 | 11.180 | 12.183 | 16.510 |
| growth   | Sureau            | dual        | 26.702 | 22.414 | 24.041 | 21.116 | 20.494 | 20.467 | 22.507 | 29.463 |
| fordyn   | Granier           | buckets     |  1.545 |  0.245 |  0.247 |  0.220 |  0.207 |  0.227 |  0.320 |  0.304 |
| fordyn   | Granier           | single      |  1.783 |  0.292 |  0.296 |  0.270 |  0.246 |  0.271 |  0.370 |  0.400 |
| fordyn   | Granier           | dual        |  8.168 |  6.134 |  5.821 |  5.656 |  6.158 |  5.728 |  6.456 |  9.432 |
| fordyn   | Sperry            | buckets     | 18.099 | 16.756 | 15.161 | 14.597 | 14.693 | 16.044 | 16.247 | 25.371 |
| fordyn   | Sperry            | single      | 18.595 | 16.675 | 15.030 | 14.640 | 14.524 | 16.274 | 16.349 | 27.511 |
| fordyn   | Sperry            | dual        | 27.134 | 25.285 | 21.272 | 21.205 | 20.392 | 23.035 | 23.592 | 39.191 |
| fordyn   | Sureau            | buckets     | 18.832 | 17.116 | 14.224 | 13.731 | 13.312 | 11.135 | 10.675 | 18.308 |
| fordyn   | Sureau            | single      | 18.604 | 16.090 | 14.312 | 14.000 | 13.712 | 11.076 | 11.050 | 18.947 |
| fordyn   | Sureau            | dual        | 26.987 | 28.045 | 22.149 | 21.737 | 22.640 | 19.623 | 18.702 | 33.148 |
