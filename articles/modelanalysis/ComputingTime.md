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

| function | transpirationMode | soilDomains | 4.4.0 | 4.7.0 | 4.8.0 | 4.8.1 | 4.8.3 | 4.9.0 | 5.0.0 |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| spwb | Granier | buckets | 0.382 | 0.117 | 0.136 | 0.117 | 0.097 | 0.134 | 0.018 |
| spwb | Granier | single | 0.426 | 0.143 | 0.181 | 0.155 | 0.137 | 0.223 | 0.036 |
| spwb | Granier | dual | 6.850 | 5.090 | 5.954 | 5.509 | 5.158 | 8.914 | 4.263 |
| spwb | Sperry | buckets | 9.531 | 6.681 | 7.139 | 7.256 | 6.620 | 13.478 | 2.102 |
| spwb | Sperry | single | 10.208 | 6.726 | 6.689 | 7.210 | 6.789 | 13.529 | 2.122 |
| spwb | Sperry | dual | 18.375 | 12.623 | 12.412 | 13.090 | 12.790 | 23.679 | 6.765 |
| spwb | Sureau | buckets | 9.730 | 6.679 | 5.806 | 6.077 | 2.680 | 5.857 | 0.607 |
| spwb | Sureau | single | 10.565 | 6.555 | 5.947 | 6.481 | 2.718 | 6.052 | 0.717 |
| spwb | Sureau | dual | 19.013 | 15.222 | 13.566 | 14.561 | 10.639 | 18.987 | 6.720 |
| growth | Granier | buckets | 1.576 | 0.209 | 0.179 | 0.177 | 0.174 | 0.257 | 0.037 |
| growth | Granier | single | 1.762 | 0.254 | 0.219 | 0.246 | 0.236 | 0.363 | 0.054 |
| growth | Granier | dual | 8.468 | 5.805 | 5.305 | 5.500 | 5.735 | 9.360 | 4.155 |
| growth | Sperry | buckets | 17.586 | 15.088 | 14.140 | 15.051 | 14.942 | 24.665 | 9.802 |
| growth | Sperry | single | 16.947 | 15.025 | 14.281 | 15.679 | 14.657 | 24.559 | 9.847 |
| growth | Sperry | dual | 25.108 | 21.341 | 20.170 | 22.366 | 21.458 | 35.205 | 14.492 |
| growth | Sureau | buckets | 17.009 | 14.063 | 13.986 | 14.535 | 10.783 | 16.685 | 8.020 |
| growth | Sureau | single | 18.023 | 14.200 | 14.426 | 13.553 | 11.180 | 17.075 | 8.498 |
| growth | Sureau | dual | 26.702 | 22.414 | 24.041 | 21.116 | 20.467 | 36.420 | 14.379 |
| fordyn | Granier | buckets | 1.545 | 0.245 | 0.247 | 0.220 | 0.227 | 0.420 | 0.090 |
| fordyn | Granier | single | 1.783 | 0.292 | 0.296 | 0.270 | 0.271 | 0.545 | 0.098 |
| fordyn | Granier | dual | 8.168 | 6.134 | 5.821 | 5.656 | 5.728 | 11.312 | 4.155 |
| fordyn | Sperry | buckets | 18.099 | 16.756 | 15.161 | 14.597 | 16.044 | 24.664 | 9.790 |
| fordyn | Sperry | single | 18.595 | 16.675 | 15.030 | 14.640 | 16.274 | 24.484 | 9.781 |
| fordyn | Sperry | dual | 27.134 | 25.285 | 21.272 | 21.205 | 23.035 | 35.603 | 14.392 |
| fordyn | Sureau | buckets | 18.832 | 17.116 | 14.224 | 13.731 | 11.135 | 16.288 | 8.248 |
| fordyn | Sureau | single | 18.604 | 16.090 | 14.312 | 14.000 | 11.076 | 16.784 | 8.535 |
| fordyn | Sureau | dual | 26.987 | 28.045 | 22.149 | 21.737 | 19.623 | 31.448 | 13.996 |
