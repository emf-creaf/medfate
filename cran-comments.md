This is a re-submission of package 'medfate', which was formely present in CRAN but got
archived due to an unsolved compilation warning (https://cran.r-project.org/web/packages/medfate/index.html) under clang-UBSAN platform. 

We solved that issue and would like to re-submit the package.

## Tested environments

* local R installation (Arch Linux), R 4.2.2
* windows-latest (on github actions), R release
* macOS-latest (on github actions), R release
* ubuntu-latest (on github actions), R release
* ubuntu-latest (on github actions), R devel
* win-builder (release)
* win-builder (devel)
* win-builder (old-rel)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (on rhub)
* Fedora Linux, R-devel, clang, gfortran (on rhub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (on rhub)

## R CMD check results

In all CI tests 1 NOTE is generated about the size of 'libs', due to the large
amount of compiled code:

     installed size is 41.7Mb
     sub-directories of 1Mb or more:
       libs  40.4Mb

In some platforms (fedora, R-level, clang) we also get 1 NOTE about elapsed time for examples:

    Examples with CPU (user + system) or elapsed time > 5s
               user system elapsed
    growth    4.111  0.004  11.074
    spwb      2.151  0.005   6.484
    plot.spwb 2.100  0.000   6.131
