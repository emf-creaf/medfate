## Major change submission

This submission (v4.8.3) incorporates some new functions and bug corrections.

## Tested environments

* local R installation (Arch Linux), R 4.5.1
* windows-latest (on github actions), R release
* macOS-latest (on github actions), R release
* ubuntu-latest (on github actions), R release
* ubuntu-latest (on github actions), R devel
* ubuntu-latest (on github actions), R oldrel-1
* valgrind (via rhub package), R-devel, Fedora Linux 38
* clang18 (via rhub package), R-devel, Ubuntu 22.04.5 LTS
* clang-asan (via rhub package), R-devel, Ubuntu 22.04.5 LTS

## R CMD check results

In all CI tests only NOTEs are generated

## Reverse/Downstream dependencies

`medfate` has a reverse dependency with `medfateland`, under our responsibility. 
A new submission of `medfateland` will be made upon acceptance of the new `medfate` release on CRAN.
