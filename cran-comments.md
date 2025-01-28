## Major change submission

This submission (v4.8.0) incorporates code optimization.

## Tested environments

* local R installation (Arch Linux), R 4.4.1
* windows-latest (on github actions), R release
* macOS (on github actions), macos-13 on GitHub
* macOS-latest (on github actions), R release
* ubuntu-latest (on github actions), R release
* ubuntu-latest (on github actions), R devel
* ubuntu-latest (on github actions), R oldrel-1
* valgrind (on github actions), R-devel, Fedora Linux 38
* clang18 (on github actions), R-devel, Ubuntu 22.04.4 LTS
* clang-asan (on github actions), R-devel, Ubuntu 22.04.4 LTS

## R CMD check results

In all CI tests only NOTEs are generated

## Reverse/Downstream dependencies

`medfate` has a reverse dependency with `medfateland`, under our responsibility. 
A new submission of `medfateland` will be made upon acceptance of the new `medfate` release on CRAN.
