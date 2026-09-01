## Minor change submission

This submission (ver. 5.1.0) some improvements and corrections of installation issues of ver 5.0.0.

## Tested environments

* local R installation (Arch Linux), R 4.6.1
* windows-latest (on github actions), R release
* macOS-latest (on github actions), R release
* ubuntu-latest (on github actions), R release
* ubuntu-latest (on github actions), R devel
* ubuntu-latest (on github actions), R oldrel-1
* valgrind (via rhub package), R-devel, Fedora Linux 38
* clang18 (via rhub package), R-devel, Ubuntu 22.04.5 LTS
* clang-asan (via rhub package), R-devel, Ubuntu 22.04.5 LTS

## R CMD check results

In all CI tests only NOTEs are generated, except in clang18 flavor were installation problems occurred. Some memory leaks have been detected with valgrind, which we will try to solve.

## Reverse/Downstream dependencies

`medfate` has a reverse dependency with `medfateland`, under our responsibility, whose current CRAN version should pass checks with `medfate` v. 5.1.0.
