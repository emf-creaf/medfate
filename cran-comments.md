## Minor change submission

This submission (v5.0.0) several improvements, including new dependencies.

## Tested environments

* local R installation (Arch Linux), R 4.6.0
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

`medfate` has a reverse dependency with `medfateland`, under our responsibility, whose current CRAN version will not pass checks with `medfate` v. 5.0.0. 
A new version of `medfateland` complying with `medfate` ver. 5.0.0 will be submitted immediately after this submission is accepted.
