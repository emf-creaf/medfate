## Error fixing submission

This submission (v3.1.4) tries to fix memory access issues arisen from previous submission (v3.1.3) under ASAN and valgrind flavors. We debugged the code using valgrind and we expect the bug to be corrected. Note, however, that we could not test the package on platform Debian Linux, R-devel, GCC ASAN/UBSAN due to pre-processing errors in Rhub.


## Tested environments

* local R installation (Arch Linux), R 4.3.1
* windows-latest (on github actions), R release
* macOS-latest (on github actions), R release
* ubuntu-latest (on github actions), R release
* ubuntu-latest (on github actions), R devel
* ubuntu-latest (on github actions), R oldrel-1
* win-builder (release)
* win-builder (devel)
* win-builder (old-rel)

## R CMD check results

In all CI tests 1 NOTE is generated about the size of 'libs', due to the large
amount of compiled code (e.g.):

     installed size is 41.7Mb
     sub-directories of 1Mb or more:
       libs  40.4Mb

