# RNWMStat
R package for the National Water Model hydrology evaluation

# Purpose
The purpose of this R package is to evaluate the streamflow simulated 
by the National Water Model (NWM)/WRF-Hydro and to carry out some basic 
statistical analyses.  These tools are both free and open-source, 
just like R, which should help make them accessible and popular.
# Installaton

Installation of RNWMStat package (not on [CRAN](http://cran.r-project.org/)) can be achieved by using devtools package (on CRAN), so devtools should be installed
first. The following is done for the initial install or to update the
rwrfhydro package.

``` r
install.packages("devtools")
devtools::install_github("NCAR/RNWMStat")
```

# Developing and bug reports

Bugs are to be reported
[here](https://github.com/NCAR/RNWMStats/issues). If you want to help
solve bugs and fixes into the code, please continue reading about
developing.

