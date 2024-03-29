
<!-- README.md is generated from README.Rmd. Please edit that file -->

# extreme.trawl

<!-- badges: start -->

[![codecov](https://codecov.io/gh/valcourgeau/extreme-trawl/branch/master/graph/badge.svg?token=DHLBP4BLP0)](https://codecov.io/gh/valcourgeau/extreme-trawl)
[![R-CMD-check](https://github.com/valcourgeau/extreme-trawl/workflows/R-CMD-check/badge.svg)](https://github.com/valcourgeau/extreme-trawl/actions)
<!-- badges: end -->

The goal of `extreme.trawl` is to provide statistical inference tools
for univariate extreme time series. It comprises pairwise and
generalised method of moments inference strategies described in Courgeau
and Veraart (2021) <doi:10.2139/ssrn.3527739> and simulation schemes
based on Gamma-distributed trawl processes.

## Installation

You can install the released version of extreme.trawl from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("extreme.trawl")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("valcourgeau/extreme-trawl")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(extreme.trawl)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```
