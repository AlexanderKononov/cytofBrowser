
<!-- README.md is generated from README.Rmd. Please edit that file -->
cytofCore
=========

<!-- badges: start -->
<!-- badges: end -->
Description: The package focuses on analysis of CyTOF proteomic data. It is of the pipeline which starts from fcs files and goes to markers and cell type abundance correlation analysis through marker filtratio and cell clusterisation steps. Despite that, there are functions available for interactive R session, it is assumed that cytofCore are run with a graphical interface as R Shiny app GUI.

Installation
------------

You can install the released version of cytofCore from [Bioconductor](https://www.bioconductor.org) with:

``` r
install.packages("cytofCore")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(cytofCore)
#> Warning: replacing previous import 'dplyr::filter' by 'flowCore::filter' when
#> loading 'cytofCore'
#> Warning: replacing previous import 'dplyr::count' by 'matrixStats::count' when
#> loading 'cytofCore'
#> Warning: replacing previous import 'flowCore::alias' by 'stats::alias' when
#> loading 'cytofCore'
#> Warning: replacing previous import 'flowCore::filter' by 'stats::filter' when
#> loading 'cytofCore'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'cytofCore'
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

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

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!
