
<!-- README.md is generated from README.Rmd. Please edit that file -->
cytofCore
=========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/AlexanderKononov/cytofCore.svg?branch=master)](https://travis-ci.org/AlexanderKononov/cytofCore) <!-- badges: end -->

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
