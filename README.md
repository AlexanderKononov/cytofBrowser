
<!-- README.md is generated from README.Rmd. Please edit that file -->
cytofCore
=========

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/AlexanderKononov/cytofCore.svg?branch=master)](https://travis-ci.com/AlexanderKononov/cytofCore) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/AlexanderKononov/cytofCore?branch=master&svg=true)](https://ci.appveyor.com/project/AlexanderKononov/cytofCore) [![Codecov test coverage](https://codecov.io/gh/AlexanderKononov/cytofCore/branch/master/graph/badge.svg)](https://codecov.io/gh/AlexanderKononov/cytofCore?branch=master) <!-- badges: end -->

Description: The package focuses on analysis of CyTOF proteomic data. It is of the pipeline which starts from fcs files and goes to markers and cell type abundance correlation analysis through marker filtratio and cell clusterisation steps. Despite that, there are functions available for interactive R session, it is assumed that cytofCore are run with a graphical interface as R Shiny app GUI.

Installation
------------

Presently you can instal the develop version of this programme:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("AlexanderKononov/cytofCore")
```

On the current stage of the project develop there are few dependency packages which don't instal and deploy automatically, such as (flowCore, FlowSOM, ConsensusClusterPlus). I advise to instal they manually:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("flowCore")
BiocManager::install("FlowSOM")
BiocManager::install("ConsensusClusterPlus")
```

Complete list of dependencies is listed in the DESCRIPTION and in the file Library\_louncher.R

In a future, you can install the released version of cytofCore from [Bioconductor](https://www.bioconductor.org) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cytofCore")
```

Example
-------

This is a basic example which shows you how to run the cytofCore:

``` r
library(cytofCore)
#> Warning: replacing previous import 'dplyr::filter' by 'flowCore::filter' when
#> loading 'cytofCore'
#> Warning: replacing previous import 'dplyr::count' by 'matrixStats::count' when
#> loading 'cytofCore'
```

``` r
cytofCoreGUI()
```
