
<!-- README.md is generated from README.Rmd. Please edit that file -->
![picture](img/interaction_logo1.jpg)

cytofCore
---------

Description: The package focuses on analysis of CyTOF proteomic data. It is of the pipeline which starts from fcs files and goes to markers and cell type abundance correlation analysis through marker filtratio and cell clusterisation steps. Despite that, there are functions available for interactive R session, it is assumed that cytofCore are run with a graphical interface as R Shiny app GUI.

Installation
------------

Presently you can instal the develop version of this programme:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("AlexanderKononov/cytofCore")
```

Complete list of dependencies is listed in the DESCRIPTION and in the file Library\_louncher.R

In a future, you can install the released version of cytofCore from [Bioconductor](https://www.bioconductor.org) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cytofCore")
```

Run the cytofCore
-----------------

I assumed that this programme will be used with GUI as Shiny App.

``` r
cytofCoreGUI()
```

Example
-------

This is a basic example which shows you how to run the cytofCore:

``` r
library(cytofCore)
#> Warning: replacing previous import 'shiny::renderDataTable' by
#> 'DT::renderDataTable' when loading 'cytofCore'
```

``` r
cytofCoreGUI()
```
