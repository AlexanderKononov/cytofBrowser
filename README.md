
<!-- README.md is generated from README.Rmd. Please edit that file -->

![picture](img/interaction_logo1.jpg)

## cytofBrowser

The package focuses on analysis of CyTOF proteomic data by graphical
interface. It is of the pipeline which starts from FCS files which you
get from the CyTOF machine. It allows preprocessing and transforms the
data then getting the specific cell populations by markers as well as
clustering cells, estimate the abundance of cell populations and find
the correlations and patterns by high-throughput hypothesis testing
algorithm. CytofBrowser supplies the tool to visualize explore and
interact with data on each step of the pipeline.
![picture](img/cytofBrowser_hallmarks.jpg)

## Graphical interface

Despite that, there are functions available for interactive R session,
it is assumed that cytofBrowser are run with a graphical interface as R
Shiny app GUI.

Since I have not devoted time to create detail documentation yet, I note
the link to the video which I sent to my colleague and mentor to get the
feedback. It can look redundantly detail or in an odd format. But there
are enough explanations of work and examples of a graphical interface
using. P.S. Don’t judge my video skills strictly=)

<https://youtu.be/rWSD7y0Gaic>

## Installation

Presently you can instal the develop version of this programme:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("AlexanderKononov/cytofBrowser")
```

## Run the cytofBrowser

I assumed that this programme will be used with GUI as Shiny App.

``` r
cytofBrowser::cytofBrowserGUI()
```

## trubelshutingTroubleshooting

If you can’t install the app correctly: Check the dependencies of
packages. List of demanded packages is waited in file DESCRIPTION. Also
the same list of the packages mentioned in additional file
Library\_louncher.R This files can be used as a script to upload all
dependencies in the running R session manually.

If you faced with a problem of updating any package from dependencies:
Check the R version and Rstudio version (if you use it). After R 4.0
update some packages can not update correctly by the previous R version
and can show “non-zero checksum” error. Some basic packages can have
problem with updating within Rstudio and show the error “can’t remove
the previous version of package”. The most straightforward way to fix
this problem is to delete the target package manually from the directory
with R package libraries (for windows usually Documents/R/win-library)
and there instal it again.

If no one button does not react or collapses the app: One previous
version of shinyFiles for windows is not stable. You can try to manually
update this package.

``` r
install.packages("shinyFiles")
```

or you can try to use the previous version of the package.

``` r
lin <- "https://cran.r-project.org/src/contrib/Archive/shinyFiles/shinyFiles_0.7.5.tar.gz"
install.packages(lin, repos=NULL, type="source")
```
