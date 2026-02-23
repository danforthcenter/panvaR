# panvaR

An R package for prioritizing candidate genes from GWA studies.

## Installation

Ensure remotes package is installed.

``` r
install.packages("remotes")
```

Install panvaR:

``` r
remotes::install_github("danforthcenter/panvaR", build_vignettes = TRUE)
```

We have also implemented a GUI for visualizing results. To use the GUI
run:

``` r
library(panvaR)
panvar_gui()
```

## Getting started

To get started, check out the documentation
[here](https://danforthcenter.github.io/panvaR/).

## Main output

PanvaR visualizes specfic QTL by creating information dense plots using
a a flexible array of data sources.

![](articles/panvar_explainer_figure_larger.png)

## Reporting issues

Please let us know if there are problems with the package via the github
issues page. We are always open to improving the ease of use so please
inform us of any potential improvements. Thanks for using our tool!
