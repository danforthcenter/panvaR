# panvaR

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16578324.svg)](https://doi.org/10.5281/zenodo.16578324)

An R package for fine-mapping and prioritizing candidate genes from GWA
studies.

## Installation

Ensure remotes package is installed.

``` r

install.packages("remotes")
```

Install panvaR:

``` r

remotes::install_github("danforthcenter/panvaR")
```

We have also implemented a GUI for visualizing results. To use the GUI
run:

``` r

library(panvaR)
panvar_shiny()
```

## Getting started

To get started, check out the tutorial vignette
[here](https://danforthcenter.github.io/panvaR/articles/panvaR.html).

## Main output

PanvaR visualizes specfic QTL by creating information dense plots using
a a flexible array of data sources.

![](https://raw.githubusercontent.com/danforthcenter/panvaR/refs/heads/main/man/figures/panvar_explainer_figure_larger.png)

## Reporting issues

Please let us know if there are problems with the package via the github
issues page. We are always open to improving the ease of use and
functionality so also please inform us of any potential additions or
features. Thanks for using our tool!
