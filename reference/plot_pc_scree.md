# Calculate and plot pcs

Calculate and plot pcs

## Usage

``` r
plot_pc_scree(bedfile.path, num.pcs = 10, cumulative = F)
```

## Arguments

- bedfile.path:

  character, path to bedfile (include .bed extension)

- num.pcs:

  numeric, total number of pcs to plot

- cumulative:

  boolean, if T will plot the cumulative proportion of variance

## Value

GGplot of the given plot.

The output is sensitive to the number of PC's indicated. The program
used to calculate PCs returns singular values only for components
indicated. So the proportion can only calculate using the variance
contained in these components. Consequently, shouldn't set the num.pcs
too low, this could lead to strange results.

## Examples

``` r
# work in progress
```
