# Make scores

Make scores

## Usage

``` r
make_scores(input.df, cols, directions, weights = NULL)
```

## Arguments

- input.df:

  data.frame, table with only snps that are to be included. Must have
  "marker.ID" column

- cols:

  character, vector of column names indicating which variables to
  included in the score

- directions:

  numeric, a vector indicating which direction is to be considered more
  indicative of an association. 1 indicates higher is better, -1
  indicates lower is better. The order should correspond with the order
  in cols.

- weights:

  numeric, a vector indicating weights for the variables. These must add
  up to 1.

## Value

A table with a score for each marker based on aggregating the variables
indicated.

1.  Variables are first made negative if indicated by `directions`
    vector.

2.  Min/max normalization is applied to put all variables on same scale.

3.  A weighted mean based on the `weights` vector is taken of the
    normalized variables and output as the final score.

## Examples

``` r
# work in progress
```
