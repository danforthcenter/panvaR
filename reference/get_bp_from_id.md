# Extract position from marker.ID in the form "CHR-POS"

Extract position from marker.ID in the form "CHR-POS"

## Usage

``` r
get_bp_from_id(marker.ID, sep = "-")
```

## Arguments

- marker.ID:

  character, marker.ID with chromosome and position separated by some
  seperator

- sep:

  character, character that separates chromosome and position in
  marker.ID

## Value

position for this marker as a number

## Examples

``` r
my.marker <- "Chr01-123"
get_bp_from_id(my.marker, sep = "-")
#> [1] 123
```
