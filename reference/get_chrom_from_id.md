# Extract chromsome from marker.ID in the form "CHR-POS"

Extract chromsome from marker.ID in the form "CHR-POS"

## Usage

``` r
get_chrom_from_id(marker.ID, sep = "-", tonumeric = TRUE)
```

## Arguments

- marker.ID:

  character, marker.ID with chromosome and position separated by some
  seperator

- sep:

  character, character that separates chromosome and position in
  marker.ID

- tonumeric:

  boolean, should the function attempt to convert the chromosome to
  numeric?

## Value

chromsome for this marker as a number if you want

## Examples

``` r
my.marker <- "Chr01-123"
get_chrom_from_id(my.marker, tonumeric = TRUE)
#> [1] 1
get_chrom_from_id(my.marker, tonumeric = FALSE)
#> [1] "Chr01"
```
