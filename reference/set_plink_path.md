# Set path to plink2 executable

Set path to plink2 executable

## Usage

``` r
set_plink_path(plink_path)
```

## Arguments

- plink_path:

  character, path to plink2 executable

## Value

sets global option that for path to global executable for plink2. If
plink2 is on path and option is not supplied, will use the executable on
the path.

Accessible at: options()\$plink_path

plink V2.00a6 minimum recommended

## Examples

``` r
# set_plink_path("path_to_executable")
```
