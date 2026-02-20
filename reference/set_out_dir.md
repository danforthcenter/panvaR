# Set path for output files

Set path for output files

## Usage

``` r
set_out_dir(out_dir)
```

## Arguments

- out_dir:

  character, path to directory

## Value

sets global option "panvar_outdir" for output directory of intermediate
panvar files. If not supplied will be sent to output of 'tempdir()'.

## Examples

``` r
# set_out_dir("/home/user/output_directory")
```
