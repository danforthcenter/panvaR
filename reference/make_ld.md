# use plink to calculate ld

use plink to calculate ld

## Usage

``` r
make_ld(plink.path, snp.name, window, bedfile, in.dir, out.dir)
```

## Arguments

- plink.path:

  path to plink2 executable

- snp.name:

  name of snp

- window:

  total window size in KB, all variants within .5 \* window are
  calculated

- bedfile:

  bedfile prefix, no .bed

- in.dir:

  path to location of bed file

- out.dir:

  path to output temp file

## Value

ld values in file named 'ld_out_temp.vcor'
