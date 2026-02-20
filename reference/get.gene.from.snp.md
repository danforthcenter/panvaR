# get genes that a snp is on

get genes that a snp is on

## Usage

``` r
get.gene.from.snp(bp, gene.df, gene.buffer = 0)
```

## Arguments

- bp:

  numeric, physical position of snp

- gene.df:

  data.frame, gene annotation table with start, end and geneID columns.

- gene.buffer:

  numeric, kilobases to add to gene start and end to include genes that
  are close but not in gene. Snpeff uses 5 KB by default to call a snp
  "upstream"/"downstream" variant

## Value

gene ids that a snp is in
