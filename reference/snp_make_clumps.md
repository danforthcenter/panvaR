# Assign snps to clumps based on LD after gwas

Assign snps to clumps based on LD after gwas

## Usage

``` r
snp_make_clumps(
  geno.bed.filename,
  geno.bed.dir,
  gwas.res,
  pvals.in.log = T,
  window = 500,
  ld.thresh = 0.5,
  plink.path = NULL,
  out.dir = NULL
)
```

## Arguments

- geno.bed.filename:

  character, filename of genotype bedfile, no .bed extension

- geno.bed.dir:

  character, path to directory where bed/bim/fam files exist, include
  trailing "/"

- gwas.res:

  data.frame, table of gwas results with columns (marker.ID, CHR, POS,
  PVAL)

- pvals.in.log:

  boolean, are pvalues in `gwas.res` in -log10(p) or not

- window:

  integer, window in kilobases either side of snp to look for snps in
  LD,

- ld.thresh:

  numeric, R2 threshold above which snps will be grouped

- plink.path:

  path to plink 2 executable

- out.dir:

  character, path to a temporary directory to output some intermediate
  files

## Value

table with columns marker.ID and clump_num. Clump_num indicates
groupings, numberings start from the larges pvalue to smallest. May want
to reassign afterwards to be along the genome.

## Examples

``` r
# work in progress
```
