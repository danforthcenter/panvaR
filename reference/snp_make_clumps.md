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
tag.snp <- "Chr_05-6857045"
gwas.df <- read.csv(system.file(
    "extdata",
    "PanvarExample_GLM_GWASresults.csv",
    package = "panvaR"))
annotation.table <- read.csv(system.file(
    "extdata",
    "Setaria_shattering_annotation.csv",
    package = "panvaR"))
plink.path <- bigsnpr::download_plink2()
temp.dir <- file.path(tempdir(), "panvar_ex")
dir.create(temp.dir, showWarnings = FALSE)
geno.bed.filename <- "Setaria_shattering_example_pruned.bed"
geno.bed.directory <- system.file("extdata", package="panvaR")

# look at only significant snps 
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
gwas.df_sub <- gwas.df %>% 
  filter(LOGPVAL > 8)

# get clumps
clump.table <- snp_make_clumps(
  geno.bed.filename = geno.bed.filename,
  geno.bed.dir = geno.bed.directory,
  gwas.res = gwas.df_sub,
  pvals.in.log = FALSE,
  window = 500,
  ld.thresh = .5,
  plink.path = plink.path,
  out.dir = temp.dir)
#> Creating clumps...
#>   |                                                                              |=====================================                                 |  52%  |                                                                              |===================================================                   |  73%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
  
 head(clump.table)
#>        marker.ID clump_num
#> 1 Chr_05-6419601         1
#> 2 Chr_05-6463617         1
#> 3 Chr_05-6466490         1
#> 4 Chr_05-6467430         1
#> 5 Chr_05-6474403         1
#> 6 Chr_05-6474887         1

# clean up
unlink(temp.dir, recursive = TRUE)
```
