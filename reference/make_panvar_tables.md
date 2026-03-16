# Make some standardized tables from gwas and annotation tables

Make some standardized tables from gwas and annotation tables

## Usage

``` r
make_panvar_tables(
  gwas.res,
  qtl.df = NULL,
  tag.snp = NULL,
  annotation.table,
  plink.path = NULL,
  pvals.in.log = T,
  geno.bed.filename,
  geno.bed.directory = "/.",
  temp.dir = tempdir(),
  window,
  snp.to.gene.vars = "LD",
  snp.to.gene.buffer = 0,
  compute.scores = F,
  score.vars = NULL,
  score.dirs = NULL,
  score.weights = NULL
)
```

## Arguments

- gwas.res:

  data.frame of all gwas results, should contain columns (CHR, POS,
  PVAL), corresponding to (chromosome, physical position, and pvalue).

- qtl.df:

  data.frame, table that includes list of snps to calculate LD to with
  columns (CHR, POS, LOGPVAL), corresponding to (chromosome, physical
  position, and -log10(p-value)). QTL are typically defined as hits
  grouped by LD by something like `plink --clump`. See
  [get_ld_in_window](https://danforthcenter.github.io/panvaR/reference/get_ld_in_window.md)

- tag.snp:

  character, marker.ID of snp around which to calculate LD. In the form
  'CHR-POS'

- annotation.table:

  table with annotations with columns (geneID, CHR, start, end,
  annotation). start and end correspond to base-pair coordinates of
  start and end of gene. CHR is chromosome of gene.

- plink.path:

  character, optional, path to plink2 executable. Will overide option
  set by
  [set_plink_path](https://danforthcenter.github.io/panvaR/reference/set_plink_path.md).

- pvals.in.log:

  boolean, if TRUE PVAL column has already been converted to
  -log10(pvalue)

- geno.bed.filename:

  character, prefix of genotype files in plink (bed/bim/fam) format. Do
  not include ".bed" extension.

- geno.bed.directory:

  character, directory where genotype files are located

- temp.dir:

  character, where to output some temporary files.

- window:

  numeric, total window size in KB, all variants within .5 \* window are
  calculated.

- snp.to.gene.vars:

  character, numeric variables in gwas.res to aggregate by gene. For
  each gene, snps with a physical position with the start and end of the
  gene are considered. The maximum value for all snps within the gene is
  returned. Special values, `DIST`, `LD` and `LOGPVAL` can be included
  in addition to any user supplied variables.

- snp.to.gene.buffer:

  numeric, kilobases to add to gene start and end to include genes that
  are close but not in gene. Snpeff uses 5 KB by default to call a snp
  "upstream"/"downstream" variant. default is 0.

- compute.scores:

  boolean, if TRUE, snp scores will be computed. See details for more
  info.

- score.vars:

  character, vector of column names indicating which variables to
  included in the score. If compute.scores is TRUE and score.vars is
  NULL, the default score will use equally weighted variables: "DIST",
  "LOGPVAL", "LD".

- score.dirs:

  numeric, a vector indicating which direction is to be considered more
  indicative of an association. 1 indicates higher is better, -1
  indicates lower is better. The order should correspond with the order
  in score.vars

- score.weights:

  numeric, a vector indicating weights for the variables. These must add
  up to 1.

## Value

A named list with the following entries:

- gwas: formatted gwas results.

- anno: formatted annotation results.

- key.snp: tag.snp or in the qtl.df case the highest p-value snp
  supplied to the function for downstream use.

- qtl.df: qtl.df supplied to the program if used.

## Details

Scores: Scores are simple scaled and weighted averages of some
variables. First variables are normalized using min/max normalization.
The variables are then made negative if they need to be reversed to
indicate a larger value as a more desirable value.

For example, distance from the key snp should be reversed as a small
distance is more desirable. A log-pvalue is already of this form 'bigger
is better' so does not need to be altered.

Finally, a weighted average is taken based on user defined weights. The
default weights all variables equally. The outcome is a score from 0-1
that ranks the snps based on these variables. see:
[make_scores](https://danforthcenter.github.io/panvaR/reference/make_scores.md)

## Examples

``` r
# organize options
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

# run function
tables <- make_panvar_tables(
  gwas.res = gwas.df,
  tag.snp = tag.snp,
  annotation.table = annotation.table,
  plink.path = plink.path,
  pvals.in.log = F,
  geno.bed.filename = geno.bed.filename,
  geno.bed.directory = geno.bed.directory,
  window = 25,
  temp.dir = temp.dir,
  compute.scores = FALSE,
  snp.to.gene.buffer = 0)
#> Calculating LD
#> Generating snp to gene correspondence

# snp level results
head(tables$gwas)
#> # A tibble: 6 × 12
#> # Rowwise: 
#>   marker.ID   CHR     POS A1    A2      MAF   EFF    SE     PVAL LOGPVAL     LD
#>   <chr>     <dbl>   <dbl> <chr> <chr> <dbl> <dbl> <dbl>    <dbl>   <dbl>  <dbl>
#> 1 5-6833177     5 6833177 G     A     0.267 1.11  0.145 6.98e-13  12.2   0.923 
#> 2 5-6833238     5 6833238 A     G     0.314 0.769 0.118 5.58e-10   9.25  0.716 
#> 3 5-6833253     5 6833253 C     T     0.230 0.816 0.128 1.24e- 9   8.90  0.0762
#> 4 5-6834607     5 6834607 C     T     0.423 0.757 0.258 3.76e- 3   2.42  0.278 
#> 5 5-6834854     5 6834854 G     A     0.472 0.202 0.202 3.19e- 1   0.496 0.321 
#> 6 5-6835212     5 6835212 C     T     0.402 1.26  0.302 4.46e- 5   4.35  0.268 
#> # ℹ 1 more variable: genes_near_snp <chr>
# gene level results
head(tables$anno)
#> # A tibble: 6 × 7
#> # Rowwise: 
#>     CHR geneID           start     end annotation          dist.from.snp      LD
#>   <int> <chr>            <int>   <int> <chr>                       <dbl>   <dbl>
#> 1     5 Sevir.5G085300 6829932 6832531 (1 of 2) PTHR20961…         24514 NA     
#> 2     5 Sevir.5G085350 6837639 6838969 No gene descriptio…         18076  0.953 
#> 3     5 Sevir.5G085800 6867108 6873803 (1 of 1) KOG4467 -…         10063  0.429 
#> 4     5 Sevir.5G085400 6847970 6850236 (1 of 1) PTHR10641…          6809  0.0960
#> 5     5 Sevir.5G085700 6866196 6868255 (1 of 1) PTHR34543…          9151  0.0471
#> 6     5 Sevir.5G085500 6859612 6862290 (1 of 2) PTHR33146…          2567  0.119 

# clean up
unlink(temp.dir, recursive = TRUE)
```
