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
  columns (CHR, POS, PVAL), corresponding to (chromosome, physical
  position, and p-value). QTL are typically defined as hits grouped by
  LD by something like `plink --clump`. See
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
  returned. Special values, `EFF`, `DIST`, `LD` and `LOGPVAL` can be
  included in addition to any user supplied variables. For `EFF` the
  largest absolute value snp effect is returned, retaining the sign.

- snp.to.gene.buffer:

  numeric, kilobases to add to gene start and end to include genes that
  are close but not in gene. Snpeff uses 5 KB by to call a snp
  "upstream"/"downstream" variant. This function's default is 0.

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

Snp to gene: This functionality attempts to group gwas results at the
gen level using maximum value of a given metric by default. Special
cases include: `EFF` (effect size of SNP), `DIST` (distance to tag snp),
`LD` (linkage disequilibrium with tag snp), and `LOGPVAL`. With these
metrics the most impactful direction is chosen instead of just the
maximum. So higher LD and -log(pvalue), smaller distance and largest
absolute value for effect as effect sizes can be negative. “

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
#> Warning: No SNPs found close to the provided tag snp or qtl.df SNPs. Maybe increase window or check gwas results table.

# snp level results
head(tables$gwas)
#>   marker.ID CHR     POS A1 A2       MAF       EFF        SE         PVAL
#> 1 5-6833177   5 6833177  G  A 0.2674419 1.1088102 0.1449118 6.975650e-13
#> 2 5-6833238   5 6833238  A  G 0.3139535 0.7688717 0.1182201 5.582764e-10
#> 3 5-6833253   5 6833253  C  T 0.2302326 0.8164636 0.1284159 1.244639e-09
#> 4 5-6834607   5 6834607  C  T 0.4232558 0.7573602 0.2584663 3.760119e-03
#> 5 5-6834854   5 6834854  G  A 0.4720930 0.2020393 0.2023104 3.191033e-01
#> 6 5-6835212   5 6835212  C  T 0.4023256 1.2573160 0.3015772 4.463099e-05
#>      LOGPVAL        LD
#> 1 12.1564153 0.9228230
#> 2  9.2531507 0.7155770
#> 3  8.9049566 0.0761779
#> 4  2.4247985 0.2776270
#> 5  0.4960687 0.3209310
#> 6  4.3503635 0.2678820
# gene level results
head(tables$anno)
#> # A tibble: 6 × 6
#> # Rowwise: 
#>     CHR geneID           start     end annotation                  dist.from.snp
#>   <int> <chr>            <int>   <int> <chr>                               <dbl>
#> 1     5 Sevir.5G085300 6829932 6832531 (1 of 2) PTHR20961//PTHR20…         24514
#> 2     5 Sevir.5G085350 6837639 6838969 No gene description.                18076
#> 3     5 Sevir.5G085800 6867108 6873803 (1 of 1) KOG4467 - Unchara…         10063
#> 4     5 Sevir.5G085400 6847970 6850236 (1 of 1) PTHR10641//PTHR10…          6809
#> 5     5 Sevir.5G085700 6866196 6868255 (1 of 1) PTHR34543//PTHR34…          9151
#> 6     5 Sevir.5G085500 6859612 6862290 (1 of 2) PTHR33146:SF2 - E…          2567

# clean up
unlink(temp.dir, recursive = TRUE)
```
