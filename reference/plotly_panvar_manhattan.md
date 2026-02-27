# Make sideways manhattan plot for building locus zoom. Receives output from a single gwas model.

Make sideways manhattan plot for building locus zoom. Receives output
from a single gwas model.

## Usage

``` r
plotly_panvar_manhattan(
  panvar.table.list = NULL,
  gwas.res = NULL,
  ld.list = NULL,
  pvals.in.log = TRUE,
  plot.r2.thresh = 0.2,
  remove.low.ld.points = FALSE,
  window,
  sig.line,
  qualitative.annotation = NULL,
  qualitative.shape.scale = NULL,
  quantitative.annotation = NULL,
  quantitative.fill.scale = NULL
)
```

## Arguments

- panvar.table.list:

  list, output from
  [make_panvar_tables](https://danforthcenter.github.io/panvaR/reference/make_panvar_tables.md).
  Provide either this list or both gwas.res and ld.list.

- gwas.res:

  data.frame of all gwas results, should contain columns (CHR, POS,
  PVAL), corresponding to (chromosome, physical position, and pvalue).

- ld.list:

  list, output of
  [get_ld_in_window](https://danforthcenter.github.io/panvaR/reference/get_ld_in_window.md)

- pvals.in.log:

  boolean, are pvalues in input data.frames in -log10(p)?

- plot.r2.thresh:

  minimum LD with qtl snps to plot snps colored by LD

- remove.low.ld.points:

  boolean, if TRUE, points below `plot.r2.thresh` will not be plotted.

- window:

  numeric, kilobases on either side of top QTL snp to plot

- sig.line:

  numeric, -log10(p) value to draw line on plot

- qualitative.annotation:

  character, column in gwas.res that contains qualitative annotations.
  For example impact grades from snpeff. See
  [format_snpeff_annotations](https://danforthcenter.github.io/panvaR/reference/format_snpeff_annotations.md).
  Will be plotted as shapes. Only accepts up to 5 classes. "IMPACT" and
  "IMPACT_PLUS" are special cases that will have a pre-assigned scale
  used if supplied here.

- qualitative.shape.scale:

  ggplot scale, an object with a stored call to
  ggplot2::scale_shape_manual. More often an output of the function
  [make_consistent_scale](https://danforthcenter.github.io/panvaR/reference/make_consistent_scale.md).

- quantitative.annotation:

  character, column in gwas.res that contains quantitative annotations.
  For example, variant effect scores. Will be plotted as fill to points.

- quantitative.fill.scale:

  character or scale object, either a character indicating the `option`
  parameter passed to ggplot2::scale_fill_viridis_b that alters the
  color scale used. Or a previous call to a ggplot2 fill scale for
  example ggplot2::scale_fill_stepsn.

## Value

Plotly of manhattan plot with points colored by R2. Accepts input

## Examples

``` r
# Work in progress
```
