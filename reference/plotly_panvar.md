# Make a plot that lines up manhattan and gene locations

Make a plot that lines up manhattan and gene locations

## Usage

``` r
plotly_panvar(
  panvar.table.list = NULL,
  gwas.res = NULL,
  ld.list = NULL,
  annotation.table = NULL,
  pvals.in.log = T,
  plot.r2.thresh = 0.2,
  remove.low.ld.points = FALSE,
  window,
  sig.line,
  qualitative.annotation = NULL,
  qualitative.shape.scale = NULL,
  quantitative.annotation = NULL,
  quantitative.fill.scale = NULL,
  plot.title = "",
  annotation.point.variable = "LD",
  annotation.point.scale = NULL
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

- annotation.table:

  table with annotations with columns (geneID, CHR, start, end,
  annotation). start and end correspond to base-pair coordinates of
  start and end of gene. CHR is chromosome of gene.

- pvals.in.log:

  boolean, if TRUE PVAL column has already been converted to
  -log10(pvalue)

- plot.r2.thresh:

  minimum LD with qtl snps to plot snps colored by LD

- remove.low.ld.points:

  boolean, if TRUE, points below `plot.r2.thresh` will not be plotted.

- window:

  numeric, total window size in KB, all variants within .5 \* window are
  calculated.

- sig.line:

  numeric, -log10(p) value to draw line on plot

- qualitative.annotation:

  character, column in `gwas.res` that contains qualitative annotations.
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

  character, column in `gwas.res` that contains quantitative
  annotations. For example, variant effect scores. Will be plotted as
  fill to points.

- quantitative.fill.scale:

  character or scale object, either a character indicating the `option`
  parameter passed to ggplot2::scale_fill_viridis_b that alters the
  color scale used. Or a previous call to a ggplot2 fill scale for
  example ggplot2::scale_fill_stepsn.

- plot.title:

  character, a title for the plot

- annotation.point.variable:

  character, variable in `annotation.table` that indicates how to color
  points plotted next to gene descriptions. If not supplied, no points
  are plotted. The input "LD" is reserved and will use LD.

- annotation.point.scale:

  ggplot2 scale object, a color scale to customize how point.color is
  displayed.

## Value

ggplot2 object of plot with manhattan plot alongside genes for a given
genomic window.

## Examples

``` r
# work in progress
```
