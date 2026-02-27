# Plot genes and their locations

Use a table of genes, descriptions and their start and end points to
plot their locations and annotations along the genome.

## Usage

``` r
plotly_gene_annotation(
  panvar.table.list = NULL,
  annotation.table = NULL,
  middle.snp = NULL,
  window,
  gene.color = "blue",
  point.color = NULL,
  point.fill.scale = NULL
)
```

## Arguments

- panvar.table.list:

  list, output from
  [make_panvar_tables](https://danforthcenter.github.io/panvaR/reference/make_panvar_tables.md).
  Provide either this list or both annotation.table and middle.snp.

- annotation.table:

  table with annotations with columns (geneID, CHR, start, end,
  annotation). start and end correspond to base-pair coordinates of
  start and end of gene. CHR is chromosome of gene.

- middle.snp:

  character, SNP name in form "CHR-POS" center of window. Often the
  key.snp output of
  [get_ld_in_window](https://danforthcenter.github.io/panvaR/reference/get_ld_in_window.md)

- window:

  integer, kilobases on either side of middle.snp to plot

- gene.color:

  character, color to plot genes

- point.color:

  character, variable in annotation.table that indicates how to color
  points plotted next to gene descriptions. If not supplied, no points
  are plotted. The input "LD" is reserved to give functionality to
  [plot_panvar](https://danforthcenter.github.io/panvaR/reference/plot_panvar.md).
  If used, legend will not be displayed.

- point.fill.scale:

  ggplot2 scale object, a fill scale to customize how point.color is
  displayed.

## Value

GGplot object

## Examples

``` r
# Work in progress
```
