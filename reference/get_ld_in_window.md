# Get linkage disequilibrium

Get LD from a group of snps or a single tag snp to all snps in a window.

## Usage

``` r
get_ld_in_window(
  qtl.df = NULL,
  tag.snp = NULL,
  window,
  plink.path = NULL,
  geno.bed,
  in.dir,
  out.dir = NULL,
  verbose = TRUE
)
```

## Arguments

- qtl.df:

  data.frame, table that includes list of snps to calculate LD to with
  columns (CHR, POS, LOGPVAL), corresponding to (chromosome, physical
  position, and -log10(p-value)). QTL are typically defined as hits
  grouped by LD by something like `plink --clump`

- tag.snp:

  character, marker.ID of snp around which to calculate LD. In the form
  'CHR-POS'

- window:

  numeric, a physical distance to determine which snps to calculate LD
  for. Snps within this distance to the tag.snp or qtl.df snps will be
  calculated. The entire size of the region centered on the tag snp
  within which LD will be calculated is 2\*window.

- plink.path:

  character, optional, path to plink2 executable. Will overide option
  set by
  [set_plink_path](https://danforthcenter.github.io/panvaR/reference/set_plink_path.md).

- geno.bed:

  character, prefix of genotype files in plink (bed/bim/fam) format. Do
  not include ".bed" extension.

- in.dir:

  character, directory where genotype files are located

- out.dir:

  character, where to output some temporary files.

- verbose:

  boolean, if TRUE, output some status reports

## Value

Named list with 2 items

- table: table with marker.IDs (CHR-POS) and maximum LD in R2 for each
  snp to the snps in the qtl.df or LD to tag.snp

- key.snp: marker.ID corresponding to the middle of the window, either
  the max(LOGPVAL) of qtl.df or tag.snp. useful to retain for downstream
  functions.

- qtl.snps: marker.ID corresponding to the markers in the qtl.df. useful
  to retain for downstream functions.

## Details

When supplying a group a snps, the ld retained for a given snp is the
maximum LD of the given snp to the snps in the group.

## Examples

``` r
# work in progress
```
