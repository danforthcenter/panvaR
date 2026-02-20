# Quality control a snp file using plink2

Does some standard qc of snp files before gwas. Mainly maf and missing
filter. Also filters for single allele snps. Extra calls to plink2 can
be included to filter in more ways.

## Usage

``` r
snp_qc_plink(
  genotype.path,
  min.maf = 0.05,
  max.missing.snp = 0.1,
  sample.list.path = NULL,
  plink.path = NULL,
  out.dir = NULL,
  out.prefix = NULL,
  extra.options = NULL
)
```

## Arguments

- genotype.path:

  character, path to genotype file, supported types: '.bed', .'vcf',
  '.vcf.gz'.

- min.maf:

  numeric, filtering cutoff for minor allele frequency, snps are removed
  if they have maf less than this value. To ignore set to 0.

- max.missing.snp:

  numeric, filtering cutoff for missing rate of snps, snps are removed
  if they have a missing rate higher than this. To ignore set to 1.

- sample.list.path:

  character, optional, path to a list of samples. Samples in file will
  be included. Sample filtering happens before other filtering per
  plink's order of operations.

- plink.path:

  character, optional, path to plink2 executable. If not provided, will
  default to option set by
  [set_plink_path](https://danforthcenter.github.io/panvaR/reference/set_plink_path.md).

- out.dir:

  character, optional, path to output files. If not provided, will
  default to option set by
  [set_out_dir](https://danforthcenter.github.io/panvaR/reference/set_out_dir.md)

- out.prefix:

  character, optional, prefix for files output.

- extra.options:

  character, a vector of options to include in call to plink2. Should be
  a vector with plink2 arguments and their values as separate elements
  of vector. E.G. c("–max-maf", ".95", "–max-alleles", "2")

## Value

filtered bed/bim/bam files stored in out.dir. character string of path
to bed file.

## Examples

``` r
# work in progress
```
