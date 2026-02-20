# Make standard inputs for panvaR

Does some filtering of the genotype file for samples that have a
phenotype and also minor allele frequency and snps with high missing
rates. Will also calculate prinicpal components and optionally kinship
matrix of the genotype file for use downstream in gwas.

## Usage

``` r
make_panvar_inputs(
  genotype.path,
  phenotype.path,
  min.maf = 0.05,
  max.missing.snp = 0.1,
  calc.kinship = F,
  plink.path = NULL,
  out.dir = NULL,
  out.prefix = NULL,
  extra.plink.options = NULL
)
```

## Arguments

- genotype.path:

  character, path to genotype file, supported types: '.bed', .'vcf',
  '.vcf.gz'.

- phenotype.path:

  character, path to table of phenotype to test. Expects samples (lines)
  in column 1 and phenotype in column 2. This is used to determine the
  set of samples (lines) to use in the analysis.

- min.maf:

  numeric, filtering cutoff for minor allele frequency, snps are removed
  if they have maf less than this value. To ignore set to 0.

- max.missing.snp:

  numeric, filtering cutoff for missing rate of snps, snps are removed
  if they have a missing rate higher than this. To ignore set to 1.

- calc.kinship:

  boolean, optional, if TRUE, the kinship matrix will be calculated for
  use in mixed linear model gwas.

- plink.path:

  character, optional, path to plink2 executable. Will overide option
  set by
  [set_plink_path](https://danforthcenter.github.io/panvaR/reference/set_plink_path.md).

- out.dir:

  character, optional, path to store output. Will overide option set by
  [set_out_dir](https://danforthcenter.github.io/panvaR/reference/set_out_dir.md).

- out.prefix:

  character, optional, a prefix for output files. Will overide option
  set by
  [set_panvar_prefix](https://danforthcenter.github.io/panvaR/reference/set_panvar_prefix.md).

- extra.plink.options:

  character, a vector of options to include in call to plink2. Should be
  a vector with plink2 arguments and their values as separate elements
  of vector. E.G. c("–max-maf", ".95", "–max-alleles", "2"). see
  [snp_qc_plink](https://danforthcenter.github.io/panvaR/reference/snp_qc_plink.md)

## Value

Input files to be used for downstream panvaR functions. Stored in
`out.dir` or the option set in
[set_out_dir](https://danforthcenter.github.io/panvaR/reference/set_out_dir.md).
Runs
[snp_qc_plink](https://danforthcenter.github.io/panvaR/reference/snp_qc_plink.md)
to filter for maf and missing using plink2 and then
[rMVP::MVP.Data](https://rdrr.io/pkg/rMVP/man/MVP.Data.html) to prepare
data for GWAS.

## Examples

``` r
# work in progress
```
