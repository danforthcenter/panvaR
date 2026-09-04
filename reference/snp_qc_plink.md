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

## References

Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015)
Second-generation PLINK: rising to the challenge of larger and richer
datasets. GigaScience, 4.

## Examples

``` r
# get some inputs
plink.path <- bigsnpr::download_plink2()
temp.dir <- file.path(tempdir(), "panvar_ex")
dir.create(temp.dir, showWarnings = FALSE)
genotype.path <- system.file("extdata", "Setaria_shattering_example_pruned.bed", package="panvaR")

# run function
snp_qc_plink(
  genotype.path = genotype.path,
  plink.path = plink.path,
  out.dir = temp.dir,
  out.prefix = "Example")
#> [1] "/tmp/RtmpyxWe1J/panvar_ex/Example_PlinkQC_maf0.05_missing0.1"
#> PLINK v2.0.0-a.7.4LM AVX2 Intel (18 Aug 2026)       cog-genomics.org/plink/2.0/
#> (C) 2005-2026 Shaun Purcell, Christopher Chang    GNU General Public License v3
#> Logging to /tmp/RtmpyxWe1J/panvar_ex/Example_PlinkQC_maf0.05_missing0.1.log.
#> Options in effect:
#>   --allow-extra-chr
#>   --bfile /home/runner/work/_temp/Library/panvaR/extdata/Setaria_shattering_example_pruned
#>   --geno 0.1
#>   --maf 0.05
#>   --make-bed
#>   --out /tmp/RtmpyxWe1J/panvar_ex/Example_PlinkQC_maf0.05_missing0.1
#>   --set-all-var-ids @-#
#> 
#> Start time: Fri Sep  4 19:10:04 2026
#> 15989 MiB RAM detected, ~14580 available; reserving 7994 MiB for main
#> workspace.
#> Using up to 4 compute threads.
#> 598 samples (0 females, 0 males, 598 ambiguous; 598 founders) loaded from
#> /home/runner/work/_temp/Library/panvaR/extdata/Setaria_shattering_example_pruned.fam.
#> Note: 1 nonstandard chromosome code present.
#> 7715 variants loaded from
#> /home/runner/work/_temp/Library/panvaR/extdata/Setaria_shattering_example_pruned.bim.
#> Note: No phenotype data present.
#> Calculating allele frequencies... 0%done.
#> --geno: 0 variants removed due to missing genotype data.
#> 2354 variants removed due to allele frequency threshold(s)
#> (--maf/--max-maf/--mac/--max-mac).
#> 5361 variants remaining after main filters.
#> Writing /tmp/RtmpyxWe1J/panvar_ex/Example_PlinkQC_maf0.05_missing0.1.fam ...
#> done.
#> Writing /tmp/RtmpyxWe1J/panvar_ex/Example_PlinkQC_maf0.05_missing0.1.bim ...
#> done.
#> Writing /tmp/RtmpyxWe1J/panvar_ex/Example_PlinkQC_maf0.05_missing0.1.bed ...
#> 0%done.
#> End time: Fri Sep  4 19:10:04 2026
#> QC was successful, output stored at /tmp/RtmpyxWe1J/panvar_ex/Example_PlinkQC_maf0.05_missing0.1
#> [1] "/tmp/RtmpyxWe1J/panvar_ex/Example_PlinkQC_maf0.05_missing0.1"

# see what we did
list.files(temp.dir)
#> [1] "Example_PlinkQC_maf0.05_missing0.1.bed"
#> [2] "Example_PlinkQC_maf0.05_missing0.1.bim"
#> [3] "Example_PlinkQC_maf0.05_missing0.1.fam"
#> [4] "Example_PlinkQC_maf0.05_missing0.1.log"

# clean up
unlink(temp.dir, recursive = TRUE)
```
