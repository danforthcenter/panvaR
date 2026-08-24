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

# specify some paths
genotype.path <- system.file(
  "extdata",
  "Setaria_shattering_example_pruned.bed",
  package = "panvaR")

phenotype.path <- system.file(
  "extdata",
  "Setaria_shattering_example_phenotype.tsv",
  package = "panvaR")

plink.path <- bigsnpr::download_plink2()

out.dir <- file.path(tempdir(), "panvar_ex")
dir.create(out.dir, showWarnings = FALSE)

# run the function
make_panvar_inputs(
  plink.path = plink.path,
  genotype.path = genotype.path,
  phenotype.path = phenotype.path,
  out.prefix = "example",
  out.dir = out.dir)
#> Removed 0 samples due to NA values in phenotype.
#> [1] "/tmp/RtmppTSlo7/panvar_ex/example_PlinkQC_maf0.05_missing0.1"
#> PLINK v2.0.0-a.7.4LM AVX2 Intel (18 Aug 2026)       cog-genomics.org/plink/2.0/
#> (C) 2005-2026 Shaun Purcell, Christopher Chang    GNU General Public License v3
#> Logging to /tmp/RtmppTSlo7/panvar_ex/example_PlinkQC_maf0.05_missing0.1.log.
#> Options in effect:
#>   --allow-extra-chr
#>   --bfile /home/runner/work/_temp/Library/panvaR/extdata/Setaria_shattering_example_pruned
#>   --geno 0.1
#>   --keep /tmp/RtmppTSlo7/panvar_ex/Panvar_list.of.samples.with.phenotype_shattering.txt
#>   --maf 0.05
#>   --make-bed
#>   --out /tmp/RtmppTSlo7/panvar_ex/example_PlinkQC_maf0.05_missing0.1
#>   --set-all-var-ids @-#
#> 
#> Start time: Mon Aug 24 15:51:00 2026
#> 15989 MiB RAM detected, ~14459 available; reserving 7994 MiB for main
#> workspace.
#> Using up to 4 compute threads.
#> 598 samples (0 females, 0 males, 598 ambiguous; 598 founders) loaded from
#> /home/runner/work/_temp/Library/panvaR/extdata/Setaria_shattering_example_pruned.fam.
#> Note: 1 nonstandard chromosome code present.
#> 7715 variants loaded from
#> /home/runner/work/_temp/Library/panvaR/extdata/Setaria_shattering_example_pruned.bim.
#> Note: No phenotype data present.
#> --keep: 215 samples remaining.
#> 215 samples (0 females, 0 males, 215 ambiguous; 215 founders) remaining after
#> main filters.
#> Calculating allele frequencies... 0%done.
#> --geno: 0 variants removed due to missing genotype data.
#> 2557 variants removed due to allele frequency threshold(s)
#> (--maf/--max-maf/--mac/--max-mac).
#> 5158 variants remaining after main filters.
#> Writing /tmp/RtmppTSlo7/panvar_ex/example_PlinkQC_maf0.05_missing0.1.fam ...
#> done.
#> Writing /tmp/RtmppTSlo7/panvar_ex/example_PlinkQC_maf0.05_missing0.1.bim ...
#> done.
#> Writing /tmp/RtmppTSlo7/panvar_ex/example_PlinkQC_maf0.05_missing0.1.bed ...
#> 0%done.
#> End time: Mon Aug 24 15:51:00 2026
#> QC was successful, output stored at /tmp/RtmppTSlo7/panvar_ex/example_PlinkQC_maf0.05_missing0.1
#> Using rMVP to calculate PC's.
#> Preparing data for MVP...
#> Reading file...
#> Preparation for MAP data is done within 0s 
#> inds: 215    markers: 5158
#> Loading genotype at a step of 10000...
#> Preparation for GENOTYPE data is done within 0s 
#> 215 common individuals between phenotype and genotype. 
#> Preparation for PHENOTYPE data is Done within 0s 
#> MVP data prepration accomplished successfully!

# created some input files
list.files(out.dir)
#>  [1] "Panvar_list.of.samples.with.phenotype_shattering.txt"
#>  [2] "example.geno.bin"                                    
#>  [3] "example.geno.desc"                                   
#>  [4] "example.geno.ind"                                    
#>  [5] "example.geno.map"                                    
#>  [6] "example.phe"                                         
#>  [7] "example_PlinkQC_maf0.05_missing0.1.bed"              
#>  [8] "example_PlinkQC_maf0.05_missing0.1.bim"              
#>  [9] "example_PlinkQC_maf0.05_missing0.1.fam"              
#> [10] "example_PlinkQC_maf0.05_missing0.1.log"              

# clean up
unlink(out.dir, recursive = TRUE)
```
