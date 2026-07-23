# Use rMVP to run gwas

Designed to be used with
[make_panvar_inputs](https://danforthcenter.github.io/panvaR/reference/make_panvar_inputs.md).
Will read in files and run GWAS using
[rMVP::MVP](https://rdrr.io/pkg/rMVP/man/MVP.html). Can also supply
inputs as matrices.

## Usage

``` r
panvar_mvp_gwas(
  inputs.dir = NULL,
  in.prefix = NULL,
  npcs = NULL,
  gwas.model = c("GLM", "MLM"),
  output.manhattan = FALSE,
  pheno.mat = NULL,
  geno.mat = NULL,
  map.mat = NULL,
  pcs.mat = NULL,
  kin.mat = NULL,
  out.dir = NULL,
  out.prefix = in.prefix
)
```

## Arguments

- inputs.dir:

  character, directory to find input files created by
  [make_panvar_inputs](https://danforthcenter.github.io/panvaR/reference/make_panvar_inputs.md).
  Consider using `options()$panvar_outdir`.

- in.prefix:

  character, prefix of files in input directory. Same as supplied to
  [make_panvar_inputs](https://danforthcenter.github.io/panvaR/reference/make_panvar_inputs.md).
  Consider using
  [set_panvar_prefix](https://danforthcenter.github.io/panvaR/reference/set_panvar_prefix.md).

- npcs:

  numeric, number of principal components to be included in the gwas
  model.

- gwas.model:

  character, one of "GLM" or "MLM" to refer to a generalized linear
  model and mixed-linear model respectively.

- output.manhattan:

  boolean, if TRUE, visualizations of gwas results will be output to
  out.dir

- pheno.mat:

  matrix, optional, object to use for phenotype

- geno.mat:

  big.matrix, optional, object to use for genotype

- map.mat:

  matrix, optional, object to use a genotype map file

- pcs.mat:

  big.matrix, optional, object to use as principal component matrix

- kin.mat:

  big.matrix, optional, object to use as kinship matrix

- out.dir:

  character, optional, path to store output. Will overide option set by
  [set_out_dir](https://danforthcenter.github.io/panvaR/reference/set_out_dir.md)

- out.prefix:

  character, optional, a prefix for output files. Will overide option
  set by
  [set_panvar_prefix](https://danforthcenter.github.io/panvaR/reference/set_panvar_prefix.md).

## Value

outputs table of gwas results and optionally visualizations produced by
[rMVP::MVP](https://rdrr.io/pkg/rMVP/man/MVP.html)

## References

Lilin Yin, Haohao Zhang, Zhenshuang Tang, Jingya Xu, Dong Yin, Zhiwu
Zhang, Xiaohui Yuan, Mengjin Zhu, Shuhong Zhao, Xinyun Li, Xiaolei Liu,
rMVP: A Memory-Efficient, Visualization-Enhanced, and
Parallel-Accelerated Tool for Genome-Wide Association Study, Genomics,
Proteomics & Bioinformatics, Volume 19, Issue 4, August 2021, Pages
619–628, https://doi.org/10.1016/j.gpb.2020.10.007

## Examples

``` r
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

# make some inputs
make_panvar_inputs(
  plink.path = plink.path,
  genotype.path = genotype.path,
  phenotype.path = phenotype.path,
  out.prefix = "example",
  out.dir = out.dir)
#> Removed 0 samples due to NA values in phenotype.
#> [1] "/tmp/RtmpTrYP6u/panvar_ex/example_PlinkQC_maf0.05_missing0.1"
#> PLINK v2.0.0-a.7.1LM AVX2 Intel (4 May 2026)        cog-genomics.org/plink/2.0/
#> (C) 2005-2026 Shaun Purcell, Christopher Chang    GNU General Public License v3
#> Logging to /tmp/RtmpTrYP6u/panvar_ex/example_PlinkQC_maf0.05_missing0.1.log.
#> Options in effect:
#>   --allow-extra-chr
#>   --bfile /home/runner/work/_temp/Library/panvaR/extdata/Setaria_shattering_example_pruned
#>   --geno 0.1
#>   --keep /tmp/RtmpTrYP6u/panvar_ex/Panvar_list.of.samples.with.phenotype_shattering.txt
#>   --maf 0.05
#>   --make-bed
#>   --out /tmp/RtmpTrYP6u/panvar_ex/example_PlinkQC_maf0.05_missing0.1
#>   --set-all-var-ids @-#
#> 
#> Start time: Thu Jul 23 15:57:38 2026
#> 15989 MiB RAM detected, ~14504 available; reserving 7994 MiB for main
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
#> Writing /tmp/RtmpTrYP6u/panvar_ex/example_PlinkQC_maf0.05_missing0.1.fam ...
#> done.
#> Writing /tmp/RtmpTrYP6u/panvar_ex/example_PlinkQC_maf0.05_missing0.1.bim ...
#> done.
#> Writing /tmp/RtmpTrYP6u/panvar_ex/example_PlinkQC_maf0.05_missing0.1.bed ...
#> 0%done.
#> End time: Thu Jul 23 15:57:38 2026
#> QC was successful, output stored at /tmp/RtmpTrYP6u/panvar_ex/example_PlinkQC_maf0.05_missing0.1
#> Using rMVP to calculate PC's.
#> Preparing data for MVP...
#> Reading file...
#> Preparation for MAP data is done within 0s 
#> inds: 215    markers: 5158
#> Loading genotype at a step of 10000...
#> Preparation for GENOTYPE data is done within 0s 
#> 215 common individuals between phenotype and genotype. 
#> Preparation for PHENOTYPE data is Done within 0s 
#> No NA in genotype, imputation has been skipped.
#> MVP data prepration accomplished successfully!
  
# run gwas
panvar_mvp_gwas(
  inputs.dir = out.dir,
  in.prefix = "example",
  npcs = 2,
  gwas.model = "GLM",
  output.manhattan = FALSE,
  out.dir = out.dir,
  out.prefix = "GWAS.Example")
#> Searching for prefix: GWAS.Example in directory: /tmp/RtmpTrYP6u/panvar_ex
#> Found the following files: 
#> example.geno.bin, 
#> example.geno.desc, 
#> example.geno.ind, 
#> example.geno.map, 
#> example.phe, 
#> example_PlinkQC_maf0.05_missing0.1.bed, 
#> example_PlinkQC_maf0.05_missing0.1.bim, 
#> example_PlinkQC_maf0.05_missing0.1.fam, 
#> example_PlinkQC_maf0.05_missing0.1.log
#> Reading in inputs.
#> Checking phenotype table.
#> Running GWAS
#> ======================== Welcome to MVP =========================
#>     an R package for Memory-efficient, Visualization-enhanced    
#>      and Parallel-accelerated genome-wide association study      
#>                       __  __  __   __  ___                       
#>                      |  \/  | \ \ / / | _ \                      
#>                      | |\/| |  \ V /  |  _/                      
#>                      |_|  |_|   \_/   |_|        Version: 1.4.6  
#>   Design & Maintain: Lilin Yin, Haohao Zhang, and Xiaolei Liu    
#>   Contributors: Zhenshuang Tang, Jingya Xu, Dong Yin, Zhiwu      
#>   Zhang, Xiaohui Yuan, Mengjin Zhu, Shuhong Zhao, Xinyun Li      
#>   Mailto: xiaoleiliu@mail.hzau.edu.cn, ylilin@mail.hzau.edu.cn   
#> =================================================================
#> Start: 2026-07-23 15:57:38 UTC 
#> Input data has 215 individuals and 5158 markers 
#> Markers are detected to be stored by column 
#> Analyzed trait: shattering 
#> Number of threads used: 4 
#> OpenBLAS Library is detected, nice job!
#> Calculate allele frequency... 
#> Computing GRM for 215 individuals using 5158 markers with a step of 10000 
#> Deriving relationship matrix successfully 
#> Eigen Decomposition on GRM 
#> Deriving PCs successfully 
#> Number of PCs included in GLM: 2 
#> -------------------------GWAS Start------------------------- 
#> General Linear Model (GLM) Start... 
#> scanning...
#> Genomic inflation factor (lambda): 10.2215 
#> Significant level: 9.69e-06 
#> End: 2026-07-23 15:57:38 UTC 
#> Total running time: 0s 
#> ===================== MVP ACCOMPLISHED ===================== 
  
# see what we did 
list.files(out.dir)
#>  [1] "GWAS.Example_GLM_GWASresults.csv"                    
#>  [2] "Panvar_list.of.samples.with.phenotype_shattering.txt"
#>  [3] "example.geno.bin"                                    
#>  [4] "example.geno.desc"                                   
#>  [5] "example.geno.ind"                                    
#>  [6] "example.geno.map"                                    
#>  [7] "example.phe"                                         
#>  [8] "example_PlinkQC_maf0.05_missing0.1.bed"              
#>  [9] "example_PlinkQC_maf0.05_missing0.1.bim"              
#> [10] "example_PlinkQC_maf0.05_missing0.1.fam"              
#> [11] "example_PlinkQC_maf0.05_missing0.1.log"              

# clean up 
unlink(out.dir, recursive = TRUE)
```
