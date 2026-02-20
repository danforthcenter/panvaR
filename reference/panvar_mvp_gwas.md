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

## Examples

``` r
# work in progress
```
