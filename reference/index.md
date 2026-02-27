# Package index

## Formatting Inputs

“Functions for standardizing inputs and organizing files”

- [`set_out_dir()`](https://danforthcenter.github.io/panvaR/reference/set_out_dir.md)
  : Set path for output files
- [`set_panvar_prefix()`](https://danforthcenter.github.io/panvaR/reference/set_panvar_prefix.md)
  : sets global option "panvar_prefix" for output directory of
  intermediate panvar files. If not supplied will be sent to output of
  'tempdir()'.
- [`set_plink_path()`](https://danforthcenter.github.io/panvaR/reference/set_plink_path.md)
  : Set path to plink2 executable
- [`make_panvar_inputs()`](https://danforthcenter.github.io/panvaR/reference/make_panvar_inputs.md)
  : Make standard inputs for panvaR
- [`snp_qc_plink()`](https://danforthcenter.github.io/panvaR/reference/snp_qc_plink.md)
  : Quality control a snp file using plink2
- [`format_snpeff_annotations()`](https://danforthcenter.github.io/panvaR/reference/format_snpeff_annotations.md)
  : Extracts snpeff annotations from vcf

## GWAS

“Functions for generating information about snps and integrating
results”

- [`get_ld_in_window()`](https://danforthcenter.github.io/panvaR/reference/get_ld_in_window.md)
  : Get linkage disequilibrium
- [`make_panvar_tables()`](https://danforthcenter.github.io/panvaR/reference/make_panvar_tables.md)
  : Make some standardized tables from gwas and annotation tables
- [`make_scores()`](https://danforthcenter.github.io/panvaR/reference/make_scores.md)
  : Make scores
- [`panvar_mvp_gwas()`](https://danforthcenter.github.io/panvaR/reference/panvar_mvp_gwas.md)
  : Use rMVP to run gwas
- [`snp_make_clumps()`](https://danforthcenter.github.io/panvaR/reference/snp_make_clumps.md)
  : Assign snps to clumps based on LD after gwas
- [`get_chrom_from_id()`](https://danforthcenter.github.io/panvaR/reference/get_chrom_from_id.md)
  : Extract chromsome from marker.ID in the form "CHR-POS"
- [`get_bp_from_id()`](https://danforthcenter.github.io/panvaR/reference/get_bp_from_id.md)
  : Extract position from marker.ID in the form "CHR-POS"

## Plotting

“Functions for visualizing results”

- [`plot_effect()`](https://danforthcenter.github.io/panvaR/reference/plot_effect.md)
  : Make volcano style effect size vs pvalue plot
- [`plot_gene_annotation()`](https://danforthcenter.github.io/panvaR/reference/plot_gene_annotation.md)
  : Plot genes and their locations
- [`plot_panvar()`](https://danforthcenter.github.io/panvaR/reference/plot_panvar.md)
  : Make a plot that lines up manhattan and gene locations
- [`plot_panvar_manhattan()`](https://danforthcenter.github.io/panvaR/reference/plot_panvar_manhattan.md)
  : Make sideways manhattan plot for building locus zoom. Receives
  output from a single gwas model.
- [`plot_pc_scree()`](https://danforthcenter.github.io/panvaR/reference/plot_pc_scree.md)
  : Calculate and plot pcs
- [`plotly_gene_annotation()`](https://danforthcenter.github.io/panvaR/reference/plotly_gene_annotation.md)
  : Plot genes and their locations
- [`plotly_panvar()`](https://danforthcenter.github.io/panvaR/reference/plotly_panvar.md)
  : Make a plot that lines up manhattan and gene locations
- [`plotly_panvar_manhattan()`](https://danforthcenter.github.io/panvaR/reference/plotly_panvar_manhattan.md)
  : Make sideways manhattan plot for building locus zoom. Receives
  output from a single gwas model.
- [`make_consistent_scale()`](https://danforthcenter.github.io/panvaR/reference/make_consistent_scale.md)
  : Generate discrete scale

## Shiny

“Functions for running the shiny app”

- [`panvar_shiny()`](https://danforthcenter.github.io/panvaR/reference/panvar_shiny.md)
  : panvar_shiny Run the panvaR shiny implementation. Will launch in a
  default browser or can navigate to IP address shown.
