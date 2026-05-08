# Extracts snpeff annotations from vcf

Will take a vcf with snpeff annotations and output a list of two
dataframes.

## Usage

``` r
format_snpeff_annotations(vcfpath)
```

## Arguments

- vcfpath:

  path to vcf file

## Value

A list of two dataframes:

- final_snp_impacts: one row per snp, has the highest impact of any
  transcript for that snp indicated

- formatted_snpeff_annotations: one row per transcript. snps will have
  multiple transcripts identified by snpeff.

The `IMPACT_PLUS` column has had modifier snps broken out into
intergenic and coding modifiers.

## References

"Using Drosophila melanogaster as a model for genotoxic chemical
mutational studies with a new program, SnpSift", Cingolani, P., et. al.,
Frontiers in Genetics, 3, 2012.

## Examples

``` r
# work in progress 
```
