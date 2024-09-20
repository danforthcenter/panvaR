#' A function that automates the Panvar analysis when supplies genotype data and phenotype data.
#' 
#' @param phenotype_data_path Path to the phenotype data
#' @param vcf_file_path Path to the VCF file
#' @param r2_threshold The r2 threshold
#' Defaults to 0.6
#' @param maf The minor Allele Frequency
#' Defaults to 0.05
#' @param missing_rate The missing rate filter for your genotype data
#' Defaults to 0.1
#' @param window The window around the tag snp
#' Defaults to 500000
#' @param pc_min (optional) What is the minimum number of PCs that should be included in GWAS?
#' Defaults to 5
#' @param pc_max (optional) What is the maximum number of PCs that should be included in GWAS?
#' Defaults to 5
#' @param dynamic_correlation (optional) Should the PCs, beyond minimum, be calculated dynamically?
#' @param all_impacts (optional) Should all impacts be included in the report?
#' Defaults to FALSE - in which case only "MODERATES" and "HIGH" impacts will be included
#' @examples
#' panvar_pheno_geno_snpeff("<path_to_phenotype_data>", "<path_to_vcf_file>", chrom = "<chorm>", bp = <bp_value>, r2_threshold = 0.6)
#' @export
panvar_pheno_geno_snpeff <- function(phenotype_data_path,vcf_file_path,chrom,bp, r2_threshold = 0.6, maf = 0.05, missing_rate = 0.10, window = 500000,pc_min = 5,pc_max = 5, dynamic_correlation = FALSE, all.impacts = FALSE){

    # Check if the vcf_file has a tbi file
    proper_tbi(vcf_file_path)

    # The end goal of this function is to convieneintly make
    # 1. The plot from Panvar
    # 2. The report table -
    # The report table should have - Gwas pvalues,

    # Run GWAS for the user
    gwas_table_denovo <- panvar_gwas(
        phentotype_path = phenotype_data_path,
        genotype_data = vcf_file_path,
        pc_min = pc_min,
        pc_max = pc_max,
        dynamic_correlation = dynamic_correlation,
        maf = maf,
        missing_rate = missing_rate
    )

    gwas_table <- check_gwas_table(gwas_table_denovo)

    # get the tag snp from the gwas results
    tag_snp <- tag_snp_func(gwas_table_denovo)

    # TODO - this needs some clean up from here
    # but for now we will just use 
    bp = tag_snp$tag_snp_bp
    chrom = tag_snp$tag_snp_chromosome

    # Now we can just call the `panvar_gwas_tagsnp_snpeff` -
    # function to the rest of the analysis.
    panvar_run <-  panvar_gwas_tagsnp_snpeff(
        gwas_table_path = gwas_table,
        vcf_file_path = vcf_file_path,
        chrom = chrom,
        bp = bp,
        r2_threshold = r2_threshold,
        maf = maf,
        missing_rate = missing_rate,
        window = window,
        all.impacts = all.impacts
    )

    return(panvar_run)
}
