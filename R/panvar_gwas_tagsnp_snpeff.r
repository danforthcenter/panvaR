#' panvar_gwas_tagsnp_snpeff
#' 
#' @param gwas_table_path Path to the GWAS table
#' @param vcf_file_path Path to the VCF file
#' @param chrom The chromosome alphanumeric of the tag SNP
#' @param bp The base pair loci of the tag SNP
#' @param r2_threshold The r2 threshold
#' Defaults to 0.6
#' @param maf The minor Allele Frequency
#' Defaults to 0.05
#' @param missing_rate The missing rate filter for your genotype data
#' Defaults to 0.1
#' @param window The window around the tag snp. 
#' Either supplied with as kilo basese as "#kb" or as bp values.
#' If you want to supply 250 Kilobases you should supply "250kb" or 250000 - whatever is your preference.
#' Defaults to "500kb".
#' @param all_impacts (optional) Should all impacts be included in the report?
#' Defaults to FALSE - in which case only "MODERATES" and "HIGH" impacts will be included
#' 
#' @import tidyverse
#' @import data.table
#' @import sys
#' @import parallel
#' @import bigsnpr
#' @import modelr
#'
#' @export
#' 
#' @examples
#' panvar_gwas_tagsnp_snpeff("<path_to_gwas_table>", "<path_to_vcf_file>", chrom = "<chorm>", bp = <bp_value>, r2_threshold = 0.6)
#' panvar_gwas_tagsnp_snpeff(
#'    "<path_to_gwas_table>", 
#'    "<path_to_vcf_file>", 
#'    chrom = chr_05, 
#'    bp = 54675,
#'    r2_threshold = 0.6
#' )

panvar_gwas_tagsnp_snpeff <- function(gwas_table_path,vcf_file_path,chrom,bp, r2_threshold = 0.6, maf = 0.05, missing_rate = 0.10, window = "500kb", all.impacts = FALSE){

  # convert the window into bp values
  window_bp <- window_unit_func(window)

  # Check if the gwas table has a tbi file
    # Check if the vcf_file has a tbi file
    proper_tbi(vcf_file_path)

    # The end goal of this function is to convieneintly make
    # 1. The plot from Panvar
    # 2. The report table -
    # The report table should have - Gwas pvalues,

    # Read the gwas table
    # TODO - This need to go to input verification
    gwas_table <- check_gwas_table(gwas_table_path)

    # convert the vcf file to plink format
    in_plink_format <- vcf_to_plink2(vcf_file_path)

    # clean up the supplied vcf file
    cleaned_up <- bed_file_clean_up(in_plink_format$bed, maf = maf, missing_rate = missing_rate)

    # subset your genotype data around the tag snp
    subset_genotype_data <- subset_around_tag(cleaned_up,chrom = chrom, bp = bp, window = window_bp)

    # using ld get the list of bps to keep
    table <- ld_filtered_snp_list(subset_genotype_data,chrom = chrom, bp = bp, r2_threshold = r2_threshold)

    # Make the LD table
	# This table has the CHROM, BP and LD values and re-adds the tag SNP
    # This is one of the tables that will be left_joined later
    ld_table <- ld_table_maker(table)

    # Convert the table into the list of SNPs to keep
    # such that BCFtools can filter them out of a the original VCF file
    keep_snp_list <- snps_to_keep(table)

    # The way plink2 converts vcf files may lead to - 
    # a situation where the names in plink2 -
    # and BCFtools have different chromosomes names. -
    # This needs to be accounted for.

	plink2_bcf_dictionary <- plink2_bcftools_chroms_dictionary(vcf_file_path,in_plink_format$bim)

	if(!is.null(plink2_bcf_dictionary)){
		ld_table_checked <- apply_dict(plink2_bcf_dictionary, ld_table)

		snp_keep_list_checked <- apply_dict(plink2_bcf_dictionary, keep_snp_list)
	} else{
		
		ld_table_checked <-  ld_table

		snp_keep_list_checked <- keep_snp_list
	}
	
    # Sanitize the table 
    keep_table_path <- keep_table_sanitizer(snp_keep_list_checked)

    # Send the keep table to bcftools
    filtered_vcf_table <- filter_vcf_file(vcf_file_path = vcf_file_path, keep_table_path)

    # split the SNP such that they have one annotation per line
    split_table_path <- split_vcf_eff(filtered_vcf_table)

    # Send the split table to SnpSift
    # This returns a list that has a $path and -
    # a $table
    snpeff_table <- execute_snpsift(split_table_path)

    # Read the output produced by SnpSift
    snpsift_table <- snpeff_table$table

    # We start the pipeline with a tag_SNP. By the time we reach the "final_table", - 
    # we might not retain the tag_SNP if its impact factor is not HIGH or MODERATE. -
    # The potential issue of the tag_SNP being dropped can be fixed by creating an -
    # exception in the conditional here.

    # For the final table, we create a new field that indicates if a loci is a Tag - 
    # SNP or a Candidate gene. This will be recorded in TYPE.
    
    if(all.impacts){
        snpsift_table_impacts <- snpsift_table
    } else {
        snpsift_table_impacts <- snpsift_table %>% 
            filter(IMPACT %in% c("HIGH","MODERATE") | BP == bp ) # The OR condition lets us retain the tag SNP which might be dropped if the IMPACT factor is not HIGH or MODERATE
    }

    # This table should have 
    # 1. The impact factor
    # 2. The tag SNP
    # 3. The Pvalues from GWAS
    # 4. The LD with the tag_snp included
    pvalues_impact_ld_table <- snpsift_table_impacts %>%
        left_join(gwas_table, by = c("CHROM","BP")) %>%
        left_join(ld_table_checked, by = c("CHROM","BP"))

    # Annotate the tag_snp
    pvalues_impact_ld_colors_table <- pvalues_impact_ld_table %>% mutate(
        Type = case_when(
            BP == bp ~ "tag_snp",
            BP != bp ~ "Candidate"
        )
    )

    # Call the weight function
    final_reports_table <- overall_weight_func(pvalues_impact_ld_colors_table, bp = bp)

    plot <- panvar_plot(final_reports_table, nrow(gwas_table))

    return(list(plot = plot, table = final_reports_table))
}