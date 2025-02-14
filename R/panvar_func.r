#' panvar_func
#' 
#' @param phenotype_data_path Path to the phenotype data
#' 
#' @param vcf_file_path Path to the VCF file
#' 
#' @param r2_threshold The r2 threshold
#' 
#' @param tag_snps The tag SNPs that you want to pass through panvar - either a single SNP or a list of SNPs.
#' A single SNP should be formatted as "<Chr>:<BP>" - For example, "Chr_09:12456" or "Chr08:14587".
#' To supply multiple tag SNPs supply a vector of tag SNPs. For example tag_snps = c("Chr_09:12456","Chr08:14587").
#' 
#' @param r2_threshold The r2 threshold
#' Defaults to 0.6
#' 
#' @param maf The minor Allele Frequency
#' Defaults to 0.05
#' 
#' @param missing_rate The missing rate filter for your genotype data
#' Defaults to 0.1
#' 
#' @param window PanvaR determines the Linkage Disequilibrium (LD) between the tag Single Nucleotide Polymorphism (SNP) and every other SNP within a genome segment. This segment is centered on the tag SNP, extending up to a specified window size in both directions.
#' Defaults to 500000
#' 
#' @param pc_min (optional) What is the minimum number of PCs that should be included in GWAS?
#' Defaults to 5
#' 
#' @param pc_max (optional) What is the maximum number of PCs that should be included in GWAS?
#' Defaults to 5
#' 
#' @param dynamic_correlation (optional) Should the PCs, beyond minimum, be calculated dynamically?
#' 
#' @param all_impacts (optional) Should all impacts be included in the report?
#' Defaults to FALSE - in which case only "MODERATES" and "HIGH" impacts will be included
#' 
#' @examples
#' panvar("<path_to_phenotype_data>", "<path_to_vcf_file>", tag_snps = c("Chr_09:12456","Chr08:14587"), r2_threshold = 0.6)
#' 
#' @import tidyverse
#' @import data.table
#' @import sys
#' @import parallel
#' @import bigsnpr
#' @import modelr
#'
#' @export
panvar_func <- function(phenotype_data_path,vcf_file_path,tag_snps = NULL, r2_threshold = 0.6, maf = 0.05, missing_rate = 0.10, window = 500000,pc_min = 5,pc_max = 5, dynamic_correlation = FALSE, all.impacts = FALSE){

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

    # convert the window into bp values
    window_bp <- window_unit_func(window)

    # The end goal of this function is to convieneintly make
    # 1. The plot from Panvar
    # 2. The report table -
    # The report table should have - Gwas pvalues,

    # convert the vcf file to plink format
    in_plink_format <- vcf_to_plink2(vcf_file_path)

    # clean up the supplied vcf file
    cleaned_up <- bed_file_clean_up(in_plink_format$bed, maf = maf, missing_rate = missing_rate)

    gwas_table <- check_gwas_table(gwas_table_denovo)
    # if the user has not supplied a tag SNP -
    # get the tag snp from the gwas results
    if(is.null(tag_snps)){

        denovo_tag_snp <- tag_snp_func(gwas_table_denovo)

        print("Note: you did not specify a tag snp - so the tag SNP will be inferred from the GWAS results")

        bp = denovo_tag_snp$tag_snp_bp
        
        chrom = denovo_tag_snp$tag_snp_chromosome

        panvar_result <- panvar_convienience_function( # Visit function body for documentation.
            chrom = chrom,
            bp = bp,
            cleaned_up = cleaned_up,
            vcf_file_path = vcf_file_path,
            gwas_table = gwas_table,
            in_plink_format = in_plink_format,
            r2_threshold = r2_threshold,
            window_bp = window_bp,
            all.impacts = all.impacts
        )

    } else if(length(tag_snps) == 1){
        user_tag_snp <- tag_snp_splitter(tag_snps)

        chrom = user_tag_snp$chrom

        bp = user_tag_snp$bp

        panvar_result <- panvar_convienience_function( # Visit function body for documentation.
            chrom = chrom,
            bp = bp,
            cleaned_up = cleaned_up,
            vcf_file_path = vcf_file_path,
            gwas_table = gwas_table,
            in_plink_format = in_plink_format,
            r2_threshold = r2_threshold,
            window_bp = window_bp,
            all.impacts = all.impacts
        )

    } else if(length(tag_snps) > 1){

        list_of_tag_snps <- lapply(tag_snps, tag_snp_splitter)

        panvar_result <- lapply(list_of_tag_snps, function(x){ # To split the list_of_tag_snps into their constituet parts.
            panvar_convienience_function(# Visit function body for documentation.
                chrom = x$chrom,
                bp = x$bp,
                cleaned_up = cleaned_up,
                vcf_file_path = vcf_file_path,
                gwas_table = gwas_table,
                in_plink_format = in_plink_format,
                r2_threshold = r2_threshold,
                window_bp = window_bp,
                all.impacts = all.impacts
            )
        })
    }

    return(panvar_result)
}

# panvar_convenience_function is an automation tool for the panvar logic, - 
# following input verification and atomization. For instance, it can convert - 
# a phenotype table into a GWAS table. Atomization involves splitting bundled - 
# parameters, such as multiple tag SNPs, for individual processing.

panvar_convienience_function <- function(
    chrom,
    bp,
    cleaned_up,
    vcf_file_path,
    gwas_table,
    in_plink_format,
    r2_threshold = 0.6, 
    window_bp = 500000,
    all.impacts = FALSE 
)
{
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

	# If the dictionary is not empty then use it.
	if(!is.null(plink2_bcf_dictionary)){
		ld_table_checked <- apply_dict(plink2_bcf_dictionary, ld_table)

		snp_keep_list_checked <- apply_dict(plink2_bcf_dictionary, keep_snp_list)

        gwas_table_dicted <- apply_dict(plink2_bcf_dictionary, gwas_table)
	} else{
		
		ld_table_checked <-  ld_table

		snp_keep_list_checked <- keep_snp_list

        gwas_table_dicted <- gwas_table
	}
	
    # Sanitize the table 
	# TODO:- Which table?
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
        left_join(gwas_table_dicted, by = c("CHROM","BP")) %>%
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