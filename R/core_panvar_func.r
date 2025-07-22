#' panvar_func
#'
#' @param phenotype_data Path to the phenotype data file (character string) OR a data.table object containing phenotype data. The first column must contain genotype identifiers matching the VCF file.
#' @param vcf_file_path Path to the VCF file
#' @param annotation_table_path (Optional) Path to the annotation table file (TSV/CSV). Must contain 'GENE' and 'Annotation' columns. Defaults to NULL. # <<< MODIFIED DESCRIPTION
#' @param r2_threshold The r2 threshold Defaults to 0.6
# ... [rest of the parameters remain the same] ...
#' @param all.impacts (optional) Should all impacts be included in the report? Defaults to FALSE - in which case only "MODERATES" and "HIGH" impacts will be included
#'
#' @examples
#' # Using phenotype file path and annotation table
#' panvar_func("<path_to_phenotype_data>", "<path_to_vcf_file>", annotation_table_path = "<path_to_annotation_table>", tag_snps = c("Chr_09:12456","Chr08:14587"), r2_threshold = 0.6) # <<< EXAMPLE UPDATED
#' # Using phenotype data.table without annotation
#' pheno_dt <- data.table::fread("<path_to_phenotype_data>")
#' panvar_func(pheno_dt, "<path_to_vcf_file>", tag_snps = c("Chr_09:12456","Chr08:14587"), r2_threshold = 0.6)
#'
#' @import tidyverse
#' @import data.table
#' @import sys
#' @import parallel
#' @import bigsnpr
#' @import modelr
#' @importFrom methods is
#'
#' @export
panvar_func <- function(phenotype_data, vcf_file_path, annotation_table_path = NULL, tag_snps = NULL, r2_threshold = 0.6, maf = 0.05, missing_rate = 0.10, window = 500000,pc_min = 5,pc_max = 5, specific_pcs = NULL,dynamic_correlation = FALSE, all.impacts = FALSE){ 
  
  # --- Start: Input Validation ---
  if(!file.exists(vcf_file_path)){
    stop("The genotype file path that you provided is not accessible. Either you supplied a wrong path or the file does not exist.")
  }
  
  # Check phenotype input type
  if (is.character(phenotype_data)) {
    if (!file.exists(phenotype_data)) {
      stop("The phenotype file path that you provided is not accessible. Either you supplied a wrong path or the file does not exist: ", phenotype_data)
    }
  } else if (!is(phenotype_data, "data.frame")) {
    stop("The phenotype_data argument must be a file path (character) or a data.frame/data.table object.")
  }
  
  # <<< Start Change: Annotation Table Reading and Validation >>>
  annotation_table <- NULL # Initialize as NULL
  if (!is.null(annotation_table_path)) {
    if (!file.exists(annotation_table_path)) {
      warning("Annotation table path provided but file not found: ", annotation_table_path, ". Proceeding without annotation.")
    } else {
      tryCatch({
        annotation_table <- data.table::fread(annotation_table_path)
        # --- MODIFIED: Check for GENE and Annotation columns ---
        required_annot_cols <- c("GENE", "Annotation")
        if (!all(required_annot_cols %in% names(annotation_table))) {
          missing_cols <- setdiff(required_annot_cols, names(annotation_table))
          # --- MODIFIED: Updated warning message ---
          warning("Annotation table is missing required columns: ", paste(missing_cols, collapse=", "), ". It must contain 'GENE' and 'Annotation'. Proceeding without annotation.")
          annotation_table <- NULL # Reset to NULL if columns are missing
        } else {
          # --- REMOVED: BP numeric check is no longer needed for joining ---
          # Ensure GENE is character for joining robustness (optional, depends on data)
          if (!is.character(annotation_table$GENE)) {
            annotation_table[, GENE := as.character(GENE)]
            warning("Coerced 'GENE' column in annotation table to character type for joining.")
          }
          print("Annotation table loaded successfully.")
        }
      }, error = function(e) {
        warning("Error reading annotation table: ", e$message, ". Proceeding without annotation.")
        annotation_table <- NULL # Ensure it's NULL on error
      })
    }
  }
  # <<< End Change >>>
  
  
  # Check if the vcf_file has a tbi file
  proper_tbi(vcf_file_path)
  
  # Run GWAS for the user
  gwas_table_denovo <- panvar_gwas(
    phenotype_input = phenotype_data,
    genotype_data = vcf_file_path,
    specific_PCs = specific_pcs,
    pc_min = pc_min,
    pc_max = pc_max,
    dynamic_correlation = dynamic_correlation,
    maf = maf,
    missing_rate = missing_rate
  )
  
  # convert the window into bp values
  window_bp <- window_unit_func(window)
  
  # convert the vcf file to plink format
  in_plink_format <- vcf_to_plink2(vcf_file_path)
  
  # clean up the supplied vcf file
  cleaned_up <- bed_file_clean_up(in_plink_format$bed, maf = maf, missing_rate = missing_rate)
  
  gwas_table <- check_gwas_table(gwas_table_denovo)
  
  # Determine tag SNPs and call convenience function
  if(is.null(tag_snps)){
    denovo_tag_snp <- tag_snp_func(gwas_table_denovo)
    print("Note: you did not specify a tag snp - so the tag SNP will be inferred from the GWAS results")
    bp = denovo_tag_snp$tag_snp_bp
    chrom = denovo_tag_snp$tag_snp_chromosome
    
    panvar_result <- panvar_convienience_function(
      chrom = chrom,
      bp = bp,
      cleaned_up = cleaned_up,
      vcf_file_path = vcf_file_path,
      gwas_table = gwas_table,
      in_plink_format = in_plink_format,
      r2_threshold = r2_threshold,
      window_bp = window_bp,
      all.impacts = all.impacts,
      annotation_table = annotation_table # <<< PASS ANNOTATION TABLE
    )
    
  } else if(length(tag_snps) == 1){
    user_tag_snp <- tag_snp_splitter(tag_snps)
    chrom = user_tag_snp$chrom
    bp = user_tag_snp$bp
    
    panvar_result <- panvar_convienience_function(
      chrom = chrom,
      bp = bp,
      cleaned_up = cleaned_up,
      vcf_file_path = vcf_file_path,
      gwas_table = gwas_table,
      in_plink_format = in_plink_format,
      r2_threshold = r2_threshold,
      window_bp = window_bp,
      all.impacts = all.impacts,
      annotation_table = annotation_table # <<< PASS ANNOTATION TABLE
    )
    
  } else if(length(tag_snps) > 1){
    list_of_tag_snps <- lapply(tag_snps, tag_snp_splitter)
    
    panvar_result <- lapply(list_of_tag_snps, function(x){
      panvar_convienience_function(
        chrom = x$chrom,
        bp = x$bp,
        cleaned_up = cleaned_up,
        vcf_file_path = vcf_file_path,
        gwas_table = gwas_table,
        in_plink_format = in_plink_format,
        r2_threshold = r2_threshold,
        window_bp = window_bp,
        all.impacts = all.impacts,
        annotation_table = annotation_table # <<< PASS ANNOTATION TABLE
      )
    })
  }
  
  return(panvar_result)
}

# --- Modify panvar_convienience_function ---

panvar_convienience_function <- function(
    chrom,
    bp,
    cleaned_up,
    vcf_file_path,
    gwas_table,
    in_plink_format,
    r2_threshold = 0.6,
    window_bp = 500000,
    all.impacts = FALSE,
    annotation_table = NULL
)
{
  # subset your genotype data around the tag snp
  subset_genotype_data <- subset_around_tag(cleaned_up,chrom = chrom, bp = bp, window = window_bp)
  
  # using ld get the list of bps to keep
  table <- ld_filtered_snp_list(subset_genotype_data,chrom = chrom, bp = bp, r2_threshold = r2_threshold)
  
  # Make the LD table
  ld_table <- ld_table_maker(table)
  
  # Convert the table into the list of SNPs to keep
  keep_snp_list <- snps_to_keep(table)
  
  # Handle chromosome name differences
  plink2_bcf_dictionary <- plink2_bcftools_chroms_dictionary(vcf_file_path,in_plink_format$bim)
  
  if(!is.null(plink2_bcf_dictionary)){
    ld_table_checked <- apply_dict(plink2_bcf_dictionary, ld_table)
    snp_keep_list_checked <- apply_dict(plink2_bcf_dictionary, keep_snp_list)
    gwas_table_dicted <- apply_dict(plink2_bcf_dictionary, gwas_table)
  } else {
    ld_table_checked <-  ld_table
    snp_keep_list_checked <- keep_snp_list
    gwas_table_dicted <- gwas_table
  }
  
  # Sanitize the keep table
  keep_table_path <- keep_table_sanitizer(snp_keep_list_checked)
  
  # Filter VCF
  filtered_vcf_table <- filter_vcf_file(vcf_file_path = vcf_file_path, keep_table_path)
  
  # Split annotations
  split_table_path <- split_vcf_eff(filtered_vcf_table)
  
  # Run SnpSift
  snpeff_table <- execute_snpsift(split_table_path)
  snpsift_table <- snpeff_table$table
  
  # Filter by impact
  if(all.impacts){
    snpsift_table_impacts <- snpsift_table
  } else {
    # Ensure GENE column exists before filtering (it should from snpsift)
    if (!"GENE" %in% names(snpsift_table)) {
      stop("SnpSift output table is missing the required 'GENE' column.")
    }
    snpsift_table_impacts <- snpsift_table %>%
      filter(IMPACT %in% c("HIGH","MODERATE") | BP == bp )
  }
  
  # Join GWAS and LD results
  pvalues_impact_ld_table <- snpsift_table_impacts %>%
    left_join(gwas_table_dicted, by = c("CHROM","BP")) %>%
    left_join(ld_table_checked, by = c("CHROM","BP"))
  
  # Annotate tag vs candidate
  pvalues_impact_ld_colors_table <- pvalues_impact_ld_table %>% mutate(
    Type = case_when(
      BP == bp ~ "tag_snp",
      BP != bp ~ "Candidate"
    )
  )
  
  # Calculate weights
  final_reports_table <- overall_weight_func(pvalues_impact_ld_colors_table, bp = bp)
  
  # <<< Start Change: Optional Annotation Join by GENE >>>
  if (!is.null(annotation_table)) {
    # Ensure GENE columns have compatible types before joining
    # final_reports_table$GENE comes from snpsift
    # annotation_table$GENE was coerced to character earlier if necessary
    if ("GENE" %in% names(final_reports_table) && "GENE" %in% names(annotation_table)) {
      # Coerce final_reports_table$GENE to character just in case for safety
      if (!is.character(final_reports_table$GENE)) {
        final_reports_table <- final_reports_table %>% mutate(GENE = as.character(GENE))
      }
      
      print("Joining with annotation table by GENE...")
      # Perform the left join using GENE
      # Select only GENE and Annotation from the annotation table to avoid duplicate columns
      final_reports_table <- final_reports_table %>%
        left_join(annotation_table %>% select(GENE, Annotation), by = "GENE") # <<< MODIFIED JOIN CONDITION
      
    } else {
      # --- MODIFIED: Update warning message ---
      warning("Cannot join annotation table because 'GENE' column is missing in either the main results or the annotation table.")
    }
  }
  # <<< End Change >>>
  
  # Generate plot
  plot <- panvar_plot(final_reports_table, nrow(gwas_table))
  
  return(list(plot = plot, table = final_reports_table))
}