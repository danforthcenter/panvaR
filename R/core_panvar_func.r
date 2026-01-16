#' panvar_func
#'
#' The main panvar function, runs the panvar program and optionally runs gwas if the user does not provide their own gwas output. 
#'
#' @param vcf_file_path Path to VCF file as a character string.
#' @param phenotype_data (Optional) Path (character string) to the phenotype table (TSV/CSV)  OR a data.table object. Required unless `gwas_data` is provided. The first column must contain genotype identifiers, last column should be phenotype to test.
#' @param gwas_data (Optional) Path (character string) to a pre-computed GWAS results table (TSV/CSV). If provided, `phenotype_data` and PC-related arguments are ignored. Must contain 'CHROM', 'BP', and 'Pvalues' columns.
#' @param annotation_table_path (Optional) Path to the annotation table file (TSV/CSV). Must contain 'GENE' and 'Annotation' columns. Defaults to NULL.
#' @param tag_snps SNP or SNP's as character string or vector respectively
#' @param r2_threshold LD threshold, SNP's in LD with tag SNP below this threshold are excluded
#' @param maf Minor allele frequency to filter VCF before conducting GWAS. default = .05.
#' @param missing_rate SNP's with missing rate above this threshold will be removed before GWAS. default = .1.
#' @param window Number of base-pairs to add and subtract to create a region of physical space around the tag snp. Can be supplied as an integer of base pairs or as a charcter string with the suffix "kb" denoting kilobases e.g. "500kb".
#' @param pc_min Minimum number of principal components for use in the GWAS model. Integer. 
#' @param pc_max Maximum number of principal components for use in the GWAS model. Integer. 
#' @param specific_pcs Specific prinicipal component numbers to include in the GWAS model. Numeric vector. 
#' @param dynamic_correlation Should dynamic correlation of PC's be used? If true, PC's between the minimum and maximum number of PC's are included in the GWAS model if they are siginificantly correlated with the phenotype at the \eqn{\alpha} = .001 level.
#' @param all.impacts Low impact snp's are omitted for better visualization by default. If all impacts is true, low impact snp's are included in the output.
#'
#' @returns 
#' A named list
#' - `plot`, a plot of the results of panvar.
#' - `table`, a table of the output. P-values returned are -log10(pvalues).
#'
#' @examples
#' # Using a pre-computed GWAS table
#' \dontrun{
#' # panvar_func(vcf_file_path = "<path_to_vcf_file>", 
#'               gwas_data = "<path_to_gwas_table>", 
#'               tag_snps = c("Chr_09:12456"))
#' 
#' # Using phenotype file path for de-novo GWAS
#' # panvar_func(vcf_file_path = "<path_to_vcf_file>", 
#'               phenotype_data = "<path_to_phenotype_data>", 
#'               tag_snps = c("Chr_09:12456"))
#'  }     
#'          
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
panvar_func <- function(vcf_file_path, phenotype_data = NULL, gwas_data = NULL, annotation_table_path = NULL, tag_snps = NULL, r2_threshold = 0.6, maf = 0.05, missing_rate = 0.10, window = 500000,pc_min = 5,pc_max = 5, specific_pcs = NULL,dynamic_correlation = FALSE, all.impacts = FALSE){ 
  
  # --- Start: Input Validation ---
  if (!is.null(gwas_data) && !is.null(phenotype_data)) {
    stop("Conflict: Please provide either 'phenotype_data' for de-novo GWAS or 'gwas_data' for pre-computed results, but not both.")
  }
  if (is.null(gwas_data) && is.null(phenotype_data)) {
    stop("Input needed: Please provide either 'phenotype_data' for de-novo GWAS or 'gwas_data'.")
  }

  if(!file.exists(vcf_file_path)){
    stop("The genotype file path that you provided is not accessible. Either you supplied a wrong path or the file does not exist.")
  }
  
  annotation_table <- NULL # Initialize as NULL
  if (!is.null(annotation_table_path)) {
    if (!file.exists(annotation_table_path)) {
      warning("Annotation table path provided but file not found: ", annotation_table_path, ". Proceeding without annotation.")
    } else {
      tryCatch({
        annotation_table <- data.table::fread(annotation_table_path)
        required_annot_cols <- c("GENE", "Annotation")
        if (!all(required_annot_cols %in% names(annotation_table))) {
          missing_cols <- setdiff(required_annot_cols, names(annotation_table))
          warning("Annotation table is missing required columns: ", paste(missing_cols, collapse=", "), ". It must contain 'GENE' and 'Annotation'. Proceeding without annotation.")
          annotation_table <- NULL
        } else {
          if (!is.character(annotation_table$GENE)) {
            annotation_table <- annotation_table %>% 
              mutate(GENE = as.character(.data$GENE))
            warning("Coerced 'GENE' column in annotation table to character type for joining.")
          }
          print("Annotation table loaded successfully.")
        }
      }, error = function(e) {
        warning("Error reading annotation table: ", e$message, ". Proceeding without annotation.")
        annotation_table <- NULL
      })
    }
  }

  # Check if the vcf_file has a tbi file
  proper_tbi(vcf_file_path)

  # --- GWAS Step: Either load pre-computed or run de-novo ---
  gwas_table_denovo <- NULL
  if (!is.null(gwas_data)) {
    message(">> Loading pre-computed GWAS results from: ", gwas_data)
    gwas_table_denovo <- data.table::fread(gwas_data)
  } else {
    message(">> No GWAS table provided, running de-novo GWAS analysis.")
    # Check phenotype input type
    if (is.character(phenotype_data)) {
      if (!file.exists(phenotype_data)) {
        stop(">> The phenotype file path that you provided is not accessible: ", phenotype_data)
      }
    } else if (!is(phenotype_data, "data.frame")) {
      stop(">> The phenotype_data argument must be a file path (character) or a data.frame/data.table object.")
    }
    
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
  }
  
  # convert the window into bp values
  window_bp <- window_unit_func(window)
  
  # convert the vcf file to plink format
  in_plink_format <- vcf_to_plink2(vcf_file_path)
  
  # clean up the supplied vcf file
  message("~~~~~~~~~~~~~~~ Applying QC filters to genotype file ~~~~~~~~~~~~~~~")
  cleaned_up <- bed_file_clean_up(in_plink_format$bed, maf = maf, missing_rate = missing_rate)
  
  # Validate the final GWAS table
  gwas_table <- check_gwas_table(gwas_table_denovo)
  
  # Determine tag SNPs and call convenience function
  if(is.null(tag_snps)){
    denovo_tag_snp <- tag_snp_func(gwas_table)
    print(">> Note: you did not specify a tag snp - so the tag SNP will be inferred from the GWAS results")
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
      annotation_table = annotation_table
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
      annotation_table = annotation_table
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
        annotation_table = annotation_table
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
  # Handle chromosome name differences BEFORE LD filtering
  plink2_bcf_dictionary <- plink2_bcftools_chroms_dictionary(vcf_file_path,in_plink_format$bim)
  
  # Convert chrom to PLINK format if dictionary exists
  chrom_for_plink <- chrom
  if(!is.null(plink2_bcf_dictionary)){
    # Convert VCF chromosome name to PLINK format
    chrom_match <- plink2_bcf_dictionary %>% filter(.data$vcf == chrom)
    if(nrow(chrom_match) > 0){
      chrom_for_plink <- chrom_match$plink2[1]
      message(paste0(">> Converting chromosome name from '", chrom, "' to '", chrom_for_plink, "' for PLINK compatibility"))
    }
  }
  
  # subset your genotype data around the tag snp
  subset_genotype_data <- subset_around_tag(cleaned_up,chrom = chrom_for_plink, bp = bp, window = window_bp)
  
  # using ld get the list of bps to keep
  table <- ld_filtered_snp_list(subset_genotype_data,chrom = chrom_for_plink, bp = bp, r2_threshold = r2_threshold)
  
  # Validate LD table is not empty
  if(is.null(table) || length(table) == 0 || all(is.na(table)) || (!is.data.frame(table)) || (is.data.frame(table) && nrow(table) == 0)){
    error_msg <- paste0(
      "\\n========================================\\n",
      "ERROR: No SNPs found in LD with tag SNP\\n",
      "========================================\\n",
      "Tag SNP: chr=", chrom, ", bp=", bp, "\\n",
      "PLINK format: chr=", chrom_for_plink, "\\n",
      "Current r² threshold: ", r2_threshold, "\\n",
      "Current window size: ", format(window_bp, big.mark=","), " bp\\n",
      "\\nPossible causes:\\n",
      "  1. r² threshold too high - no SNPs meet the linkage threshold\\n",
      "  2. Window size too small - try increasing to 1,000,000 bp or more\\n",
      "  3. Tag SNP not found in genotype data\\n",
      "  4. Chromosome naming mismatch between GWAS and VCF\\n",
      "\\nRecommended solutions:\\n",
      "  • LOWER r² threshold to 0.3 or 0.4 (currently: ", r2_threshold, ")\\n",
      "  • INCREASE window size to 1000000 or larger (currently: ", format(window_bp, big.mark=","), ")\\n",
      "  • VERIFY tag SNP exists in your VCF file using: bcftools view -H ", vcf_file_path, " | grep '", bp, "'\\n",
      "  • CHECK chromosome naming: GWAS uses '", chrom, "', PLINK uses '", chrom_for_plink, "'\\n",
      "========================================\\n"
    )
    stop(error_msg)
  }
  
  message(paste0(">> Found ", nrow(table), " SNPs in LD with tag SNP"))
  
  # clean up ld table a little
  ld_table <- ld_table_maker(table)
  
  # extract vector of snp names returned from ld filtering
  keep_snp_list <- snps_to_keep(table)
  
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
  
  # Debug: check keep table
  message(paste0(">> Number of SNPs to keep: ", nrow(snp_keep_list_checked)))
  
  # Filter VCF for just these snps
  filtered_vcf_table <- filter_vcf_file(vcf_file_path = vcf_file_path, keep_table_path)
  
  # Split annotations
  split_table_path <- split_vcf_eff(filtered_vcf_table)
  
  # Run SnpSift
  snpeff_table <- execute_snpsift(split_table_path)
  snpsift_table <- snpeff_table$table
  
  # Validate SnpSift output
  if(nrow(snpsift_table) == 0){
    stop("SnpSift returned an empty table. This likely means no SNPs were found in the filtered VCF file. Check chromosome naming compatibility.")
  }
  
  message(paste0(">> SnpSift extracted ", nrow(snpsift_table), " SNP annotations"))
  
  # Filter by impact
  if(all.impacts){
    snpsift_table_impacts <- snpsift_table
  } else {
    # Ensure GENE column exists before filtering (it should from snpsift)
    if (!"GENE" %in% names(snpsift_table)) {
      stop("SnpSift output table is missing the required 'GENE' column.")
    }
    snpsift_table_impacts <- snpsift_table %>%
      filter(.data$IMPACT %in% c("HIGH","MODERATE") | .data$BP == bp )
  }
  
  # Join GWAS and LD results
  pvalues_impact_ld_table <- snpsift_table_impacts %>%
    left_join(gwas_table_dicted, by = c("CHROM","BP")) %>%
    left_join(ld_table_checked, by = c("CHROM","BP"))
  
  # Annotate tag vs candidate
  pvalues_impact_ld_colors_table <- pvalues_impact_ld_table %>% mutate(
    Type = case_when(
      .data$BP == bp ~ "tag_snp",
      .data$BP != bp ~ "Candidate"
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
        final_reports_table <- final_reports_table %>% mutate(GENE = as.character(.data$GENE))
      }
      
      print("Joining with annotation table by GENE...")
      # Perform the left join using GENE
      # Select only GENE and Annotation from the annotation table to avoid duplicate columns
      final_reports_table <- final_reports_table %>%
        left_join(annotation_table %>% select(.data$GENE, .data$Annotation), by = "GENE") # <<< MODIFIED JOIN CONDITION
      
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