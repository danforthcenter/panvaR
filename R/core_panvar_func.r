#' panvar_func
#'
#' @param phenotype_data Path to the phenotype data file (character string) OR a data.table object containing phenotype data. The first column must contain genotype identifiers matching the VCF file.
#' @param vcf_file_path Path to the VCF file
#' @param annotation_table_path (Optional) Path to the annotation table file (TSV/CSV). Must contain 'GENE' and 'Annotation' columns. Defaults to NULL.
#' @param r2_threshold The r2 threshold Defaults to 0.6
# ... [rest of the parameters remain the same] ...
#' @param all.impacts (optional) Should all impacts be included in the report? Defaults to FALSE - in which case only "MODERATES" and "HIGH" impacts will be included
#' @param auto_generate_tbi (Optional) If TRUE, automatically generates the .tbi index file for the VCF if it is missing. This is useful for non-interactive scripts. Defaults to FALSE.
#'
#' @examples
#' # Using phenotype file path and annotation table
#' # panvar_func("<path_to_phenotype_data>", "<path_to_vcf_file>", annotation_table_path = "<path_to_annotation_table>", tag_snps = c("Chr_09:12456","Chr08:14587"), r2_threshold = 0.6, auto_generate_tbi = TRUE)
#' # Using phenotype data.table without annotation
#' # pheno_dt <- data.table::fread("<path_to_phenotype_data>")
#' # panvar_func(pheno_dt, "<path_to_vcf_file>", tag_snps = c("Chr_09:12456","Chr08:14587"), r2_threshold = 0.6)
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
panvar_func <- function(phenotype_data, vcf_file_path, annotation_table_path = NULL, tag_snps = NULL, r2_threshold = 0.6, maf = 0.05, missing_rate = 0.10, window = 500000,pc_min = 5,pc_max = 5, specific_pcs = NULL,dynamic_correlation = FALSE, all.impacts = FALSE, auto_generate_tbi = FALSE){

  # --- Start: Profiling Initialization ---
  run_timestamp_start_overall <- Sys.time()
  profiling_records <- list()
  cores_available <- good_core_count() # Assuming good_core_count() is defined elsewhere and works

  # Helper to get memory stats from gc()
  get_gc_stats_df <- function() {
    gc_res <- gc(verbose = FALSE) # gc(reset = TRUE) could be used for max_used stats relative to this point
    df <- as.data.frame(gc_res)
    # Ensure column names are consistent for older/newer R versions if they differ
    # For simplicity, assuming "used (Mb)" and "max used (Mb)" exist or similar.
    # If not, this part might need adjustment based on exact gc() output on target R.
    # Fallback if specific columns are missing
    ncells_used_mb <- if ("used (Mb)" %in% colnames(df)) df["Ncells", "used (Mb)"] else NA
    vcells_used_mb <- if ("used (Mb)" %in% colnames(df)) df["Vcells", "used (Mb)"] else NA
    ncells_max_used_mb <- if ("max used (Mb)" %in% colnames(df)) df["Ncells", "max used (Mb)"] else NA
    vcells_max_used_mb <- if ("max used (Mb)" %in% colnames(df)) df["Vcells", "max used (Mb)"] else NA
    
    return(data.frame(
      Ncells_used_Mb = ncells_used_mb,
      Vcells_used_Mb = vcells_used_mb,
      Ncells_max_used_Mb = ncells_max_used_mb,
      Vcells_max_used_Mb = vcells_max_used_mb
    ))
  }
  
  initial_memory_stats_df <- get_gc_stats_df()
  
  profiling_records[[length(profiling_records) + 1]] <- list(
    event_name = "panvar_func_start",
    timestamp = as.character(run_timestamp_start_overall),
    cores_available = cores_available,
    Ncells_used_Mb = initial_memory_stats_df$Ncells_used_Mb,
    Vcells_used_Mb = initial_memory_stats_df$Vcells_used_Mb,
    Ncells_max_used_Mb = initial_memory_stats_df$Ncells_max_used_Mb,
    Vcells_max_used_Mb = initial_memory_stats_df$Vcells_max_used_Mb
  )
  # --- End: Profiling Initialization ---

  # --- Start: Input Validation ---
  if(!file.exists(vcf_file_path)){
    stop("The genotype file path that you provided is not accessible. Either you supplied a wrong path or the file does not exist.")
  }

  if (is.character(phenotype_data)) {
    if (!file.exists(phenotype_data)) {
      stop("The phenotype file path that you provided is not accessible. Either you supplied a wrong path or the file does not exist: ", phenotype_data)
    }
  } else if (!is(phenotype_data, "data.frame")) {
    stop("The phenotype_data argument must be a file path (character) or a data.frame/data.table object.")
  }

  annotation_table <- NULL
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
            annotation_table[, GENE := as.character(GENE)]
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

  # --- Profile: proper_tbi ---
  pt_start_time <- Sys.time()
  proper_tbi(vcf_file_path, auto_generate_tbi = auto_generate_tbi)
  pt_end_time <- Sys.time()
  profiling_records[[length(profiling_records) + 1]] <- list(
    event_name = "proper_tbi",
    start_time = as.character(pt_start_time),
    end_time = as.character(pt_end_time),
    duration_seconds = as.numeric(difftime(pt_end_time, pt_start_time, units = "secs")),
    cores_available = cores_available
  )
  # --- End Profile: proper_tbi ---
  
  # --- Profile: panvar_gwas ---
  gwas_start_time <- Sys.time()
  gwas_table_denovo <- panvar_gwas(
    phenotype_input = phenotype_data,
    genotype_data = vcf_file_path,
    specific_PCs = specific_pcs,
    pc_min = pc_min,
    pc_max = pc_max,
    dynamic_correlation = dynamic_correlation,
    maf = maf,
    missing_rate = missing_rate,
    auto_generate_tbi = auto_generate_tbi
  )
  gwas_end_time <- Sys.time()
  profiling_records[[length(profiling_records) + 1]] <- list(
    event_name = "panvar_gwas",
    start_time = as.character(gwas_start_time),
    end_time = as.character(gwas_end_time),
    duration_seconds = as.numeric(difftime(gwas_end_time, gwas_start_time, units = "secs")),
    cores_available = cores_available
  )
  # --- End Profile: panvar_gwas ---

  window_bp <- window_unit_func(window)

  # --- Profile: vcf_to_plink2 ---
  vcf2plink_start_time <- Sys.time()
  in_plink_format <- vcf_to_plink2(vcf_file_path, auto_generate_tbi = auto_generate_tbi)
  vcf2plink_end_time <- Sys.time()
  profiling_records[[length(profiling_records) + 1]] <- list(
    event_name = "vcf_to_plink2",
    start_time = as.character(vcf2plink_start_time),
    end_time = as.character(vcf2plink_end_time),
    duration_seconds = as.numeric(difftime(vcf2plink_end_time, vcf2plink_start_time, units = "secs")),
    cores_available = cores_available
  )
  # --- End Profile: vcf_to_plink2 ---

  # --- Profile: bed_file_clean_up ---
  bedclean_start_time <- Sys.time()
  cleaned_up <- bed_file_clean_up(in_plink_format$bed, maf = maf, missing_rate = missing_rate)
  bedclean_end_time <- Sys.time()
  profiling_records[[length(profiling_records) + 1]] <- list(
    event_name = "bed_file_clean_up",
    start_time = as.character(bedclean_start_time),
    end_time = as.character(bedclean_end_time),
    duration_seconds = as.numeric(difftime(bedclean_end_time, bedclean_start_time, units = "secs")),
    cores_available = cores_available
  )
  # --- End Profile: bed_file_clean_up ---

  gwas_table <- check_gwas_table(gwas_table_denovo)

  panvar_result <- NULL # Initialize

  if(is.null(tag_snps)){
    denovo_tag_snp <- tag_snp_func(gwas_table_denovo)
    print("Note: you did not specify a tag snp - so the tag SNP will be inferred from the GWAS results")
    bp = denovo_tag_snp$tag_snp_bp
    chrom = denovo_tag_snp$tag_snp_chromosome
    
    # --- Profile: panvar_convienience_function (single call) ---
    conv_start_time <- Sys.time()
    convenience_output <- panvar_convienience_function(
      chrom = chrom, bp = bp, cleaned_up = cleaned_up, vcf_file_path = vcf_file_path,
      gwas_table = gwas_table, in_plink_format = in_plink_format, r2_threshold = r2_threshold,
      window_bp = window_bp, all.impacts = all.impacts, annotation_table = annotation_table,
      cores_available_for_profiling = cores_available, # Pass cores for internal records
      auto_generate_tbi = auto_generate_tbi
    )
    conv_end_time <- Sys.time()
    
    panvar_result <- list(plot = convenience_output$plot, table = convenience_output$table)
    
    profiling_records[[length(profiling_records) + 1]] <- list(
      event_name = paste0("panvar_convienience_function_call_for_tag_", chrom, "_", bp),
      start_time = as.character(conv_start_time),
      end_time = as.character(conv_end_time),
      duration_seconds = as.numeric(difftime(conv_end_time, conv_start_time, units = "secs")),
      cores_available = cores_available
    )
    if (length(convenience_output$profiling) > 0) {
       profiling_records <- c(profiling_records, convenience_output$profiling)
    }
    # --- End Profile ---
    
  } else if(length(tag_snps) == 1){
    user_tag_snp <- tag_snp_splitter(tag_snps)
    chrom = user_tag_snp$chrom
    bp = user_tag_snp$bp
    
    # --- Profile: panvar_convienience_function (single call) ---
    conv_start_time <- Sys.time()
    convenience_output <- panvar_convienience_function(
      chrom = chrom, bp = bp, cleaned_up = cleaned_up, vcf_file_path = vcf_file_path,
      gwas_table = gwas_table, in_plink_format = in_plink_format, r2_threshold = r2_threshold,
      window_bp = window_bp, all.impacts = all.impacts, annotation_table = annotation_table,
      cores_available_for_profiling = cores_available,
      auto_generate_tbi = auto_generate_tbi
    )
    conv_end_time <- Sys.time()

    panvar_result <- list(plot = convenience_output$plot, table = convenience_output$table)

    profiling_records[[length(profiling_records) + 1]] <- list(
      event_name = paste0("panvar_convienience_function_call_for_tag_", chrom, "_", bp),
      start_time = as.character(conv_start_time),
      end_time = as.character(conv_end_time),
      duration_seconds = as.numeric(difftime(conv_end_time, conv_start_time, units = "secs")),
      cores_available = cores_available
    )
    if (length(convenience_output$profiling) > 0) {
       profiling_records <- c(profiling_records, convenience_output$profiling)
    }
    # --- End Profile ---
    
  } else if(length(tag_snps) > 1){
    list_of_tag_snps <- lapply(tag_snps, tag_snp_splitter)
    
    # --- Profile: lapply block for panvar_convienience_function ---
    lapply_start_time <- Sys.time()
    panvar_result_components <- lapply(list_of_tag_snps, function(x){
      iter_start_time <- Sys.time()
      convenience_output <- panvar_convienience_function(
        chrom = x$chrom, bp = x$bp, cleaned_up = cleaned_up, vcf_file_path = vcf_file_path,
        gwas_table = gwas_table, in_plink_format = in_plink_format, r2_threshold = r2_threshold,
        window_bp = window_bp, all.impacts = all.impacts, annotation_table = annotation_table,
        cores_available_for_profiling = cores_available,
        auto_generate_tbi = auto_generate_tbi
      )
      iter_end_time <- Sys.time()
      
      outer_call_profile_record <- list(
        event_name = paste0("panvar_convienience_function_call_for_tag_", x$chrom, "_", x$bp),
        start_time = as.character(iter_start_time),
        end_time = as.character(iter_end_time),
        duration_seconds = as.numeric(difftime(iter_end_time, iter_start_time, units = "secs")),
        cores_available = cores_available
      )
      
      list(
        main_output = list(plot = convenience_output$plot, table = convenience_output$table),
        internal_profiling = convenience_output$profiling,
        outer_call_profiling = outer_call_profile_record
      )
    })
    lapply_end_time <- Sys.time()
    
    profiling_records[[length(profiling_records) + 1]] <- list(
      event_name = "lapply_panvar_convienience_function_multiple_tags",
      start_time = as.character(lapply_start_time),
      end_time = as.character(lapply_end_time),
      duration_seconds = as.numeric(difftime(lapply_end_time, lapply_start_time, units = "secs")),
      cores_available = cores_available
    )
    
    panvar_result <- lapply(panvar_result_components, `[[`, "main_output")
    
    for (comp in panvar_result_components) {
      profiling_records[[length(profiling_records) + 1]] <- comp$outer_call_profiling
      if (length(comp$internal_profiling) > 0) {
        profiling_records <- c(profiling_records, comp$internal_profiling)
      }
    }
    # --- End Profile ---
  }
  
  # --- Start: Profiling Finalization ---
  run_timestamp_end_overall <- Sys.time()
  total_duration_sec <- as.numeric(difftime(run_timestamp_end_overall, run_timestamp_start_overall, units = "secs"))
  
  final_memory_stats_df <- get_gc_stats_df()

  profiling_records[[length(profiling_records) + 1]] <- list(
    event_name = "panvar_func_total_execution",
    start_time = as.character(run_timestamp_start_overall),
    end_time = as.character(run_timestamp_end_overall),
    duration_seconds = total_duration_sec,
    cores_available = cores_available,
    initial_Ncells_used_Mb = initial_memory_stats_df$Ncells_used_Mb,
    initial_Vcells_used_Mb = initial_memory_stats_df$Vcells_used_Mb,
    initial_Ncells_max_used_Mb = initial_memory_stats_df$Ncells_max_used_Mb,
    initial_Vcells_max_used_Mb = initial_memory_stats_df$Vcells_max_used_Mb,
    final_Ncells_used_Mb = final_memory_stats_df$Ncells_used_Mb,
    final_Vcells_used_Mb = final_memory_stats_df$Vcells_used_Mb,
    final_Ncells_max_used_Mb = final_memory_stats_df$Ncells_max_used_Mb,
    final_Vcells_max_used_Mb = final_memory_stats_df$Vcells_max_used_Mb
  )

  if (length(profiling_records) > 0) {
    # Ensure all list elements are named lists for rbindlist
    valid_records <- lapply(profiling_records, function(x) {
      if(is.list(x) && !is.null(names(x))) return(x)
      return(NULL) # Or handle error/warning
    })
    valid_records <- Filter(Negate(is.null), valid_records)

    if(length(valid_records) > 0) {
        profiling_dt <- data.table::rbindlist(valid_records, fill = TRUE, use.names = TRUE)
        profile_filename <- paste0(format(run_timestamp_start_overall, "%Y%m%d_%H%M%S"), "_panvar_run_profile.tsv")
        tryCatch({
          data.table::fwrite(profiling_dt, file = profile_filename, sep = "\t", na = "NA", quote = FALSE)
          message("Profiling data written to: ", profile_filename)
        }, error = function(e) {
          warning("Failed to write profiling data: ", e$message)
        })
    } else {
        warning("No valid profiling records to write.")
    }
  }
  # --- End: Profiling Finalization ---
  
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
    annotation_table = NULL,
    cores_available_for_profiling = 1, # New argument for profiling consistency
    auto_generate_tbi = FALSE
)
{
  # --- Profiling Init for convenience function ---
  internal_profiling_records <- list()
  
  # --- Profile: subset_around_tag ---
  sat_start_time <- Sys.time()
  subset_genotype_data <- subset_around_tag(cleaned_up,chrom = chrom, bp = bp, window = window_bp)
  sat_end_time <- Sys.time()
  internal_profiling_records[[length(internal_profiling_records) + 1]] <- list(
    event_name = "subset_around_tag",
    start_time = as.character(sat_start_time),
    end_time = as.character(sat_end_time),
    duration_seconds = as.numeric(difftime(sat_end_time, sat_start_time, units = "secs")),
    cores_available = cores_available_for_profiling
  )
  # --- End Profile ---
  
  # --- Profile: ld_filtered_snp_list ---
  ldfilt_start_time <- Sys.time()
  table <- ld_filtered_snp_list(subset_genotype_data,chrom = chrom, bp = bp, r2_threshold = r2_threshold)
  ldfilt_end_time <- Sys.time()
  internal_profiling_records[[length(internal_profiling_records) + 1]] <- list(
    event_name = "ld_filtered_snp_list",
    start_time = as.character(ldfilt_start_time),
    end_time = as.character(ldfilt_end_time),
    duration_seconds = as.numeric(difftime(ldfilt_end_time, ldfilt_start_time, units = "secs")),
    cores_available = cores_available_for_profiling
  )
  # --- End Profile ---
  
  ld_table <- ld_table_maker(table)
  keep_snp_list <- snps_to_keep(table)
  
  plink2_bcf_dictionary <- plink2_bcftools_chroms_dictionary(vcf_file_path,in_plink_format$bim, auto_generate_tbi = auto_generate_tbi)
  
  if(!is.null(plink2_bcf_dictionary)){
    ld_table_checked <- apply_dict(plink2_bcf_dictionary, ld_table)
    snp_keep_list_checked <- apply_dict(plink2_bcf_dictionary, keep_snp_list)
    gwas_table_dicted <- apply_dict(plink2_bcf_dictionary, gwas_table)
  } else {
    ld_table_checked <-  ld_table
    snp_keep_list_checked <- keep_snp_list
    gwas_table_dicted <- gwas_table
  }
  
  keep_table_path <- keep_table_sanitizer(snp_keep_list_checked)
  
  # --- Profile: filter_vcf_file ---
  fvcf_start_time <- Sys.time()
  filtered_vcf_table <- filter_vcf_file(vcf_file_path = vcf_file_path, keep_table_path, auto_generate_tbi = auto_generate_tbi)
  fvcf_end_time <- Sys.time()
  internal_profiling_records[[length(internal_profiling_records) + 1]] <- list(
    event_name = "filter_vcf_file",
    start_time = as.character(fvcf_start_time),
    end_time = as.character(fvcf_end_time),
    duration_seconds = as.numeric(difftime(fvcf_end_time, fvcf_start_time, units = "secs")),
    cores_available = cores_available_for_profiling
  )
  # --- End Profile ---
  
  # --- Profile: split_vcf_eff ---
  sve_start_time <- Sys.time()
  split_table_path <- split_vcf_eff(filtered_vcf_table)
  sve_end_time <- Sys.time()
  internal_profiling_records[[length(internal_profiling_records) + 1]] <- list(
    event_name = "split_vcf_eff",
    start_time = as.character(sve_start_time),
    end_time = as.character(sve_end_time),
    duration_seconds = as.numeric(difftime(sve_end_time, sve_start_time, units = "secs")),
    cores_available = cores_available_for_profiling
  )
  # --- End Profile ---
  
  # --- Profile: execute_snpsift ---
  ess_start_time <- Sys.time()
  snpeff_table <- execute_snpsift(split_table_path)
  ess_end_time <- Sys.time()
  internal_profiling_records[[length(internal_profiling_records) + 1]] <- list(
    event_name = "execute_snpsift",
    start_time = as.character(ess_start_time),
    end_time = as.character(ess_end_time),
    duration_seconds = as.numeric(difftime(ess_end_time, ess_start_time, units = "secs")),
    cores_available = cores_available_for_profiling
  )
  # --- End Profile ---
  snpsift_table <- snpeff_table$table
  
  if(all.impacts){
    snpsift_table_impacts <- snpsift_table
  } else {
    if (!"GENE" %in% names(snpsift_table)) {
      stop("SnpSift output table is missing the required 'GENE' column.")
    }
    snpsift_table_impacts <- snpsift_table %>%
      filter(IMPACT %in% c("HIGH","MODERATE") | BP == bp )
  }
  
  pvalues_impact_ld_table <- snpsift_table_impacts %>%
    left_join(gwas_table_dicted, by = c("CHROM","BP")) %>%
    left_join(ld_table_checked, by = c("CHROM","BP"))
  
  pvalues_impact_ld_colors_table <- pvalues_impact_ld_table %>% mutate(
    Type = case_when(
      BP == bp ~ "tag_snp",
      BP != bp ~ "Candidate"
    )
  )
  
  # --- Profile: overall_weight_func ---
  owf_start_time <- Sys.time()
  final_reports_table <- overall_weight_func(pvalues_impact_ld_colors_table, bp = bp)
  owf_end_time <- Sys.time()
  internal_profiling_records[[length(internal_profiling_records) + 1]] <- list(
    event_name = "overall_weight_func",
    start_time = as.character(owf_start_time),
    end_time = as.character(owf_end_time),
    duration_seconds = as.numeric(difftime(owf_end_time, owf_start_time, units = "secs")),
    cores_available = cores_available_for_profiling
  )
  # --- End Profile ---
  
  if (!is.null(annotation_table)) {
    if ("GENE" %in% names(final_reports_table) && "GENE" %in% names(annotation_table)) {
      if (!is.character(final_reports_table$GENE)) {
        final_reports_table <- final_reports_table %>% mutate(GENE = as.character(GENE))
      }
      print("Joining with annotation table by GENE...")
      final_reports_table <- final_reports_table %>%
        left_join(annotation_table %>% select(GENE, Annotation), by = "GENE")
    } else {
      warning("Cannot join annotation table because 'GENE' column is missing in either the main results or the annotation table.")
    }
  }
  
  # --- Profile: panvar_plot (can be skipped if plotting is not a performance concern) ---
  pp_start_time <- Sys.time()
  plot <- panvar_plot(final_reports_table, nrow(gwas_table))
  pp_end_time <- Sys.time()
  internal_profiling_records[[length(internal_profiling_records) + 1]] <- list(
    event_name = "panvar_plot",
    start_time = as.character(pp_start_time),
    end_time = as.character(pp_end_time),
    duration_seconds = as.numeric(difftime(pp_end_time, pp_start_time, units = "secs")),
    cores_available = cores_available_for_profiling
  )
  # --- End Profile ---
  
  # Return main results and profiling data
  return(list(plot = plot, table = final_reports_table, profiling = internal_profiling_records))
}