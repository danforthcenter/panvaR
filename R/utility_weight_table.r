#' Calculate Panvar Weights
#'
#' Generates weight scores for SNPs based on distance, LD, and P-values relative
#' to a specified tag SNP. Handles the edge case where only the tag SNP is present.
#'
#' @param current_table A data frame containing SNP data. Must include columns:
#'   'CHROM', 'BP' (numeric), 'Pvalues' (numeric), 'LD' (numeric or coercible).
#' @param bp The base pair position (numeric) of the tag SNP.
#'
#' @return A data frame with original columns plus calculated columns:
#'   'abs_dist', 'normalized_dist', 'normalized_LD', 'normalized_Pvalues',
#'   'final_weight'. Rows are sorted by 'final_weight' descending.
#'   If only the tag SNP exists in the input, the added columns will contain NA,
#'   and a warning will be issued.
#'
#' @importFrom dplyr %>% filter mutate select arrange bind_rows if_else rowwise ungroup case_when any_of across
#' @importFrom stats na.omit setNames
#' @importFrom methods is
#'
overall_weight_func <- function(current_table, bp) {
  
  # --- Input Validation (Recommended) ---
  required_cols <- c("CHROM", "BP", "Pvalues", "LD")
  if (!all(required_cols %in% names(current_table))) {
    missing_cols <- setdiff(required_cols, names(current_table))
    stop("Input table `current_table` is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  # Coerce BP and Pvalues early, LD will be handled specifically later
  if (!is.numeric(current_table$BP)) {
    current_table <- current_table %>% mutate(BP = as.numeric(BP))
    warning("Coerced BP column to numeric.")
  }
  if (!is.numeric(current_table$Pvalues)) {
    current_table <- current_table %>% mutate(Pvalues = as.numeric(Pvalues))
    warning("Coerced Pvalues column to numeric.")
  }
  # We'll check LD numeric property *after* handling empty strings
  
  if (!is.numeric(bp) || length(bp) != 1) {
    stop("Input `bp` must be a single numeric value.")
  }
  if (nrow(current_table) == 0) {
    warning("Input `current_table` is empty. Returning an empty table with expected columns.")
    # Define expected output columns even for empty input
    expected_names <- c(names(current_table), "abs_dist", "normalized_dist",
                        "normalized_LD", "normalized_Pvalues", "final_weight")
    # Create an empty tibble with these names and appropriate types (guessing types)
    out_tbl <- dplyr::tibble(
      !!!stats::setNames(lapply(expected_names, function(x) logical()), expected_names),
      .rows = 0 # Ensure it's zero rows
    ) %>%
      # Attempt to coerce to likely types based on names
      mutate(across(any_of(c("BP", "Pvalues", "LD", "abs_dist", "normalized_dist",
                             "normalized_LD", "normalized_Pvalues", "final_weight")),
                    as.numeric))
    
    return(out_tbl)
  }
  
  # --- FIX: Ensure LD column is numeric and empty strings become NA BEFORE calculations ---
  if ("LD" %in% names(current_table)) {
    # Coerce to numeric. This will turn "" into NA with a likely warning.
    # Suppress warnings temporarily during coercion.
    current_table <- current_table %>%
      mutate(LD = suppressWarnings(as.numeric(as.character(LD))))
    
    # Now check if LD column is actually numeric after coercion attempt
    if (!is.numeric(current_table$LD)) {
      stop("LD column could not be coerced to numeric. Check for non-numeric values other than empty strings.")
    }
  } else {
    stop("Input table `current_table` is missing required column: LD") # Should have been caught earlier, but safety check
  }
  # --- END LD FIX ---
  
  
  # --- Check for Unique BP Values ---
  # Use na.omit to handle potential NAs in BP column before checking uniqueness
  unique_bps <- unique(stats::na.omit(current_table$BP))
  num_unique_bps <- length(unique_bps)
  
  # Check if the only unique non-NA BP value is the tag SNP's BP
  # Also handles the case where the input table might *only* contain NAs in BP
  is_only_tag <- num_unique_bps == 1 && !is.na(unique_bps[1]) && unique_bps[1] == bp
  
  # Also consider the case where the table contains *only* the tag SNP (even if multiple rows)
  all_rows_are_tag <- all(stats::na.omit(current_table$BP) == bp)
  
  # --- Conditional Logic: Handle Edge Case ---
  if (is_only_tag || all_rows_are_tag) {
    # Case: Only the tag SNP is present (or only NAs in BP column besides tag)
    warning("Your run was not able to find any candidates. The only value in the table is the tag SNP - can't calculate weight scores.")
    
    # Add placeholder columns with NA_real_ (numeric NA)
    # Use mutate to add columns to the original table
    weight_table <- current_table %>%
      mutate(
        abs_dist = NA_real_,         # Absolute distance is not applicable
        normalized_dist = NA_real_,  # Normalization not applicable
        normalized_LD = NA_real_,    # Normalization not applicable
        normalized_Pvalues = NA_real_,# Normalization not applicable
        final_weight = NA_real_      # Final weight cannot be calculated
      )
    # No sorting is needed as weights are NA or there's effectively one group
    
  } else {
    # --- Standard Logic: Calculate Weights ---
    
    # Separate tag SNP data (could be multiple rows if input has duplicates)
    # Ensure tag_snp_data has the same columns as subject_snps will have after calculations
    tag_snp_data <- current_table %>%
      filter(BP == bp) %>%
      mutate(
        abs_dist = NA_real_,
        normalized_dist = NA_real_,
        normalized_LD = NA_real_,
        normalized_Pvalues = NA_real_,
        final_weight = NA_real_
      )
    
    
    # Filter out tag SNP and calculate weights for subject SNPs
    subject_snps <- current_table %>%
      filter(BP != bp | is.na(BP)) # Keep SNPs not matching BP, and also keep NAs if any
    
    # Check if subject_snps became empty after filtering
    # This covers cases where the table had only the tag SNP and potentially NAs
    if (nrow(subject_snps) == 0) {
      # This condition might seem redundant given the initial check, but serves as a safeguard
      warning("Filtered data contains no SNPs other than the tag SNP. Cannot calculate relative weights.")
      
      # Return the tag SNP data (already has NA weights added)
      weight_table <- tag_snp_data
      
    } else {
      # Proceed with calculations as there are other SNPs
      
      # Calculate absolute distance from the tag SNP
      subject_snps <- subject_snps %>%
        mutate(abs_dist = abs(BP - bp)) # abs() handles NA in BP gracefully (returns NA)
      
      # --- Normalization ---
      # Note: Calculations involving min/max ignore NAs by default if na.rm=TRUE (implied in dplyr's min/max)
      # However, the normalization itself might produce NaN if the range is zero.
      
      # Distance normalization: (min_dist / abs_dist) -> higher score for closer SNPs
      min_dist <- min(subject_snps$abs_dist, na.rm = TRUE)
      subject_snps <- subject_snps %>%
        mutate(
          # Check for Inf/NaN resulting from min() on empty or all-NA data
          # Check for division by zero (abs_dist == 0)
          # Check for min_dist being Inf (if abs_dist had no non-NA values)
          normalized_dist = case_when(
            is.infinite(min_dist) | is.na(min_dist) ~ NA_real_, # No valid distances found or min() returned NA
            abs_dist == 0 ~ NA_real_, # Avoid division by zero (should not happen due to filter)
            is.na(abs_dist) ~ NA_real_, # Propagate NAs
            TRUE ~ min_dist / abs_dist # Standard case
          )
        )
      
      # LD normalization (Min-Max Scaling: (x - min) / (max - min))
      min_ld <- min(subject_snps$LD, na.rm = TRUE)
      max_ld <- max(subject_snps$LD, na.rm = TRUE)
      ld_range <- max_ld - min_ld
      subject_snps <- subject_snps %>%
        mutate(
          normalized_LD = case_when(
            is.infinite(min_ld) | is.infinite(max_ld) | is.na(min_ld) | is.na(max_ld) ~ NA_real_, # No valid LD values or min/max returned NA/Inf
            is.na(LD) ~ NA_real_, # Propagate NAs
            ld_range == 0 ~ 0.5, # All non-NA LD values are the same; assign neutral middle score. Could also use NA_real_.
            TRUE ~ (LD - min_ld) / ld_range
          )
        )
      
      # P-value normalization (Min-Max Scaling: (x - min) / (max - min))
      # This formula gives HIGHER scores to HIGHER P-values (-log10 scale)
      min_pval <- min(subject_snps$Pvalues, na.rm = TRUE)
      max_pval <- max(subject_snps$Pvalues, na.rm = TRUE)
      pval_range <- max_pval - min_pval
      subject_snps <- subject_snps %>%
        mutate(
          normalized_Pvalues = case_when(
            is.infinite(min_pval) | is.infinite(max_pval) | is.na(min_pval) | is.na(max_pval) ~ NA_real_, # No valid P-values or min/max returned NA/Inf
            is.na(Pvalues) ~ NA_real_, # Propagate NAs
            pval_range == 0 ~ 0.5, # All non-NA P-values are the same; assign neutral middle score. Could also use NA_real_.
            TRUE ~ (Pvalues - min_pval) / pval_range
          )
        )
      
      # --- Final Weight Calculation ---
      # Average the available normalized scores for each row
      # Use pmax/rowMeans/apply carefully to handle NAs
      subject_snps <- subject_snps %>%
        rowwise() %>% # Process row by row for robust NA handling in average
        mutate(
          final_weight = mean(c(normalized_dist, normalized_LD, normalized_Pvalues), na.rm = TRUE)
        ) %>%
        ungroup() %>% # Important to ungroup after rowwise operation
        # Replace NaN with NA (mean of zero non-NA values is NaN)
        mutate(final_weight = if_else(is.nan(final_weight), NA_real_, final_weight))
      
      
      # --- Combine and Sort ---
      # Use dplyr::bind_rows for robust combining. Handles column mismatches by filling with NA.
      weight_table <- bind_rows(subject_snps, tag_snp_data)
      
      # Arrange by final_weight descending. NAs are typically sorted last.
      weight_table <- weight_table %>%
        arrange(desc(final_weight))
      
    } # End else block for nrow(subject_snps) > 0
  } # End main if/else block (is_only_tag or all_rows_are_tag)
  
  
  return(weight_table)
}