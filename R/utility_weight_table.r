#' Calculate Panvar Weights
#'
#' Generates weight scores for SNPs based on distance, LD, and P-values relative
#' to a specified tag SNP. Handles the edge case where only the tag SNP is present.
#'
#' @param current_table A data frame containing SNP data. Must include columns:
#'   'CHROM', 'BP' (numeric), 'Pvalues' (numeric), 'LD' (numeric).
#' @param bp The base pair position (numeric) of the tag SNP.
#'
#' @return A data frame with original columns plus calculated columns:
#'   'abs_dist', 'normalized_dist', 'normalized_LD', 'normalized_Pvalues',
#'   'final_weight'. Rows are sorted by 'final_weight' descending.
#'   If only the tag SNP exists in the input, the added columns will contain NA,
#'   and a warning will be issued.
#'
#' @importFrom dplyr %>% filter mutate select arrange bind_rows if_else row_number all_of n
#' @importFrom stats na.omit
#'
overall_weight_func <- function(current_table, bp) {
  
  # --- Input Validation (Recommended) ---
  required_cols <- c("CHROM", "BP", "Pvalues", "LD")
  if (!all(required_cols %in% names(current_table))) {
    missing_cols <- setdiff(required_cols, names(current_table))
    stop("Input table `current_table` is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(current_table$BP) ||
      !is.numeric(current_table$Pvalues) ||
      !is.numeric(current_table$LD)) {
    stop("Columns BP, Pvalues, and LD must be numeric.")
  }
  if (!is.numeric(bp) || length(bp) != 1) {
    stop("Input `bp` must be a single numeric value.")
  }
  if (nrow(current_table) == 0) {
    warning("Input `current_table` is empty. Returning an empty table with expected columns.")
    # Define expected output columns even for empty input
    expected_names <- c(names(current_table), "abs_dist", "normalized_dist",
                        "normalized_LD", "normalized_Pvalues", "final_weight")
    # Create an empty tibble with these names and appropriate types (guessing types)
    # A more robust approach might involve defining types explicitly
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
    tag_snp_data <- current_table %>%
      filter(BP == bp)
    
    # Filter out tag SNP and calculate weights for subject SNPs
    subject_snps <- current_table %>%
      filter(BP != bp | is.na(BP)) # Keep SNPs not matching BP, and also keep NAs if any
    
    # Check if subject_snps became empty after filtering
    # This covers cases where the table had only the tag SNP and potentially NAs
    if (nrow(subject_snps) == 0) {
      # This condition might seem redundant given the initial check, but serves as a safeguard
      warning("Filtered data contains no SNPs other than the tag SNP. Cannot calculate relative weights.")
      
      # Return the tag SNP data with NA weights
      weight_table <- tag_snp_data %>%
        mutate(
          abs_dist = NA_real_,
          normalized_dist = NA_real_,
          normalized_LD = NA_real_,
          normalized_Pvalues = NA_real_,
          final_weight = NA_real_
        )
      
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
            is.infinite(min_dist) ~ NA_real_, # No valid distances found
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
            is.infinite(min_ld) | is.infinite(max_ld) ~ NA_real_, # No valid LD values
            is.na(LD) ~ NA_real_, # Propagate NAs
            ld_range == 0 ~ 0.5, # All non-NA LD values are the same; assign neutral middle score? Or NA? User choice. Let's use 0.5 as a placeholder. Could also use NA_real_.
            TRUE ~ (LD - min_ld) / ld_range
          )
        )
      
      # P-value normalization (Min-Max Scaling: (x - min) / (max - min))
      # ** CRITICAL NOTE:** This formula gives HIGHER scores to HIGHER P-values.
      # Usually, LOWER p-values (more significant) are desired to have higher weight.
      # Consider INVERTING this logic if appropriate for your analysis, e.g., using:
      # (max_pval - Pvalues) / pval_range
      # OR transforming P-values first, e.g., -log10(Pvalues), and then normalizing.
      # Sticking to the original formula as requested for now.
      min_pval <- min(subject_snps$Pvalues, na.rm = TRUE)
      max_pval <- max(subject_snps$Pvalues, na.rm = TRUE)
      pval_range <- max_pval - min_pval
      subject_snps <- subject_snps %>%
        mutate(
          normalized_Pvalues = case_when(
            is.infinite(min_pval) | is.infinite(max_pval) ~ NA_real_, # No valid P-values
            is.na(Pvalues) ~ NA_real_, # Propagate NAs
            pval_range == 0 ~ 0.5, # All non-NA P-values are the same; assign neutral middle score? Or NA? Let's use 0.5. Could also use NA_real_.
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
      # This correctly adds NAs for the calculated columns in the tag_snp_data rows.
      weight_table <- bind_rows(subject_snps, tag_snp_data)
      
      # Arrange by final_weight descending. NAs are typically sorted last.
      # Using arrange(desc(final_weight)) correctly places rows with NA weights at the bottom.
      weight_table <- weight_table %>%
        arrange(desc(final_weight))
      
    } # End else block for nrow(subject_snps) > 0
  } # End main if/else block (is_only_tag or all_rows_are_tag)
  
  # Ensure consistent column order if desired (optional)
  # final_cols <- c(required_cols, "abs_dist", "normalized_dist", "normalized_LD", "normalized_Pvalues", "final_weight")
  # weight_table <- weight_table %>% select(any_of(final_cols), everything())
  
  return(weight_table)
}