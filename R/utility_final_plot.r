#' plot panvar
#' a function to make a manhattan type plot of panvar results
#'
#' @param reports_table data.frame, a results table from panvar
#' @param nrows_in_gwas number of snps from gwas, used to draw bonferroni line
#' @param pvalue_threshold alpha significance threshold for calculating bonferroni line
#' @param point_size size of points in scatter plot
#' @param alpha_base transparency of points in scatter plot
#' @param total_snps optionally provide number of snps in window, added to title
#' @param total_genes optionally provide number of genes in window, added to title
#' @param export_file filename to export, if not provided plotted to device
#' @param export_format format of export either svg (default), png or pdf
#'
#' @returns plot of panvar results
#' @export
#'
#' @examples
#' panvar_plot(panvar_results_dataframe)
panvar_plot <- function(reports_table,
                        nrows_in_gwas = NULL,
                        pvalue_threshold = 0.05, # Currently unused, but kept for potential future use
                        point_size = 3,
                        alpha_base = 0.7,
                        total_snps = NULL,
                        total_genes = NULL,
                        export_file = NULL, # New: Base filename for export (e.g., "myplot")
                        export_format = "svg" # New: "svg", "png", or "pdf"
                        ) {

  # --- Input Validation and Setup ---

  # Basic check for required columns
  required_cols <- c("BP", "Pvalues", "Type", "IMPACT", "LD", "CHROM", "GENE")
  if (!all(required_cols %in% colnames(reports_table))) {
    stop("Input 'reports_table' is missing one or more required columns: ",
         paste(setdiff(required_cols, colnames(reports_table)), collapse=", "))
  }

  # Make sure Pvalues column is numeric
  if (!is.numeric(reports_table$Pvalues)) {
     stop("'Pvalues' column must be numeric (and ideally already -log10 transformed).")
  }
   # Make sure BP column is numeric
  if (!is.numeric(reports_table$BP)) {
     stop("'BP' column must be numeric.")
  }

  # Set Bonferroni cutoff
  if(is.null(nrows_in_gwas)){
    message("nrows_in_gwas not supplied. Using 1e6 as a placeholder for Bonferroni calculation.")
    nrows_in_gwas <- 1e6
  } else if (!is.numeric(nrows_in_gwas) || nrows_in_gwas <= 0) {
     warning("Invalid nrows_in_gwas provided. Using 1e6. Please provide a positive number.")
     nrows_in_gwas <- 1e6
  }

  hline_value = -log10(0.05 / nrows_in_gwas)

  # Use provided counts or calculate if not provided
  if(is.null(total_snps)) {
    total_snps <- nrow(reports_table)
  }

  if(is.null(total_genes)) {
    total_genes <- length(unique(reports_table$GENE[!is.na(reports_table$GENE)]))
  }

  # --- Extract Tag SNP Information ---
  tag_df <- reports_table %>%
    dplyr::filter(Type == "tag_snp")

  tag_bp <- NULL
  tag_chrom <- NULL

  if(nrow(tag_df) > 0) {
    # Ensure only one unique BP value for tag snps if multiple rows exist
    unique_tag_bps <- unique(tag_df$BP)
    if (length(unique_tag_bps) > 1) {
        warning("Multiple distinct BP values found for 'tag_snp' rows. Using the first one for the vertical line.")
    }
    tag_bp <- unique_tag_bps[1] # Use the first unique BP found

    # Extract unique chromosome - guaranteed to be the same per user description
    tag_chrom <- unique(tag_df$CHROM[!is.na(tag_df$CHROM)])
    if (length(tag_chrom) == 0) {
        warning("Tag SNP(s) found, but CHROM value is missing.")
        tag_chrom <- NULL # Set back to NULL if NA was the only value
    } else if (length(tag_chrom) > 1) {
        # This case shouldn't happen based on user info, but good to handle
        warning("Multiple different CHROM values found for 'tag_snp' rows, despite expectation. Using the first one.")
        tag_chrom <- tag_chrom[1]
    }
     # If tag_chrom is still a list/vector (e.g. from multiple identical values), take the first
    if (is.list(tag_chrom) || length(tag_chrom) > 1) {
        tag_chrom <- tag_chrom[[1]]
    }

  } else {
    warning("No rows with Type == 'tag_snp' found in the data. Cannot draw tag SNP line or display chromosome.")
  }

  # --- Build Plot ---

  reports_table$IMPACT <- factor(reports_table$IMPACT, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))
  
  # Initialize the ggplot object
  plot <- ggplot2::ggplot(ggplot2::aes(x = BP, y = Pvalues), data = reports_table)

  # Add points
  # Using pch values that have separate fill and color: 21-25
  # Map IMPACT to shape, LD to fill color
  plot <- plot +
    ggplot2::geom_point(
      ggplot2::aes(shape = IMPACT, fill = LD),
      size = point_size,
      color = "black", # Outline color
      alpha = alpha_base,
      na.rm = TRUE
    )

  # Manually set the shape scale (using fillable shapes)
  shapes.to.use <- c(24, 22, 25, 23)
  names(shapes.to.use) <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  shapeScale <- scale_shape_manual(name = "Impact", values = shapes.to.use,
                                   drop = F, na.translate = F)
  
  plot <- plot +
    shapeScale

  # Set the fill scale for LD (R2)
  plot <- plot +
    viridis::scale_fill_viridis(name = bquote(R^2))

  # Add vertical line for tag SNP position (if found)
  if(!is.null(tag_bp)) {
    plot <- plot +
      ggplot2::geom_vline(xintercept = tag_bp, color = "red", linetype = "dashed", linewidth = 1.2)
  }

  # Add horizontal line for Bonferroni cutoff
  plot <- plot +
    ggplot2::geom_hline(yintercept = hline_value, color = "black", linetype = "dashed")

  # Apply theme
  plot <- plot +
    ggplot2::theme_bw()

  # Customize theme elements
  plot <- plot +
    ggplot2::theme(
      text = ggplot2::element_text(size = 14), # Slightly smaller base text
      axis.text = ggplot2::element_text(size=12),
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5), # Center title
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0), # Left-align subtitle (for chromosome)
      legend.position = "right", # Common position for legends
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  # Customize the x-axis labels and limits
  x_min <- min(reports_table$BP, na.rm = TRUE)
  x_max <- max(reports_table$BP, na.rm = TRUE)
  x_padding <- (x_max - x_min) * 0.03 # 3% padding

  plot <- plot +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      limits = c(x_min - x_padding, x_max + x_padding)
    )

  # Add labels and title/subtitle
  # Prepare subtitle for chromosome display
  chrom_subtitle <- NULL
  if (!is.null(tag_chrom)) {
      chrom_subtitle <- paste0("Chromosome: ", tag_chrom)
  } else if (nrow(tag_df) > 0 && is.null(tag_chrom)) {
      # If tag snp exists but chrom is NA
      chrom_subtitle <- "Chromosome: (Not Available for Tag SNP)"
  } # Else: No tag SNP, no subtitle needed.

  plot <- plot +
    ggplot2::labs(
      x = "Position (bp)",
      #y = expression(-log[10](italic(P)-value)),
      y = "-log(p-value)",
      title = paste0("Total SNPs: ", format(total_snps, big.mark=","),
                     " | Total Genes: ", format(total_genes, big.mark=",")),
      subtitle = chrom_subtitle # Use the subtitle for chromosome info
    )

  # --- Export Plot (Optional) ---
  if (!is.null(export_file)) {
    # Ensure export format is lowercase
    export_format <- tolower(export_format)

    # Validate format
    valid_formats <- c("svg", "png", "pdf")
    if (!export_format %in% valid_formats) {
      warning("Invalid export_format specified ('", export_format,
              "'). Must be one of: ", paste(valid_formats, collapse=", "),
              ". Defaulting to 'svg'.")
      export_format <- "svg"
    }

    # Construct filename with extension
    filename <- paste0(export_file, ".", export_format)

    message("Attempting to save plot to: ", filename)

    if (export_format == "svg") {
      # Try SVG export, fallback to PNG on error
      svg_success <- tryCatch({
        # Check if svglite is available *before* trying to use it
        if (!requireNamespace("svglite", quietly = TRUE)) {
          stop("'svglite' package needed for SVG export. Please install it.", call. = FALSE)
        }
        ggplot2::ggsave(filename = filename, plot = plot, device = "svg", width = 10, height = 6)
        TRUE # Return TRUE on success
      }, error = function(e) {
        warning("SVG export failed. Error: ", e$message,
                "\nAttempting to save as PNG instead. (Is 'svglite' package installed and working?)",
                call. = FALSE)
        # Construct PNG filename as fallback
        fallback_filename <- sub("\\.svg$", ".png", filename)
        if (!endsWith(fallback_filename, ".png")){ # Ensure it ends with .png
             fallback_filename <- paste0(export_file, ".png")
        }
        tryCatch({
             ggplot2::ggsave(filename = fallback_filename, plot = plot, device = "png", width = 10, height = 6, dpi = 300)
             message("Successfully saved plot as PNG (fallback): ", fallback_filename)
             # We could return FALSE here to indicate SVG failed, but the message suffices.
             FALSE # Indicate original SVG save failed
        }, error = function(e2){
             warning("Fallback PNG export also failed. Error: ", e2$message, call. = FALSE)
             FALSE # Indicate fallback PNG save also failed
        })
      })
      if (svg_success) {
        message("Successfully saved plot as SVG: ", filename)
      }

    } else if (export_format == "png") {
      # Export as PNG
       tryCatch({
            ggplot2::ggsave(filename = filename, plot = plot, device = "png", width = 10, height = 6, dpi = 300)
            message("Successfully saved plot as PNG: ", filename)
        }, error = function(e){
             warning("PNG export failed. Error: ", e$message, call. = FALSE)
        })

    } else if (export_format == "pdf") {
      # Export as PDF
       tryCatch({
            # Ensure device is explicitly set for potentially missing fonts etc.
            ggplot2::ggsave(filename = filename, plot = plot, device = cairo_pdf, width = 10, height = 6)
            message("Successfully saved plot as PDF: ", filename)
        }, error = function(e){
             warning("PDF export failed. Error: ", e$message,
                     "\nTrying default pdf device (may have font issues)...", call. = FALSE)
             tryCatch({
                  ggplot2::ggsave(filename = filename, plot = plot, device = "pdf", width = 10, height = 6)
                  message("Successfully saved plot as PDF (using default device): ", filename)
             }, error = function(e2){
                  warning("Default PDF export also failed. Error: ", e2$message, call. = FALSE)
             })
        })
    }
  } # End of export block

  # --- Return Plot Object ---
  return(plot)
}
