panvar_plot <- function(reports_table, nrows_in_gwas = NULL, pvalue_threshold = 0.05, 
                        point_size = 3, alpha_base = 0.7, total_snps = NULL, total_genes = NULL) {
  
  # This function works with the provided molybdenum_table.tsv format
  my.colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF")
  
  # Make sure Pvalues column is already in -log10 format (which it appears to be in your data)
  # No need to transform it if it's already in that format
  
  # Set Bonferroni cutoff
  if(is.null(nrows_in_gwas)){
    message("You did not supply a nrows for the original GWAS so making something up as a placeholder.")
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
  
  # Extract tag SNPs
  tag_df <- reports_table %>% 
    dplyr::filter(Type == "tag_snp")
  
  # Handle case where there might be no tag SNPs
  if(nrow(tag_df) > 0) {
    tag_bp <- tag_df %>%
      dplyr::pull(BP)
  } else {
    tag_bp <- NULL
    warning("No tag SNPs found in the data")
  }
  
  # Initialize the ggplot object with the base aesthetics and data
  plot <- ggplot2::ggplot(ggplot2::aes(x = BP, y = Pvalues), data = reports_table)
  
  # Add points to the plot with shape and fill aesthetics
  # Handle missing LD values by using NA in the aesthetic mapping
  plot <- plot +
    ggplot2::geom_point(
      ggplot2::aes(shape = IMPACT, fill = LD),
      size = point_size,
      color = "black",
      alpha = alpha_base,
      na.rm = TRUE  # Skip missing values
    )
  
  # Add vertical lines from tag_df data if they exist
  if(nrow(tag_df) > 0) {
    plot <- plot +
      ggplot2::geom_vline(data = tag_df, ggplot2::aes(xintercept = BP), linewidth = 1.5)
    
    # Add a vertical line at tag_bp with red dashed line
    plot <- plot +
      ggplot2::geom_vline(xintercept = tag_bp, color = "red", linetype = "dashed", linewidth = 1.5)
  }
  
  # Add a horizontal line at bonf.cut with black dashed line
  plot <- plot +
    ggplot2::geom_hline(yintercept = hline_value, color = "black", linetype = "dashed")
  
  # Manually set the shape scale with fallback option
  unique_impacts <- unique(reports_table$IMPACT)
  shape_values <- rep(c(25, 21, 22), length.out = length(unique_impacts))
  names(shape_values) <- unique_impacts
  
  plot <- plot +
    ggplot2::scale_shape_manual(values = shape_values)
  
  # Set the fill scale with custom colors and name, handling NAs
  plot <- plot +
    ggplot2::scale_fill_gradient(low = "#0000FF", high = "#FF0000", name = "R2", na.value = "grey50")
  
  # Apply a black and white theme
  plot <- plot +
    ggplot2::theme_bw()
  
  # Customize the theme elements
  plot <- plot +
    ggplot2::theme(
      text = ggplot2::element_text(size = 16),
      legend.position = "left",
      legend.text = ggplot2::element_text(size = 10)
    )
  
  # Customize the x-axis labels with number formatting and better range
  # Add some padding to the x-axis to prevent points from being cut off
  x_min <- min(reports_table$BP, na.rm = TRUE)
  x_max <- max(reports_table$BP, na.rm = TRUE)
  x_padding <- (x_max - x_min) * 0.05  # 5% padding
  
  plot <- plot +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      limits = c(x_min - x_padding, x_max + x_padding)
    )
  
  # Add labels to the x and y axes, plus title with counts
  plot <- plot +
    ggplot2::labs(
      x = "Position",
      y = expression(-log[10](p-value)),
      title = paste0("Total SNPs in Window: ", total_snps, " | Total Genes in Window: ", total_genes)
    )
  
  # Return the final plot (only once!)
  return(plot)
}