panvar_plot <- function(reports_table,nrows_in_gwas = NULL, pvalue_threshold = 0.05, point_size = 3, alpha_base = 0.7) {
  
  # This function very much assumes that you have supplied the default panvar results table.
  my.colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF")
  
  # This is the adjusted bonferroni pvalue
  if(is.null(nrows_in_gwas)){
    print("You did not supply a nrows for the original GWAS so making something up as a placeholder.")
    nrows_in_gwas <- 1e6
  }
  
  hline_value = -log10(
    0.05 / nrows_in_gwas
  )
  
  # What is the tag df
  
  tag_df <- reports_table %>% 
    filter(Type == "tag_snp")
  
  tag_bp <- tag_df %>%
    pull(BP)
  
  # Initialize the ggplot object with the base aesthetics and data
  plot <- ggplot(aes(x = BP, y = Pvalues), data = reports_table)
  
  # Add points to the plot with shape and fill aesthetics
  plot <- plot +
    geom_point(
      aes(shape = IMPACT, fill = LD),
      size = point_size,
      color = "black",
      alpha = alpha_base
    )
  
  # Add vertical lines from tag_df data
  plot <- plot +
    geom_vline(data = tag_df, aes(xintercept = BP), linewidth = 1.5)
  
  # Add a vertical line at gene.snp with red dashed line
  plot <- plot +
    geom_vline(xintercept = tag_bp, color = "red", linetype = "dashed", linewidth = 1.5)
  
  # Add a horizontal line at bonf.cut with black dashed line
  plot <- plot +
    geom_hline(aes(yintercept = hline_value), color = "black", linetype = "dashed")
  
  # Manually set the shape scale
  plot <- plot +
    scale_shape_manual(values = c(25, 21, 22))
  
  # Set the fill scale with custom colors and name
  plot <- plot +
    scale_fill_stepsn(colors = my.colors, name = "R2")
  
  # Apply a black and white theme
  plot <- plot +
    theme_bw()
  
  # Customize the theme elements
  plot <- plot +
    theme(
      text = element_text(size = 16),
      legend.position = "left",
      legend.text = element_text(size = 10)
    )
  
  # Customize the x-axis labels with number formatting
  plot <- plot +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
  
  # Add labels to the x and y axes
  plot <- plot +
    labs(
      x = "Position",
      y = expression(-log[10](p-value))
    )
  
  return(plot)
}
