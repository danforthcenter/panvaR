panvar_plot <- function(reports_table, nrows){

    # This is the adjusted bonferroni pvalue
    hline_value = -log10(
        0.05 / nrows
    )

    panvar_plots <- reports_table %>% ggplot(aes(x = BP, y = Pvalues,color = IMPACT,alpha = LD)) +
       theme_classic() +
       theme(
        aspect.ratio = 0.9,
        text = element_text(face = "bold"),
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2, "cm"),
        legend.title = element_text(size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")
        ) +
        geom_point() +
        geom_hline(aes(yintercept = hline_value , linetype = "Gwas Bonferroni"), color = "purple") 
    

    return(panvar_plots)

}