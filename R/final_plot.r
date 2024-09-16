panvar_plot <- function(reports_table, nrows){

    # This is the adjusted bonferroni pvalue
    hline_value = -log10(
        0.05 / nrows
    )

    # A farily simple plot that the user can adjust
    panvar_plots <- reports_table %>% ggplot(aes(x = BP, y = Pvalues,color = IMPACT,alpha = LD)) +
       theme_classic() +
        geom_point() +
        geom_hline(aes(yintercept = hline_value , linetype = "Gwas Bonferroni"), color = "purple") 
    

    return(panvar_plots)

}