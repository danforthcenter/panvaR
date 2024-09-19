# Currently mostly conceptual

# The goal of this function is to make a weight value
# for the Panvar analysis.
# The table should have - CHROM, BP, Pvalues, LD.
# This function will generate weight scores -
# using normalized and then averaged scores for LD, Pvalues and distance.

overall_weight_func <- function(current_table, bp){

    subject_snps <- current_table %>% 
        filter(BP != bp) %>%
        mutate(
            abs_dist = abs(BP - bp)
        ) %>%
        mutate(
            normalized_dist = ( min(abs_dist) / abs_dist)
        ) %>%
        mutate(
            normalized_LD = (LD - min(LD)) / (max(LD) - min(LD))
        ) %>%
        mutate(
            normalized_Pvalues = (Pvalues - min(Pvalues)) / (max(Pvalues) - min(Pvalues))
        )%>%
        mutate(
            final_weight = ((normalized_dist + normalized_LD + normalized_Pvalues) / 3)
        )
    
    tag_snp_data <- current_table %>%
        filter(BP == bp)

    # rbind data with the fill set to true
    weight_table <- rbind(subject_snps, tag_snp_data, fill = TRUE)

    weight_table <- weight_table %>%
        arrange(desc(final_weight))

    return(weight_table)
}

