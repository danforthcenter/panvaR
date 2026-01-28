#' Extracts snpeff annotations from vcf 
#' 
#' Will take a vcf with snpeff annotations and output a list of two dataframes.
#'
#' @param vcfpath path to vcf file
#'
#' @returns
#' A list of two dataframes: 
#' - final_snp_impacts: one row per snp, has the highest impact of any transcript for that snp indicated
#' - formatted_snpeff_annotations: one row per transcript. snps will have multiple transcripts identified by snpeff. 
#' 
#' The `IMPACT_PLUS` column has had modifier snps broken out into intergenic and coding modifiers. 
#' @export
#'
#' @examples
#' # work in progress 
format_snpeff_annotations <- function(vcfpath){
  
  path <-
    execute_snpsift2(vcfpath)
  
  message("~~~~~~~~~~~~~~~ Formatting SNP impacts ~~~~~~~~~~~~~~~")
  
  dat.long <- fread(path) %>% 
    mutate(across(starts_with("ANN"), ~ strsplit(.x, split = ","))) %>% 
    tidyr::unnest_longer(col = starts_with("ANN"))
  names(dat.long) <- stringr::str_replace_all(names(dat.long), "ANN\\[\\*\\]\\.", "")
  
  # make a new impact category that splits modifier
  sort <- dat.long %>% 
    ungroup() %>% 
    mutate(marker.ID = paste(.data$CHROM, .data$POS, sep = "-")) %>% 
    mutate(IMPACT_PLUS = case_when(.data$IMPACT == "MODIFIER" & .data$BIOTYPE == "protein_coding" ~ "MODIFIER_CODING",
                                   .data$IMPACT == "MODIFIER" & .data$BIOTYPE == "." ~ "MODIFIER_INTERGENIC",
                                   TRUE ~ .data$IMPACT)) 
  
  #  make a numeric scale to rank our preference for choosing these impacts if a snp has multiple
  scoring.key <- data.frame(IMPACT_PLUS = c("MODIFIER_INTERGENIC", "MODIFIER_CODING", "LOW", "MODERATE", "HIGH"),
                            IMPACT_score = c(1:5))
  
  # retain our favorite impact
  sort_score <- sort %>% 
    left_join(scoring.key, by = "IMPACT_PLUS") %>% 
    group_by(.data$marker.ID) %>% 
    mutate(is.max.impact.score = .data$IMPACT_score == max(.data$IMPACT_score)) %>% 
    filter(.data$is.max.impact.score) %>% 
    select("CHROM", "POS", "marker.ID", "IMPACT", "IMPACT_PLUS") %>% 
    distinct() 
  
  # keep 2 tables
  # table with 1 impact (row) per marker
  # table with 1 row per transcript
  
  if(max(table(sort_score$marker.ID)) > 1){
    warning("Snp impacts table has multiple rows per marker!")
  }
  
  return(list(final_snp_impacts = sort_score,
              formatted_snpeff_annotations = sort))
  
}
