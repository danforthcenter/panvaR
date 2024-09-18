plink2_bcftools_chroms_dictionary <- function(vcf_file_path,bim_file_path) {

    # Function to extract integers
    extract_integers <- function(x) {
      as.numeric(gsub("\\D", "", x))
    }
    
    bcf_chroms <- bcf_chroms_func(vcf_file_path)

    bim_file_chroms <- fread(bim_file_path) %>% 
        pull(1) %>% 
        unique()
    
    bcf_chroms_table <- data.table(
      chrom_names = bcf_chroms,
      chrom_ints = sapply(bcf_chroms, extract_integers)
	  )

    chrom_test <- all(bim_file_chroms %in% bcf_chroms)

    if(!chrom_test){
      plink2_chroms <- as.data.table(bim_file_chroms)

      vcf_plink2_dictionary_table <- left_join(
          plink2_chroms,
          bcf_chroms_table, 
          by = c( "bim_file_chroms" = "chrom_ints" )
      ) %>% rename(
          plink2 = bim_file_chroms,
          vcf = chrom_names
      )

      return(vcf_plink2_dictionary_table)
    } else {
      return(NULL) # NULL is a substitute for TRUE here
    }
    
}

apply_dict <- function(dictionary, table){

    # The function assumes that -
    # The dictionary table has two fields - 
    # 1. plink2, the chromosomes in plink2 -
    # 2. vcf, the chromosomes in vcf

    # The table should have the following fields -
    # CHROM

    fixed_table <- table %>% 
        left_join(dictionary, by = c("CHROM"="plink2")) %>% # fix the CHROM's by matching them to plink2 names
        select(-CHROM) %>%
        rename(CHROM = vcf) # rename vcf column to CHROMs
    
    return(fixed_table)
}