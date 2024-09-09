check_plink2_chroms <- function(vcf_file_path,bim_file_path,keep_table) {

    # Function to extract integers
    extract_integers <- function(x) {
      as.numeric(gsub("\\D", "", x))
    }
    
    bcf_chroms <- bcf_chroms_func(vcf_file_path)

    bim_file_chroms <- fread(bim_file_path) %>% 
        pull(1) %>% 
        unique()
    
    chrom_test <- all(bim_file_chroms %in% bcf_chroms)

    if(chrom_test){
        keep_table_path <- keep_table_sanitizer(keep_table)
        return(keep_table_path)
    } else{

        print("When you used Plink2 to make a bed file is forcefully changed the names of the - PanvaR will now use a heuristic to work around this.")
        print("Please read about it here: https://www.biostars.org/p/9602116/")

          bcf_chroms_table <- data.table(
            chrom_names = bcf_chroms,
            chrom_ints = sapply(bcf_chroms, extract_integers)
	        )

        plink2_chroms <- as.data.table(bim_file_chroms)
    
        vcf_plink2_dictionary_table <- left_join(plink2_chroms,bcf_chroms_table, by = c( "bim_file_chroms" = "chrom_ints" ))
    
        new_keep_table <- left_join(keep_table,vcf_plink2_dictionary_table,by = c("V1"="bim_file_chroms")) %>%
        select(chrom_names, V2)
    
        return(keep_table_sanitizer(new_keep_table))
    }
    
    
}