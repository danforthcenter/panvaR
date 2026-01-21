make_panvar_inputs <- function(genotype.path,
                               phenotype.table, # two columns, 1: linenames 2: phenotype
                               min.maf = .05,
                               max.missing.snp = .1,
                               calc.kinship = F,
                               plink.path = options()$plink_path,
                               out.dir = options()$panvar_outdir){
  
  # make output directory
  out.dir.path <- temporary_directory()
  # store phenotype name
  pheno.name <- names(phenotype.table)[2]
  
  # ~~~~ QC phenotype ~~~~
  
  # any na's? 
  pheno.sub <- phenotype.table %>% 
    filter(!is.na(.data[[pheno.name]]))
  message(paste("Removed", nrow(phenotype.table) - nrow(pheno.sub), "samples due to NA values in phenotype."))

  # store list of genotypes in study
  samples <- pheno.sub[,1]
  samples.out.file <- paste0(out.dir.path, "/Panvar_list.of.samples.with.phenotype_", pheno.name, ".txt")
  write.table(samples, samples.out.file, quote = F, header = F, row.names = F)
  
  # ~~~~ QC genotype ~~~~~
  
  # filter for:
  #   - maf
  #   - snp missing rate
  #   - pheno'd genos
  # convert to bed if in vcf
  
  # make new function
  
  tryCatch(
    {
      error_message <- tempfile()
      try <- exec_wait(
        binary_call,
        args = binary_args,
        std_out = TRUE,
        std_err = error_message
      )
    },
    error = function(e){
      # Custom error message
      print(paste("Execution attempt produced error:-", e$message))
      1 # Return 1 on error
    }
  )
  
  
  
  if(try == 0){
    
    final_output_path = paste0(output_path,".bed")
    
    print(
      paste("The Plink files are available in",final_output_path)
    )
    return(final_output_path)
  } else{
    
    print("There were errors when applying MAF and missing rate filter to your BED file.")
    print("Please read the error message and re-try")
    return(readLines(error_message))
  }
  
  # ~~~~ generate PC's and Kinship ~~~~ 
  
  
}