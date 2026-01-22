make_panvar_inputs <- function(genotype.path,
                               phenotype.table, # two columns, 1: linenames 2: phenotype
                               min.maf = .05,
                               max.missing.snp = .1,
                               calc.kinship = F,
                               plink.path = NULL,
                               out.dir = NULL){
  
  # ~~~~ Initialize ~~~~
  
  # store phenotype name
  pheno.name <- names(phenotype.table)[2]
  # check plink.path
  if (!is.null(options()$plink_path)) {
    plink.exec <- options()$plink_path
  } else {
    plink.exec <- "plink2"
  }
  # check output directory
  out.dir.path <- temporary_directory(out.dir)
  
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
  
  
  # ~~~~ generate PC's and Kinship ~~~~ 
  
  
}