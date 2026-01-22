make_panvar_inputs <- function(genotype.path,
                               phenotype.table, # two columns, 1: linenames 2: phenotype
                               min.maf = .05,
                               max.missing.snp = .1,
                               calc.kinship = F,
                               plink.path = NULL,
                               out.dir = NULL,
                               out.name = NULL){
  
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
  write.table(samples, samples.out.file, quote = F, col.names = F, row.names = F)
  
  # ~~~~ QC genotype ~~~~~
  
  # filter for:
  #   - maf
  #   - snp missing rate
  #   - pheno'd genos
  # convert to bed if in vcf
  snp_qc_path <-
    snp_qc_plink(
      genotype.path,
      min.maf = min.maf,
      max.missing.snp = max.missing.snp,
      sample.list.path = samples.out.file,
      plink.path = plink.path,
      out.dir = out.dir,
      out.name = out.name
    )

  # ~~~~ generate PC's and Kinship ~~~~ 
  if(calc.kinship){
    message("Using rMVP to calculate PC's and kinship matrix.")
  } else {
    message("Using rMVP to calculate PC's.")
  }
  
  # make our output
  if(is.null(out.name)){
    mvp.out.path <- file.path(out.dir.path, "PanvarIN")
  } else {
    mvp.out.path <- file.path(out.dir.path, out.name)
  }
  
  rMVP::MVP.Data(fileBed = snp_qc_path,
                 filePC = TRUE,
                 fileKin = calc.kinship,
                 pcs.keep = 10,
                 out = mvp.out.path)
}