#' Make standard inputs for panvaR
#'
#' Does some filtering of the genotype file for samples that have a phenotype
#' and also minor allele frequency and snps with high missing rates. Will also
#' calculate prinicpal components and optionally kinship matrix of the genotype
#' file for use downstream in gwas.
#'
#' @param genotype.path character, path to genotype file, supported types:
#'   '.bed', .'vcf', '.vcf.gz'.
#' @param phenotype.path character, path to table of phenotype to test. Expects
#'   samples (lines) in column 1 and phenotype in column 2. This is used to
#'   determine the set of samples (lines) to use in the analysis.
#' @param min.maf numeric, filtering cutoff for minor allele frequency, snps are
#'   removed if they have maf less than this value. To ignore set to 0.
#' @param max.missing.snp numeric, filtering cutoff for missing rate of snps,
#'   snps are removed if they have a missing rate higher than this. To ignore
#'   set to 1.
#' @param calc.kinship boolean, optional, if TRUE, the kinship matrix will be
#'   calculated for use in mixed linear model gwas.
#' @param plink.path character, optional, path to plink2 executable. Will
#'   overide option set by [panvaR::set_plink_path].
#' @param out.dir character, optional, path to store output. Will overide option
#'   set by [panvaR::set_out_dir].
#' @param out.prefix character, optional, a prefix for output files. Will
#'   overide option set by [panvaR::set_panvar_prefix].
#' @param extra.plink.options character, a vector of options to include in call
#'   to plink2. Should be a vector with plink2 arguments and their values as
#'   separate elements of vector. E.G. c("--max-maf", ".95", "--max-alleles",
#'   "2"). see [panvaR::snp_qc_plink]
#'
#' @returns Input files to be used for downstream panvaR functions. Stored in
#' `out.dir` or the option set in [panvaR::set_out_dir]. Runs
#' [panvaR::snp_qc_plink] to filter for maf and missing using plink2 and then
#' [rMVP::MVP.Data] to prepare data for GWAS.
#'
#' @export
#'
#' @examples
#' # work in progress
#'
#' # specify some paths
#' genotype.path <- system.file(
#'   "extdata",
#'   "Setaria_shattering_example_pruned.bed",
#'   package = "panvaR")
#' 
#' phenotype.path <- system.file(
#'   "extdata",
#'   "Setaria_shattering_example_phenotype.tsv",
#'   package = "panvaR")
#' 
#' plink.path <- bigsnpr::download_plink2()
#'
#' out.dir <- file.path(tempdir(), "panvar_ex")
#' dir.create(out.dir, showWarnings = FALSE)
#'
#' # run the function
#' make_panvar_inputs(
#'   plink.path = plink.path,
#'   genotype.path = genotype.path,
#'   phenotype.path = phenotype.path,
#'   out.prefix = "example",
#'   out.dir = out.dir)
#'
#' # created some input files
#' list.files(out.dir)
#'
#' # clean up
#' unlink(out.dir, recursive = TRUE)

make_panvar_inputs <- function(genotype.path,
                               phenotype.path, # two columns, 1: linenames 2: phenotype
                               min.maf = .05,
                               max.missing.snp = .1,
                               calc.kinship = F,
                               plink.path = NULL,
                               out.dir = NULL,
                               out.prefix = NULL,
                               extra.plink.options = NULL){
  
  # ~~~~ Initialize ~~~~
  
  # read in phenotype file
  phenotype.table <- data.table::fread(phenotype.path)
  
  # store phenotype name
  pheno.name <- names(phenotype.table)[2]
  
  # check plink.path
  if(!is.null(plink.path)){
    plink.exec <- plink.path
  } else if(!is.null(options()$plink_path)){
    plink.exec <- options()$plink_path
  } else {
    plink.exec <- "plink2"
  }
  
  # check output directory
  out.dir.path <- temporary_directory(out.dir)
  
  # check prefix
  if(!is.null(out.prefix)){
    out.prefix <- out.prefix
  } else if(!is.null(options()$panvar_prefix)){
    out.prefix <- options()$panvar_prefix
  } else {
    out.prefix <- NULL
  }
  
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
      out.prefix = out.prefix,
      extra.options = extra.plink.options
    )

  # ~~~~ generate PC's and Kinship ~~~~ 
  if(calc.kinship){
    message("Using rMVP to calculate PC's and kinship matrix.")
  } else {
    message("Using rMVP to calculate PC's.")
  }
  
  # make our output
  if(is.null(out.prefix)){
    mvp.out.path <- file.path(out.dir.path, "PanvarIN")
  } else {
    mvp.out.path <- file.path(out.dir.path, out.prefix)
  }
  
  # get phenotype file separater
  # phefile <- basename(phenotype.path)
  phe_ext <- tools::file_ext(phenotype.path)
  if(phe_ext == "tsv"){
    sep.phe <- "\t"
  } else if (phe_ext == "csv"){
    sep.phe <- ","
  } else {
    stop("Phenotype file does not have either .tsv or .csv extension.
         Only tab or comma separated tables are accepted. Please ensure proper 
         file extension is present.")
  }
  
  rMVP::MVP.Data(fileBed = snp_qc_path,
                 filePC = F,
                 fileKin = calc.kinship,
                 filePhe = phenotype.path,
                 sep.phe = sep.phe, 
                 pcs.keep = 10,
                 SNP.impute = NULL,
                 out = mvp.out.path)
}