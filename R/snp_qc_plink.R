#' Quality control a snp file using plink2
#'
#' @param genotype.path character, path to genotype file, supported types: '.bed', .'vcf', '.vcf.gz'.
#' @param min.maf numeric, filtering cutoff for minor allele frequency, snps are removed if they have maf less than this value. To ignore set to 0.
#' @param max.missing.snp numeric, filtering cutoff for missing rate of snps, snps are removed if they have a missing rate higher than this. To ignore set to 1.
#' @param sample.list.path character, optional, path to a list of samples. 
#' Samples in file will be included. Sample filtering happens before other filtering per plink's order of operations.
#' @param plink.path character, optional, path to plink2 executable. If not provided, will default to option set by [panvaR::set_plink_path].
#' @param out.dir character, optional, path to output files. If not provided, will default to option set by [panvaR::set_out_dir]
#' @param out.name character, optional, prefix for files output. 
#'
#' @returns
#' filtered bed/bim/bam files stored in out.dir. 
#' character string of path to bed file.
#' @export
#'
#' @examples
#' # work in progress
#' 
snp_qc_plink <- function(genotype.path,
                         min.maf = .05, # to turn off set to 0
                         max.missing.snp = .1, # to turn off set to 1
                         sample.list.path = NULL,
                         plink.path = NULL,
                         out.dir = NULL,
                         out.name = NULL){
  
  # ~~~~ Initialize ~~~~
  
  # check plink.path
  if (!is.null(options()$plink_path)) {
    plink.exec <- options()$plink_path
  } else {
    plink.exec <- "plink2"
  }
  
  # check output directory
  out.dir.path <- temporary_directory(out.dir)
  
  # determine filetype that was input
  filetype <- get_geno_filetype(genotype.path)
  
  # make wellformed path
  genotype.path <- normalizePath(genotype.path)
  
  # name output files
  maflabel <- paste0("maf", min.maf)
  missinglabel <- paste0("missing",max.missing.snp)
  if(is.null(out.name)){
    outfullname <- paste("PlinkQC", maflabel, missinglabel, sep = "_")
  } else {
    outfullname <- paste(out.name, "PlinkQC", maflabel, missinglabel, sep = "_")
  }
  
  # make full output path
  out_fullpath <- file.path(out.dir.path, outfullname)
  
  print(out_fullpath)
  
  # ~~~~ Store plink2 arguments ~~~~
  
  if(filetype == "bed"){
    geno.file <- str_replace(genotype.path, "\\.bed$", "") 
    current_args <- c(
      "--allow-extra-chr",
      "--bfile", geno.file,
      "--geno", max.missing.snp,
      "--maf", min.maf,
      "--make-bed",
      "--set-all-var-ids", "@-#",
      "--out", out_fullpath
    )
  } else if(filetype == "vcf"){
    # plink wants the vcf file extension but not the bed one, leave this out
    # geno.file <- str_replace(genotype.path, "\\.vcf.*$", "") 
    current_args <- c(
    "--allow-extra-chr",
    "--vcf", genotype.path,
    "--geno", max.missing.snp,
    "--maf", min.maf,
    "--make-bed",
    "--set-all-var-ids", "@-#",
    "--out", out_fullpath
    )
  } else {
    stop(paste0("Unsupported genotype file type for file: ", genotype.path, ". \n",
                "Supported filetypes: '.bed', '.vcf', '.vcf.gz' \n",
                "If file is in one of these formats, please ensure proper suffix exists in filename."))
  }
  
  # add in optional ability to subset based on samples
  if(!is.null(sample.list.path)){
    sample.list.path.norm <- normalizePath(sample.list.path)
    current_args <- c(
      current_args,
      "--keep", sample.list.path.norm
    )
  }

  # ~~~~ Run plink2 ~~~~
  
  tryCatch(
    {
      error_message <- tempfile()
      try <- exec_wait(
        cmd = plink.exec,
        args = current_args,
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
    
    final_output_path <- paste0(out_fullpath,".bed")
    
    message(
      paste("QC was successful, output stored at",final_output_path)
    )
    return(final_output_path)
  } else{
    
    message("SNP QC failed, plink produced the following error message:")
    return(readLines(error_message))
  }
}