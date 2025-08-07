#' panvar_gwas
#'
#' This function performs genome-wide association studies (GWAS) on provided phenotype and genotype data.
#'
#' @param genotype_data character. Path to the input genotype data file, either in VCF (vcf/gz) or BED format.
#' @param phenotype_input character or data.table. Path to the phenotype data file OR a data.table object containing phenotype data. The first column must contain genotype identifiers.
#' @param pc_min integer, optional. Minimum number of principal components (PCs) to include in GWAS. Default is 5.
#' @param pc_max integer, optional. Maximum number of PCs to include in GWAS. Default is 5.
#' @param maf numeric, optional. Minor allele frequency filter for the genotype data. Default is 0.05.
#' @param missing_rate numeric, optional. Missing rate filter for the genotype data. Default is 0.1.
#' @param dynamic_correlation logical, optional. Whether additional PCs beyond the minimum should be calculated dynamically. Default is FALSE.
#' @param specific_pcs Vector, optional. If you want to supply specific PCs instead of calculating them dynamically then use this to supply a vector of PCs.
#' @return A data frame containing GWAS results. Pvalues represented as -log10(pvalue). 
#'
#' @examples
#' # Using file path
#' panvar_gwas("path/to/genotype_data.vcf", "path/to/phenotype_data.csv", pc_min = 5, pc_max = 5, maf = 0.05, missing_rate = 0.1)
#' # Using data.table object
#' pheno_dt <- data.table::fread("path/to/phenotype_data.csv")
#' panvar_gwas("path/to/genotype_data.vcf", pheno_dt, pc_min = 5, pc_max = 5, maf = 0.05, missing_rate = 0.1)
#'
#' @import tidyverse
#' @import data.table
#' @import parallel
#' @importFrom bigstatsr big_randomSVD
#' @import modelr
#' @importFrom methods is
#'
#' @export
panvar_gwas <- function(genotype_data, phenotype_input, pc_min = 5, pc_max = 5, maf = 0.05, missing_rate = 0.1, dynamic_correlation = FALSE, specific_PCs = NULL) {
  
  print("~~~~~~~~~~~~~~~ Beginning GWAS! ~~~~~~~~~~~~~~~")
  
  # Get the core count from the ergonomics set of code
  core_count = good_core_count()
  
  # get the extension of the genotype_data file
  genotype_data_format <- extention_func(genotype_data)
  
  # if the format is gz or vcf send to plink2
  # else send the genotype data to be filtered for maf and missing rate filter
  if(genotype_data_format == "gz"){
    genotype_data_right_format_base <- vcf_to_plink2(genotype_data) # This returns a list of bed and bim file
    genotype_data_right_format <- genotype_data_right_format_base$bed
    # Now clean up the data
    cleaned_up_bed_file <- bed_file_clean_up(genotype_data_right_format,maf = maf, missing_rate = missing_rate)
  } else if (genotype_data_format == "vcf"){
    genotype_data_right_format_base <- vcf_to_plink2(genotype_data) # This returns a list of bed and bim file
    genotype_data_right_format <- genotype_data_right_format_base$bed
    
    # Now clean up the data
    cleaned_up_bed_file <- bed_file_clean_up(genotype_data_right_format,maf = maf, missing_rate = missing_rate)
  } else {
    
    # if the data is already in plink2 format
    # Clean up for the missing rate and the maf
    cleaned_up_bed_file <- bed_file_clean_up(genotype_data_right_format,maf = maf, missing_rate = missing_rate)
  }
  
  # make the path to the fam file for the genotype data
  cleaned_up_bed_file_path <-
    base_name_func(cleaned_up_bed_file,super_name = TRUE, include_dir = TRUE)
  
  fam_file_path <- paste0(cleaned_up_bed_file_path,".fam")
  
  # Read the bam file into using `snp_readBed`
  
  # but first check to see if the .bk file exists
  bk_base_path <- base_name_func(cleaned_up_bed_file, include_dir = TRUE)
  bk_file_path <- paste0(bk_base_path,".bk")
  
  if(!file.exists(bk_file_path)){
    snp_readBed(cleaned_up_bed_file)
  }
  
  rds_file_base <- base_name_func(cleaned_up_bed_file,include_dir = TRUE)
  
  # the path to the rds file made by `rds_file_base`
  rds_file_path <- paste0(c(rds_file_base,".rds"),collapse = "")
  
  # Read the genotype using `snp_attach`
  genotype_rds_data <- snp_attach(rds_file_path)
  
  # The genotypes
  the_genotypes <- genotype_rds_data$genotypes
  
  # The vector of chromosomes
  the_chromosomes <- genotype_rds_data$map$chromosome
  
  # The vector of input basepairs
  the_bp <- genotype_rds_data$map$physical.pos
  
  # How take the list of existing chromosomes
  # Make a table that mapes the existing chromsomes to a unqiue numerical ID
  # This saves us the trouble chars causing issues
  chromosomes_as_ints <- the_chromosomes %>%
    enframe(name = NULL, value = "chromosome") %>%
    mutate(chr_int = dense_rank(chromosome)) %>%
    pull(chr_int)
  
  # Make the exclusion index
  # Rijan: Reading material clumping and pruning here https://www.biostars.org/p/343818/
  
  # Apply a PCA using an algorithm optimized for large file backed matrices
  big_random_pca <- bigsnpr::snp_autoSVD(the_genotypes,
                                         infos.chr = chromosomes_as_ints,
                                         infos.pos = the_bp,
                                         fun.scaling = snp_scaleBinom())
  
  the_PCs <- predict(big_random_pca)
  
  the_covariates <- the_PCs[,1:pc_max]
  
  # Check if phenotype_input is a path or a data.table
  if (is.character(phenotype_input)) {
    if (!file.exists(phenotype_input)) {
      stop("The phenotype file path provided does not exist: ", phenotype_input)
    }
    tryCatch(
      expr = {
        phenotype_data <- data.table::fread(phenotype_input)
      },
      error = function(e) {
        stop("Error reading phenotype file: ", phenotype_input, "\nOriginal error: ", e)
      }
    )
  } else if (is(phenotype_input, "data.table")) {
    # Check if it's specifically a data.table, not just data.frame
    if (!inherits(phenotype_input, "data.table")) {
      warning("phenotype_input is a data.frame, converting to data.table.")
      phenotype_data <- data.table::as.data.table(phenotype_input)
    } else {
      phenotype_data <- phenotype_input
    }
  } else {
    stop("phenotype_input must be either a valid file path (character string) or a data.table object.")
  }
  
  tryCatch(
    expr = {
      # Code that might throw an error
      fam_data <- data.table::fread(fam_file_path)
    },
    error = function(e) {
      # Handle the error
      print(paste("There was difficulty reading your fam file, Please check the error message:- ", e))
      stop("Error reading FAM file: ", fam_file_path) # Stop execution if FAM file cannot be read
    }
  )
  
  fam_data <- fam_data %>%
    rename(genotype = 2)
  
  # Stringently assumes that the line names are the first field of data
  phenotype_data <- phenotype_data %>%
    rename(genotype = 1)
  
  pheno_data_for_gwas <- left_join(fam_data, phenotype_data, by = "genotype")
  
  ordered_phenotype <- pheno_data_for_gwas %>%
    slice(match(genotype_rds_data$fam$sample.ID, pheno_data_for_gwas$genotype))
  
  genotype_rds_data$fam <- cbind(genotype_rds_data$fam, ordered_phenotype)
  
  # This assumes that the last column
  # of the supplied table is phenotype data
  phenotype_scores <- genotype_rds_data$fam %>%
    pull(ncol(genotype_rds_data$fam))
  
  include_in_gwas <- which(!is.na(phenotype_scores) & complete.cases(the_covariates))
  
  if(!is.null(specific_PCs)){
    pcs_to_include = specific_PCs
  } else {
    pcs_to_include = seq(1:pc_min)
    
    # If the dynamic correlation is set to TRUE
    # then the more than the minimum number of PC will be included after they
    # have been tested for correlation
    if(dynamic_correlation == TRUE) {
      if(pc_max > pc_min){
        for (j in pc_min:pc_max) {
          # Check correlation and handle potential errors/warnings in cor.test
          cor_result <- tryCatch(cor.test(phenotype_scores[include_in_gwas], the_covariates[include_in_gwas, j]), warning = function(w) w, error = function(e) e)
          if (!inherits(cor_result, "warning") && !inherits(cor_result, "error") && cor_result$p.value < 0.001) {
            pcs_to_include = c(pcs_to_include, j)
          } else if (inherits(cor_result, "warning") || inherits(cor_result, "error")) {
            warning(paste("Correlation test for PC", j, "failed or produced a warning. Skipping PC."))
          }
        }
      } else {
        # This condition should ideally check pc_max >= pc_min for dynamic calculation range
        if (pc_max < pc_min) {
          warning("pc_max is less than pc_min. Cannot perform dynamic PC selection beyond pc_min.")
        } else {
          # If pc_max == pc_min, no dynamic selection is needed beyond the initial pcs_to_include
          print("Note: pc_max equals pc_min, dynamic PC selection range is zero.")
        }
      }
    }
  }
  
  print(paste("GWAS model will include the following PC's: ", paste(pcs_to_include, collapse = ",")))
  print(paste("Running model with", length(include_in_gwas), "genotypes and", nrow(genotype_rds_data$map), "snps."))
  
  # ind_u <- matrix(PC[genoLineIndx,pcs_to_include], ncol = length(pcs_to_include))
  
  gwas <- big_univLinReg(
    the_genotypes,
    scale(
      phenotype_scores[include_in_gwas]),
    ind.train = include_in_gwas,
    covar.train = the_PCs[, pcs_to_include, drop = FALSE][include_in_gwas, , drop = FALSE], # Added drop=FALSE for robustness
    ncores = 1
  )
  
  # Use map data directly for order, assuming it corresponds to the genotypes matrix columns
  map_data <- genotype_rds_data$map
  the_order <- order(map_data$chromosome, map_data$physical.pos) # Order by chromosome then position
  
  
  pvalues <- predict(gwas, log10 = TRUE) # Raw scores/stats from regression
  
  # Convert stats to p-values if predict doesn't return p-values directly
  # Assuming predict returns log10(p) or similar stat that needs conversion
  pvalues <- -1 * pvalues # Example conversion if predict gives -log10(p)
  
  # Create the GWAS results table using ordered map data
  the_gwas <- data.table(
    CHROM = map_data$chromosome[the_order],
    BP = map_data$physical.pos[the_order],
    Pvalues = pvalues[the_order] # Apply order to p-values
  )
  
  # Ensure correct types
  the_gwas[, `:=`(BP = as.integer(BP), Pvalues = as.numeric(Pvalues))]
  
  
  # return the gwas table in descending order per pvalues
  return_gwas <- the_gwas %>%
    arrange(CHROM, BP) 
  
  print("GWAS completed!")
  
  return(return_gwas)
}