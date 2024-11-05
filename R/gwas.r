#' panvar_gwas
#'
#' @description 
#' This function does gwas on your phenotype and genotype data.
#'
#' @param genotype_data Path to the input genotype data.
#' Either a vcf file (vcf/gz) or bed file
#' @param phenotype_data Path to your genotype data.
#' @param pc_min (optional) What is the minimum number of PCs that should be included in GWAS?
#' Defaults to 5
#' @param pc_max (optional) What is the maximum number of PCs that should be included in GWAS?
#' Defaults to 5
#' @param maf (optional) The minor allele frequency that will be applied to the genotype data.
#' Defaults to 0.05
#' @param missing_rate (optional) The missing_rate filter that will be applied to the genotype data.
#' Defaults to 0.1
#' @param dynamic_correlation (optional) Should the PCs, beyond minimum, be calculated dynamically?
#' Defaults to FALSE
#' @return GWAS results in tabular format.
#'
#' @examples
#' panvar_gwas("path/to/genotype_data",pc_min = 5,pc_max = 5, maf = 0.05, missing_rate = 0.1)
#'
#' @export

# A function to make PCAs for your genotype data.

panvar_gwas <- function(genotype_data,phentotype_path,pc_min = 5,pc_max = 5, maf = 0.05, missing_rate = 0.1, dynamic_correlation = FALSE) {

	# Get the core count from the ergonomics set of code
	core_count = good_core_count()

	# get the extention of the genotype_data file
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
	big_random_pca <- big_randomSVD(
	    the_genotypes, 
	    fun.scaling = snp_scaleBinom(),
	    ncores = core_count
	)

	the_PCs <- predict(big_random_pca)

	the_covariates <- the_PCs[,1:pc_max]

	tryCatch(
	  expr = {
	    # Code that might throw an error
	    phenotype_data <- fread(phentotype_path)
	  },
	  error = function(e) {
	    # Handle the error
	    print(paste("There was difficulty reading your phenotype file, Please check the error message:- ", e))
	  }
	)

	tryCatch(
	  expr = {
	    # Code that might throw an error
	    fam_data <- fread(fam_file_path)
	  },
	  error = function(e) {
	    # Handle the error
	    print(paste("There was difficulty reading your fam file, Please check the error message:- ", e))
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

	pcs_to_include = seq(1:pc_min)

	# If the dynamic correlation is set to TRUE 
	# then the more than the minimum number of PC will be included after they 
	# have been tested for correlation
	if(dynamic_correlation == TRUE) {
		if(pc_max > pc_min){
			for (j in pc_min:pc_max) {
			  if (cor.test(phenotype_scores, the_covariates[, j])[3] < 0.05) {
			    pcs_to_include = c(pcs_to_include, j)
			  }
			}
		}
	}
		
	gwas <- big_univLinReg(
    	the_genotypes, 
    	scale(
    	    phenotype_scores[include_in_gwas]), 
    	    ind.train = include_in_gwas,
    	    covar.train = the_PCs[, pcs_to_include][include_in_gwas, ],
    	    ncores = core_count 
    	)
	
	the_order <- order(the_chromosomes, the_chromosomes)

	pvalues <- stats::predict(gwas)[the_order]

	pvalues <- -1 * pvalues # To just get the -log10 values

	the_gwas <- as.data.table(
		cbind(
			CHROM = the_chromosomes,
			BP = the_bp,
			Pvalues = pvalues
		)
	)

	# convert bp and pvalues to numeric types -
	# to err on the side of caution
	the_gwas[, `:=`(BP = as.integer(BP), Pvalues = as.numeric(Pvalues))]


	# return the gwas table in descending order per pvalues
	return_gwas <- the_gwas %>%
		arrange(desc(Pvalues))
	
	
	return(return_gwas)
}