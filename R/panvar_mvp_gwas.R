
#' Use rMVP to run gwas 
#' 
#' Designed to be used with [panvaR::make_panvar_inputs]. Will read in files and 
#' run GWAS using [rMVP::MVP]. Can also supply inputs as matrices. 
#'
#' @param inputs.dir character, directory to find input files created by [panvaR::make_panvar_inputs]. 
#' Consider using `options()$panvar_outdir`. 
#' @param in.prefix character, prefix of files in input directory. Same as supplied to [panvaR::make_panvar_inputs].
#' Consider using [panvaR::set_panvar_prefix]. 
#' @param npcs numeric, number of principal components to be included in the gwas model. 
#' @param gwas.model character, one of "GLM" or "MLM" to refer to a generalized linear model 
#' and mixed-linear model respectively. 
#' @param output.manhattan boolean, if TRUE, visualizations of gwas results will be output to out.dir
#' @param pheno.mat matrix, optional, object to use for phenotype 
#' @param geno.mat big.matrix, optional, object to use for genotype 
#' @param map.mat matrix, optional, object to use a genotype map file
#' @param pcs.mat big.matrix, optional, object to use as principal component matrix
#' @param kin.mat big.matrix, optional, object to use as kinship matrix
#' @param out.dir character, optional, path to store output. Will overide option set by [panvaR::set_out_dir]
#' @param out.prefix character, optional, a prefix for output files. Will overide option set by [panvaR::set_panvar_prefix].
#'
#' @returns
#' outputs table of gwas results and optionally visualizations produced by 
#' [rMVP::MVP]
#' @export
#'
#' @examples
#' # work in progress
panvar_mvp_gwas <- function(inputs.dir = NULL,
                            in.prefix = NULL,
                            npcs = NULL,
                            gwas.model = c("GLM", "MLM"),
                            output.manhattan = FALSE,
                            pheno.mat = NULL,
                            geno.mat = NULL,
                            map.mat = NULL,
                            pcs.mat = NULL,
                            kin.mat = NULL,
                            out.dir = NULL,
                            out.prefix = in.prefix) {
  
  # check input directory
  if(!is.null(inputs.dir)){
    inputs.dir <- normalizePath(inputs.dir)
  }
  
  # check out.dir
  if(!is.null(out.dir)){
    out.dir <- out.dir
  } else if(!is.null(options()$panvar_outdir)){
    out.dir <- options()$panvar_outdir
  } else {
    out.dir <- tempdir()
  }
  
  # check out prefix
  if(!is.null(out.prefix)){
    out.prefix <- out.prefix
  } else if(!is.null(options()$panvar_prefix)){
    out.prefix <- options()$panvar_prefix
  } else {
    out.prefix <- NULL
  }
  
  # check in prefix
  if(!is.null(in.prefix)){
    in.prefix <- in.prefix
  } else if(!is.null(options()$panvar_prefix)){
    in.prefix <- options()$panvar_prefix
  } else {
    in.prefix <- NULL
  }
  
  # parse model type option
  gwas.model <- arg_match(gwas.model)
  
  if(any(is.null(inputs.dir), is.null(in.prefix))){
    message("Either no input directory, input prefix or neither provided. Getting inputs from user arguments.")
    user.supplied <- TRUE
  } else {
    matched.files <- list.files(path = inputs.dir, pattern = paste0("^", in.prefix))
    if(length(matched.files) == 0){
      stop(paste0("Error in finding input files using input directory: ", inputs.dir, " and input prefix ", in.prefix))
    }
    message(paste0("Searching for prefix: ", out.prefix, " in directory: ", out.dir))
    message(paste0("Found the following files: \n",
                   paste(matched.files, collapse = ", \n")))
    user.supplied <- FALSE
  }
  
  message("Reading in inputs.")
  # read in the matrices if user supplied and make sure they're all here
  if(user.supplied){
    geno <- geno.mat
    pheno <- pheno.mat
    map <- map.mat
    pcs <- pcs.mat
    mats <- list(geno.mat = geno, map.mat = map, pcs.mat = pcs, pheno.mat = pheno)
    if(gwas.model == "GLM"){
      if(any(sapply(mats, is.null))){
        stop(paste0("Must supply all of geno.mat, map.mat and pcs.mat for GLM.
             The following inputs are NULL: ", names(which(sapply(mats, is.null)))))
      }
    } else {
      kin <- kin.mat
      mats <- c(mats, list(kin.mat = kin))
      if(any(sapply(mats, is.null))){
        stop(paste0("Must supply all of geno.mat, map.mat, pcs.mat and kin.mat for MLM.
             The following inputs are NULL: ", names(which(sapply(mats, is.null)))))
      }
    }
  # find the files if they weren't user supplied 
  } else {
    # could check if the files are here first?
    matched.files <- list.files(path = inputs.dir, pattern = paste0("^", in.prefix))
    
    list.files(path = inputs.dir, pattern = paste0("^", in.prefix, ".*", ".geno.desc$"), full.names = T)
    
    geno <- attach.big.matrix(get_an_input(inputs.dir, in.prefix, "geno.desc"))
    pheno <- read.table(get_an_input(inputs.dir, in.prefix, ".phe"), header = T)
    map <- read.table(get_an_input(inputs.dir, in.prefix, "geno.map"), header = T)
    # pcs <- attach.big.matrix(get_an_input(inputs.dir, in.prefix, "pc.desc"))
    if(gwas.model == "MLM"){
      kin.path <- get_an_input(inputs.dir, in.prefix, "kin.desc")
      if(length(kin.path) < 1){
        stop(paste0("Cannot find expected kinship matrix file: ", 
                    paste0(in.prefix, "kin.desc"),
                    " in directory: ", inputs.dir, 
                    ". Was it created when running make_panvar_inputs()?"))
      }
      kin <- attach.big.matrix(get_an_input(inputs.dir, in.prefix, "kin.desc"))
    }
  }
  
  message("Checking phenotype table.")
  # check phenotype table
  
  # pheno.name <- names(phenotype.table)[2]
  # pheno.sub <- phenotype.table %>% 
  #   filter(!is.na(.data[[pheno.name]])) 
  # message(paste("Removed", nrow(phenotype.table) - nrow(pheno.sub), "samples due to NA values in phenotype."))
  # names(pheno.sub)[1] <- "Taxa"
  
  # make output parameters
  if(output.manhattan){
    out.params <- c("plot")
  } else {
    out.params <- NULL
  }
  
  message("Running GWAS")
  if(gwas.model == "GLM"){
    
    mvp.res <- 
    MVP(phe = pheno,
        geno = geno,
        map = map,
        method = "GLM",
        nPC.GLM = npcs,
        maf = NULL,
        file.output = out.params, # consider setting this to false and parsing the output a little more myself 
        outpath = out.dir,
        memo = out.prefix)
    
  } else {
    
    mvp.res <-
    MVP(phe = pheno,
        geno = geno,
        map = map,
        K = kin,
        method = "MLM",
        nPC.MLM = npcs,
        maf = NULL,
        file.output = out.params,
        outpath = out.dir,
        memo = out.prefix)
    
  }
  
  # format gwas results
  dim(mvp.res$glm.results)
  if(gwas.model == "GLM"){
    out <- as.data.frame(cbind(mvp.res$map, mvp.res$glm.results))
  } else {
    out <- as.data.frame(cbind(mvp.res$map, mvp.res$mlm.results))
  }
  
  names(out)[9] <- "PVAL"
  out <- out %>% 
    mutate(LOGPVAL = -log10(.data$PVAL)) %>% 
    rename("CHR" = "CHROM",
           "EFF" = "Effect",
           "marker.ID" = "SNP")
  
  # write out gwas results
  outfullfilename <- paste0(out.prefix, "_", gwas.model, "_GWASresults.csv")
  write.csv(out, file.path(out.dir, outfullfilename), row.names = F)
}