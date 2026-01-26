# pheno <- read.table(paste0(in.prefix, ".phe"), header = T)
# geno <- attach.big.matrix(paste0(in.prefix, ".geno.desc"))
# map <- read.table(paste0(in.prefix, ".geno.map"), header = T)
# pcs <- attach.big.matrix(paste0(in.prefix, ".pc.desc"))
# kin <- attach.big.matrix(paste0(in.prefix, ".kin.desc"))
# 
# 
# MVP(phe = pheno.sub[,c(1,2)],
#     geno = geno,
#     map = map,
#     K = kin,
#     nPC.GLM = num.pcs,
#     nPC.MLM = num.pcs,
#     nPC.FarmCPU = num.pcs,
#     p.threshold = effective.bonf,
#     QTN.threshold = .01, # when to include snps in fcpu model, default = .01
#     method.bin= "FaST-LMM",
#     priority = "speed",
#     vc.method = "EMMA",
#     # method = c("GLM", "FarmCPU"), # original
#     method = c("FarmCPU", "MLM", "GLM"), # flip order
#     outpath = paste0("/scratch/gwas_out/mvp_out/tables/", pheno.name),
#     file.output = c("pmap", "pmap.signal", "log"),
#     verbose = T,
#     # cutoff for "signals" file and for drawing and plots. plots only plotted if "plot" in file.output vector
#     # cutoff is the alpha value (gets divided by marker number), so if want specific pvalue have to do this
#     threshold = effective.bonf * nrow(geno), 
#     memo = this.memo)


panvar_mvp_gwas <- function(inputs.dir = NULL,
                            in.prefix = NULL,
                            npcs = NULL,
                            gwas.model = c("GLM", "MLM"),
                            output.manhattan = FALSE,
                            phenotype.mat = NULL,
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
    pcs <- attach.big.matrix(get_an_input(inputs.dir, in.prefix, "pc.desc"))
    if(gwas.model == "MLM"){
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