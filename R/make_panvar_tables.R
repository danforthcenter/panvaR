#' Make some standardized tables from gwas and annotation tables
#'
#' @param gwas.res data.frame of all gwas results, should contain columns (CHR,
#'   POS, PVAL), corresponding to (chromosome, physical position, and pvalue).
#' @param qtl.df data.frame, table that includes list of snps to calculate LD to
#'   with columns (CHR, POS, LOGPVAL), corresponding to (chromosome, physical
#'   position, and -log10(p-value)). QTL are typically defined as hits grouped
#'   by LD by something like `plink --clump`. See [panvaR::get_ld_in_window]
#' @param tag.snp character, marker.ID of snp around which to calculate LD. In
#'   the form 'CHR-POS'
#' @param annotation.table table with annotations with columns (geneID, CHR,
#'   start, end, annotation). start and end correspond to base-pair coordinates
#'   of start and end of gene. CHR is chromosome of gene.
#' @param plink.path character, optional, path to plink2 executable. Will
#'   overide option set by [panvaR::set_plink_path].
#' @param pvals.in.log boolean, if TRUE PVAL column has already been converted to -log10(pvalue)
#' @param geno.bed.filename character, prefix of genotype files in plink
#'   (bed/bim/fam) format. Do not include ".bed" extension.
#' @param geno.bed.directory character, directory where genotype files are
#'   located
#' @param temp.dir character, where to output some temporary files.
#' @param window numeric, total window size in KB, all variants within .5 *
#'   window are calculated.
#' @param snp.to.gene.vars character, numeric variables in gwas.res to aggregate by gene. 
#' For each gene, snps with a physical position with the start and end of the gene are considered. 
#' The maximum value for all snps within the gene is returned. Special values, `DIST`, `LD` and `LOGPVAL` can
#' be included in addition to any user supplied variables. 
#' @param snp.to.gene.buffer numeric, kilobases to add to gene start and end to include genes 
#' that are close but not in gene. Snpeff uses 5 KB by default to call a snp "upstream"/"downstream" variant. default is 0. 
#' @param compute.scores boolean, if TRUE, snp scores will be computed. See details for more info.  
#' @param score.vars character, vector of column names indicating which variables to included in the score. 
#' If compute.scores is TRUE and score.vars is NULL, the default score will use equally weighted variables: "DIST", "LOGPVAL", "LD". 
#' @param score.dirs numeric, a vector indicating which direction is to be considered more indicative
#' of an association. 1 indicates higher is better, -1 indicates lower is better. The order should correspond 
#' with the order in cols. 
#' @param score.weights numeric, a vector indicating weights for the variables. These must add up to 1. 
#'
#' @details
#' Scores:
#' Scores are simple scaled and weighted averages of some variables. First variables are 
#' normalized using min/max normalization. The variables are then made negative if 
#' they need to be reversed to indicate a larger value as a more desirable value. 
#' 
#' For example, distance from the key snp should be reversed as a small distance is more desirable.
#' A log-pvalue is already of this form 'bigger is better' so does not need to be altered. 
#' 
#' Finally, a weighted average is taken based on user defined weights. The default weights all variables equally. 
#' The outcome is a score from 0-1 that ranks the snps based on these variables. 
#' see: [panvaR::make_scores]
#' 
#' @returns
#' A named list with the following entries:
#' - gwas: formatted gwas results.
#' - anno: formatted annotation results.
#' - key.snp: tag.snp or in the qtl.df case the highest p-value snp supplied to the function for downstream use.
#' - qtl.df: qtl.df supplied to the program if used. 
#' @export
#'
#' @examples
#' # work in progress
make_panvar_tables <- function(gwas.res,
                               qtl.df = NULL,
                               tag.snp = NULL,
                               annotation.table,
                               plink.path = NULL,
                               pvals.in.log = T,
                               geno.bed.filename,
                               geno.bed.directory = "/.",
                               temp.dir = tempdir(),
                               window,
                               snp.to.gene.vars = "LD",
                               snp.to.gene.buffer = 0,
                               compute.scores = F,
                               score.vars = NULL,
                               score.dirs = NULL,
                               score.weights = NULL){
  
  
  # ------------------------------------------------------------------------\
  # calculate LD --------
  # ------------------------------------------------------------------------\
  
  message("Calculating LD")
  # check temp directory
  
  
  ld.list <- get_ld_in_window(qtl.df = qtl.df,
                              tag.snp = tag.snp,
                              window = window,
                              plink.path = plink.path,
                              geno.bed = geno.bed.filename,
                              in.dir = geno.bed.directory,
                              out.dir = temp.dir,
                              verbose = T)
  
  # ------------------------------------------------------------------------\
  # subset gwas  --------
  # ------------------------------------------------------------------------\
  
  # filter gwas df to just in window and join LD
  gwas.sub <- gwas.res %>%
    as.data.frame() %>% 
    # mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
    left_join(ld.list$table, by = c("CHR", "POS")) %>%
    select(-contains("marker.ID")) %>% 
    mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>% 
    filter(!is.na(.data$R2)) %>% 
    # rename this here to work with special case "LD" given to point variable
    rename("LD" = "R2") %>% 
    relocate("marker.ID")
  
  # add middlesnp back in, gets ignored by LD calc if using tagsnp
  if(!is.null(tag.snp)){
    gwas.sub.mid.snp <- gwas.res %>% 
      filter(.data$marker.ID == ld.list$key.snp) %>% 
      mutate(LD = 1)
    gwas.sub <- bind_rows(gwas.sub, gwas.sub.mid.snp) 
  }
  
  
  # ------------------------------------------------------------------------\
  # make scores --------
  # ------------------------------------------------------------------------\
  
  # make_scores <- function(input.df, cols, directions, weights = NULL){
  if(compute.scores){
    message("Computing scores")
    # add dist column 
    score.in <- gwas.sub %>% 
      mutate(key.snp.pos = get_bp_from_id(ld.list$key.snp)) %>% 
      mutate(DIST = abs(.data$POS - .data$key.snp.pos)) 
    
    if(!pvals.in.log){
      score.in <- score.in %>%
        mutate(LOGPVAL = -log10(.data$PVAL))
    } else {
      score.in <- score.in %>% 
        rename("LOGPVAL" = "PVAL")
    }
    
    # use default scoring if nothing specific selected
    if(is.null(score.vars)){
      message("No user defined variables for scores specified. Using equally weighted Distance, Pvalue and LD.")
      score.vars <- c("DIST", "LOGPVAL", "LD")
      score.dirs <- c(-1, 1, 1)
      weights <- c(1/3, 1/3, 1/3)
    } else {
      # use user defined stuff
      message(paste0("Using ", paste(score.vars, collapse = ", "), " to calculate."))
      if(is.null(score.dirs)){
        message("No score directions specified. Assuming all variables are of the style 'larger is better'.")
        score.dirs <- rep.int(1, length(score.vars))
      }
      if(is.null(score.weights)){
        message("No score weights specified, weighting all equally.")
        score.weights <- rep.int(1, length(score.vars)) / length(score.vars)
      }
      dirs.df <- data.frame(variable = score.vars, direction = score.dirs, weight = score.weights)
      message("Computing scores using following parameters:")
      message(paste0(capture.output(knitr::kable(dirs.df)), collapse = "\n"))
    }
    
    scores.df <- 
      make_scores(score.in,
                  cols = score.vars,
                  directions = score.dirs,
                  weights = score.weights)
    gwas.res <- gwas.res %>% 
      left_join(scores.df, by = "marker.ID")
    gwas.sub <- gwas.sub %>% 
      left_join(scores.df, by = "marker.ID")
  }
  
  # ------------------------------------------------------------------------\
  # subset annotation --------
  # ------------------------------------------------------------------------\
  
  # filter anno to just window
  this.chrom <- get_chrom_from_id(ld.list$key.snp)
  this.pos <- get_bp_from_id(ld.list$key.snp)
  if(!is.numeric(annotation.table$CHR)){
    stop("Annotation table column 'CHR' must be numeric.")
  }
  anno.sub <- annotation.table %>%
    filter(.data$CHR == this.chrom) %>%
    rowwise() %>%
    mutate(dist.from.snp = get.gene.dist.from.snp(this.pos, .data$start, .data$end)) %>%
    filter(.data$dist.from.snp <= window * 1000) 
  
  # ------------------------------------------------------------------------\
  # make snp to gene stats --------
  # ------------------------------------------------------------------------\
  
  message("Generating snp to gene correspondence")
  
  if(!is.null(snp.to.gene.vars)){
    # get the snp to gene correspondence for whatever you want
    point.color.stat <- gwas.sub %>% 
      rowwise() %>% 
      mutate(snp.in.gene = get.gene.from.snp(.data$POS, anno.sub, snp.to.gene.buffer)) %>% 
      filter(!is.null(.data$snp.in.gene)) %>% 
      unnest_longer(.data$snp.in.gene) %>% 
      group_by(.data$snp.in.gene) %>% 
      # summarize(maximum.value = max(.data[[annotation.point.variable]])) %>% 
      summarize(across(all_of(snp.to.gene.vars), ~ max(.x, na.rm = T)))  %>% 
      rename("geneID" = "snp.in.gene")
    
    anno.out <- anno.sub %>% 
      left_join(point.color.stat, by = "geneID")
  } else {
    anno.out <- anno.sub
  }
  
  
  # ------------------------------------------------------------------------\
  # output tables --------
  # ------------------------------------------------------------------------\
  
  out <- list(gwas = gwas.sub,
              anno = anno.out,
              key.snp = ld.list$key.snp,
              qtl.snps = ld.list$qtl.snps)
  
  return(out)
}