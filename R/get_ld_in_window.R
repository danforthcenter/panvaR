#' Get linkage disequilibrium 
#' 
#' Get LD from a group of snps or a single tag snp to all snps in a window.
#' 
#' When supplying a group a snps, the ld retained for a given snp is the maximum LD of the given snp to the snps in the group.  
#'
#' @param qtl.df data.frame, table that includes list of snps to calculate LD to with columns (CHR, POS, LOGPVAL), corresponding to (chromosome, physical position, and -log10(p-value)). 
#' QTL are typically defined as hits grouped by LD by something like `plink --clump`
#' @param tag.snp character, marker.ID of snp around which to calculate LD. In the form 'CHR-POS'
#' @param window numeric, total window size in KB, all variants within .5 * window are calculated. 
#' @param plink.path character, optional, path to plink2 executable. Will overide option set by [panvaR::set_plink_path].
#' @param geno.bed character, prefix of genotype files in plink (bed/bim/fam) format. Do not include ".bed" extension.
#' @param in.dir character, directory where genotype files are located
#' @param out.dir character, where to output some temporary files. 
#' @param verbose boolean, if TRUE, output some status reports
#'
#' @returns
#' Named list with 2 items
#'  - table: table with marker.IDs (CHR-POS) and maximum LD in R2 for each snp to the snps in the qtl.df or LD to tag.snp
#'  - key.snp: marker.ID corresponding to the middle of the window, either the max(LOGPVAL) of qtl.df or tag.snp. useful to retain for downstream functions.

#' @export
#'
#' @examples
#' # work in progress
get_ld_in_window <- function(qtl.df= NULL,
                             tag.snp = NULL,
                             window,
                             plink.path = NULL,
                             geno.bed,
                             in.dir,
                             out.dir = NULL,
                             verbose = TRUE){
  
  # check plink.path
  if(!is.null(plink.path)){
    plink.exec <- plink.path
  } else if(!is.null(options()$plink_path)){
    plink.exec <- options()$plink_path
  } else {
    plink.exec <- "plink2"
  }
  
  # check out.dir
  if(!is.null(out.dir)){
    out.dir <- out.dir
  } else if(!is.null(options()$panvar_outdir)){
    out.dir <- options()$panvar_outdir
  } else {
    out.dir <- tempdir()
  }
  out.dir <- normalizePath(out.dir)
  
  # check in.dir 
  in.dir <- normalizePath(in.dir)
  
  # check for .bed extension on input
  if(grepl("\\.bed$", geno.bed)){
    geno.bed <- tools::file_path_sans_ext(geno.bed)
  }
  
  # check that bedfile exists
  if(!file.exists(file.path(in.dir, paste0(geno.bed, ".bed")))){
    stop("Bed file not found in input directory.")
  }
  
  # calc ld either using 1 tag snp or multiple snps in a qtl
  if(!is.null(qtl.df) & !is.null(tag.snp)){
    stop("Cannot supply both qtl.df and tag.snp.")
  } else if(is.null(qtl.df) & is.null(tag.snp)){
    stop("Must supply one of qtl.df or tag.snp.")
  # Use tag snp
  } else if(!is.null(tag.snp)){
    this.snp.name <- tag.snp
    key.snp <- tag.snp
    make_ld(plink.path = plink.exec, 
            snp.name = this.snp.name, 
            window = window, 
            bedfile = geno.bed, 
            in.dir = in.dir,
            out.dir = out.dir)
    ld.table <- read.delim(file.path(out.dir, "ld_out_temp.vcor"), header = T)
    ld.table_out <- ld.table %>%
      select("marker.ID" = "ID_B", "R2" = "UNPHASED_R2")  
  } else {
    this.clump.df <- qtl.df %>%
      mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-"))
    
    # Might be a tie here sometimes if top snp has two traits with same pvalue, I think the [1] solves it
    key.snp <- this.clump.df$marker.ID[which.max(this.clump.df$LOGPVAL)[1]]
    # initialize loop
    this.snp.name <- key.snp
    # make empty storage frame
    ld.table_all <- data.frame()
    # decide if we want a progress bar
    do.progress <- nrow(this.clump.df) > 1 & verbose
    
    # send some nice messages
    if(verbose){
      message(paste("Calculating LD of", nrow(this.clump.df), "snps to all snps in window."))
    }
    if(do.progress){
      pb <- txtProgressBar(min = 1, max = nrow(this.clump.df), style = 3)
    }
    # start looping through rows of qtl.df
    for(i in 1:nrow(this.clump.df)){
      
      this.snp.name <- this.clump.df$marker.ID[i]
      
      make_ld(plink.path = plink.exec, 
              snp.name = this.snp.name, 
              window = window, 
              bedfile = geno.bed, 
              in.dir = in.dir,
              out.dir = out.dir)
      ld.table <- read.delim(file.path(out.dir, "ld_out_temp.vcor"), header = T)
      ld.table_sub <- ld.table %>%
        select("marker.ID" = "ID_B", "R2" = "UNPHASED_R2")  
      ld.table_all <- bind_rows(ld.table_all, ld.table_sub)
      if(do.progress){
        setTxtProgressBar(pb, i)
      }
    }
    if(do.progress){
      close(pb)
    }
    # filter for top r2 for each marker in the window
    ld.table_out <- ld.table_all %>%
      group_by(.data$marker.ID) %>%
      arrange(desc(.data$R2)) %>%
      filter(1:n() == 1)
  }
  
  # retain list of qtl snps 
  if(!is.null(tag.snp)){
    qtl.snps <- NULL
  } else {
    qtl.snps <- unique(this.clump.df$marker.ID)
  }
  
  return(list(table = ld.table_out,
              key.snp = key.snp,
              qtl.snps = qtl.snps))
}