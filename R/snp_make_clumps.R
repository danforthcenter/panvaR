#' Assign snps to clumps based on LD after gwas
#'
#' @param geno.bed.filename character, filename of genotype bedfile, no .bed extension
#' @param geno.bed.dir character, path to directory where bed/bim/fam files exist, include trailing "/"
#' @param gwas.res data.frame, table of gwas results with columns (marker.ID, CHR, POS, PVAL)
#' @param pvals.in.log boolean, are pvalues in `gwas.res` in -log10(p) or not
#' @param window integer, window in kilobases either side of snp to look for snps in LD,
#' @param ld.thresh numeric, R2 threshold above which snps will be grouped
#' @param plink.path path to plink 2 executable
#' @param out.dir character, path to a temporary directory to output some intermediate files
#'
#' @returns
#' table with columns marker.ID and clump_num. Clump_num indicates groupings, numberings start from the larges pvalue to smallest. May want to reassign afterwards to be along the genome.
#' @export
#'
#' @examples
#' # work in progress
snp_make_clumps <- function(geno.bed.filename,
                        geno.bed.dir,
                        gwas.res,
                        pvals.in.log = T,
                        window = 500,
                        ld.thresh = .5,
                        plink.path = NULL,
                        out.dir = NULL) {
  
  # ~~~~ Initialize ~~~~
  
  # check plink.path
  if(!is.null(plink.path)){
    plink.exec <- plink.path
  } else if(!is.null(options()$plink_path)){
    plink.exec <- options()$plink_path
  } else {
    plink.exec <- Sys.which("plink2")
  }
  
  # check output directory
  out.dir <- temporary_directory(out.dir)
  
  # make wellformed path
  geno.bed.dir <- normalizePath(geno.bed.dir)
  
  # check for .bed extension on input
  if(grepl("\\.bed$", geno.bed.filename)){
    geno.bed.filename <- tools::file_path_sans_ext(geno.bed.filename)
  }
  
  # check that bedfile exists
  if(!file.exists(file.path(geno.bed.dir, paste0(geno.bed.filename, ".bed")))){
    stop("Bed file not found in input directory.")
  }
  
  # convert pvals to log if we need
  if(!pvals.in.log){
    gwas.res <- gwas.res %>%
      mutate(PVAL = -log10(.data$PVAL))
  }
  
  # get window in basepairs from kilobases
  window.bp <- window * 1000
  
  # get snps sorted by logpvalue, highest to lowest
  snps.to.test <- gwas.res %>%
    arrange(desc(.data$PVAL)) %>%
    pull(.data$marker.ID) %>%
    unique()
  
  # make an output df and start a counter
  out <- data.frame()
  i <- 1
  # start a progress bar
  message("Creating clumps...")
  pb.end.value <- length(snps.to.test)
  pb <- txtProgressBar(min = 1, max = pb.end.value, style = 3)
  while(length(snps.to.test) > 0){
    
    this.snp.name <- snps.to.test[1]
    # check if any snps are in window
    this.snp.info <- filter(gwas.res, .data$marker.ID == this.snp.name) %>%
      select("SNP" = "marker.ID", "POS", "CHR") %>%
      distinct()
    
    pos.range <- c(this.snp.info$POS - window.bp, this.snp.info$POS + window.bp)
    
    snps.in.range <- gwas.res %>%
      filter(.data$CHR == this.snp.info$CHR) %>%
      filter(between(.data$POS, pos.range[1], pos.range[2]))
    
    if(nrow(snps.in.range) < 2){
      this.out <- data.frame(marker.ID = this.snp.info$SNP,
                             clump_num = i)
      out <- bind_rows(out, this.out)
      snps.to.test <- snps.to.test[-which(snps.to.test %in% this.out$marker.ID)]
      i <- i+1
      setTxtProgressBar(pb, pb.end.value - length(snps.to.test))
      
    } else {
      
      make_ld(plink.exec,
              this.snp.name,
              window,
              geno.bed.filename,
              geno.bed.dir,
              out.dir)
      
      ld.table <- read.delim(file.path(out.dir, "ld_out_temp.vcor"), header = T)
      ld.table_sub <- ld.table %>%
        select("SNP" = "ID_B", "R2" = "UNPHASED_R2") %>%
        filter(.data$R2 >= ld.thresh) %>%
        # plink 2 does not retain key snp
        add_row(SNP = this.snp.name, R2 = 1) %>% 
        filter(.data$SNP %in% snps.to.test)
      
      # ld.table <- read.delim(file.path(out.dir, "ld_out_temp.vcor"), header = T)
      # ld.table_sub <- ld.table %>%
      #   select("marker.ID" = "ID_B", "R2" = "UNPHASED_R2")  
      
      
      this.out <- data.frame(marker.ID = ld.table_sub$SNP,
                             clump_num = i)
      out <- bind_rows(out, this.out)
      snps.to.test <- snps.to.test[-which(snps.to.test %in% this.out$marker.ID)]
      i <- i+1
      # make this into a progress bar
      # message(paste(length(snps.to.test), "snps remaining."))
      setTxtProgressBar(pb, pb.end.value - length(snps.to.test))
    }
  }
  close(pb)
  
  return(out)
}



