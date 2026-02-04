# ------------------------------------------------------------------------\
# unexported --------
# ------------------------------------------------------------------------\

get.gene.from.snp <- function(bp, gene.df){
  check.vec <- data.table::between(bp, gene.df$start, gene.df$end)
  if(any(check.vec)){
    gene.id.out <- list(gene.df$geneID[check.vec])
    # gene.id.out <- paste(gene.df$geneID[check.vec], collapse = "|")
  } else {
    gene.id.out <- NA
  }
  return(gene.id.out)
}

# ------------------------------------------------------------------------\
# main function --------
# ------------------------------------------------------------------------\


#' Make a locus zoom style plot with results from a single gwas, points colored by R2 to a qtl and nearby genes annotated.
#'
#'
#' @param gwas.res table of all gwas results, should contain columns (`CHR`, `POS`, `PVAL`), corresponding to (chromosome, physical position, and pvalue).
#' @param qtl.df table with same columns that includes only significant hits in a qtl. QTL are typically defined as hits grouped by LD by something like `plink --clump`.
#' @param annotation.table table with annotations with columns (geneID, CHR, start, end, annotation). start and end correspond to base-pair coordinates of start and end of gene. CHR is chromosome of gene.
#' @param pvals.in.log boolean, are pvalues in input data.frames in -log10(p) format?
#' @param plink.path path to plink 1.9 executable.
#' @param geno.bed prefix of genotype files in plink (bed/bim/fam) format.
#' @param geno.bed.filename prefix of genotype files in plink (bed/bim/fam) format.
#' @param geno.bed.directory directory where genotype files are located.
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD. SNPs below are plotted in grey.
#' @param window kilobases on either side of top QTL snp to plot.
#' @param sig.line -log10(p) value to draw line on plot.
#' @param plot.title Optional. Title of plot.
#' @param include.gene.id boolean, include geneID in gene annotations or not
#' @param plot.effect boolean, include plot of pvalues vs effect size? If TRUE, `gwas.res` and `qtl.df` must have column `EFF` that conatins gwas effect sizes.
#'
#' @returns
#' ggplot object of plot
#' @export
#'
#' @examples
#' # work in progress
make_panvar_plot <- function(gwas.res,
                             qtl.df = NULL,
                             tag.snp = NULL,
                             annotation.table,
                             plink.path,
                             pvals.in.log = T,
                             geno.bed.filename,
                             geno.bed.directory = "/.",
                             temp.dir = tempdir(),
                             plot.r2.thresh = .2,
                             unplotted.alpha = .4,
                             window,
                             sig.line,
                             orient = "H",
                             qualitative.annotation = NULL,
                             qualitative.shape.scale = NULL,
                             quantitative.annotation = NULL,
                             quantitative.fill.scale = NULL,
                             plot.title = "",
                             include.gene.id = F,
                             highlight.gene.ids = NULL,
                             gene.highlight.color = "red",
                             annotation.point.variable = "LD",
                             annotation.point.scale = NULL,
                             plot.effect = F) {
  
  # ------------------------------------------------------------------------\
  # make LD --------
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
  # make manhattan --------
  # ------------------------------------------------------------------------\
    
  message("Making manhattan")
  man <- make_panvar_manhattan(gwas.res = gwas.res,
                               pvals.in.log = pvals.in.log,
                               plot.r2.thresh = plot.r2.thresh,
                               unplotted.alpha = unplotted.alpha,
                               ld.list = ld.list,
                               window = window,
                               sig.line = sig.line,
                               orient = orient,
                               qualitative.annotation = qualitative.annotation,
                               qualitative.shape.scale = qualitative.shape.scale,
                               quantitative.annotation = quantitative.annotation,
                               quantitative.fill.scale = quantitative.fill.scale)
  
  # ------------------------------------------------------------------------\
  # make snp to gene stats --------
  # ------------------------------------------------------------------------\
  
  if(!is.null(annotation.point.variable)){
    message("Generating snp to gene correspondence")
    
    # filter gwas df to just in window and join LD
    gwas.sub <- gwas.res %>%
      as.data.frame() %>% 
      mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
      left_join(ld.list$table, by = "marker.ID") %>%
      filter(!is.na(.data$R2)) %>% 
      rename("LD" = "R2")
    
    # filter anno to just window
    this.chrom <- get_chrom_from_id(middle.snp)
    this.pos <- get_bp_from_id(middle.snp)
    anno.sub <- annotation.table %>%
      filter(.data$CHR == this.chrom) %>%
      rowwise() %>%
      mutate(dist.from.snp = get.gene.dist.from.snp(this.pos, .data$start, .data$end)) %>%
      filter(.data$dist.from.snp <= window * 1000) 
    
    # get maximum value per gene 
    point.color.stat <- gwas.sub %>% 
      rowwise() %>% 
      mutate(snp.in.gene = get.gene.from.snp(.data$POS, anno.sub)) %>% 
      filter(!is.null(.data$snp.in.gene)) %>% 
      unnest_longer(.data$snp.in.gene) %>% 
      group_by(.data$snp.in.gene) %>% 
      summarize(maximum.value = max(.data[[annotation.point.variable]])) %>% 
      rename("geneID" = "snp.in.gene")
    
    anno.in <- annotation.table %>% 
      left_join(point.color.stat, by = "geneID")
  } else {
    anno.in <- annotation.table
  }
  
  

  # ------------------------------------------------------------------------\
  # make annotation --------
  # ------------------------------------------------------------------------\
  
  if(is.null(annotation.point.variable)){
    point.color.option <- NULL
  } else if(annotation.point.variable == "LD"){
    default.LD.fill.scale <- binned_scale(
      aesthetics = "fill",
      name = "R2 \n",
      palette = function(x)
        c("#43638E", "#88DAA0", "#DBC32D", "#B94712"),
      limits = c(plot.r2.thresh, 1),
      breaks = seq(plot.r2.thresh, 1, length.out = 5)[-c(1, 5)],
      show.limits = T,
      guide = "colorsteps",
      na.value = "grey50"
    )
    annotation.point.scale <- default.LD.fill.scale
    anno.in <- anno.in %>% 
      rename("LD" = "maximum.value")
    point.color.option <- "LD"
  } else if(!is.null(annotation.point.variable)) {
    point.color.option <- maximum.value
  } else {
    point.color.option <- NULL
  }

  anno <- make_gene_annotation_plot(annotation.table = anno.in,
                                    middle.snp = ld.list$key.snp,
                                    window = window,
                                    include.id = include.gene.id,
                                    highlight.ids = highlight.gene.ids,
                                    highlight.color = gene.highlight.color,
                                    use.arrows = F,
                                    point.color = point.color.option,
                                    point.fill.scale = annotation.point.scale)

  # ------------------------------------------------------------------------\
  # make effect --------
  # ------------------------------------------------------------------------\
    
  if (plot.effect) {
    effect.plot <- make_effect_plot(gwas.res,
                                    qtl.df,
                                    pvals.in.log,
                                    plot.r2.thresh,
                                    ld.list,
                                    window,
                                    sig.line)
    
  }
  
  # ------------------------------------------------------------------------\
  # final plot --------
  # ------------------------------------------------------------------------\
  
  chrom <- unique(qtl.df$CHR)
  
  if (plot.effect) {
    out <-
      man + theme(plot.margin = ggplot2::margin(6, 6, 6, 120, unit = "points")) +
      patchwork::inset_element(
        effect.plot,
        left = .0,
        bottom = 0,
        top = .3,
        right = .3,
        align_to = "full"
      ) +
      anno +
      patchwork::plot_layout(nrow = 1, widths = c(4, 6)) +
      patchwork::plot_annotation(title = plot.title,
                                 subtitle = ld.list$key.snp)
    
  } else {
    # # specify manhattan margin (smaller if there isn't an effect plot)
    # man <- man +
    #   theme(plot.margin = ggplot2::margin(6, 6, 6, 6, unit = "points"))
    # out <- 
    #   patchwork::wrap_plots(man, anno, nrow = 1, widths = c(4, 6))
    
    out <-
      man + theme(plot.margin = ggplot2::margin(6, 6, 6, 6, unit = "points")) +
      anno +
      patchwork::plot_layout(nrow = 1, widths = c(4, 6)) +
      patchwork::plot_annotation(title = plot.title,
                                 subtitle = ld.list$key.snp)
  }
  
  return(out)
}
