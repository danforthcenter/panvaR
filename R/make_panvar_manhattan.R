#' Make sideways manhattan plot for building locus zoom. Receives output from a single gwas model.
#'
#' @param gwas.res table of all gwas results, should contain columns (CHR, POS, PVAL), corresponding to (chromosome, physical position, and pvalue).
#' @param qtl.df table with same columns that includes only significant hits in a qtl. QTL are typically defined as hits grouped by LD by something like `plink --clump`
#' @param pvals.in.log boolean, are pvalues in input data.frames in -log10(p)?
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD
#' @param ld.list output of [luebbert::get_ld_in_window]
#' @param window kilobases on either side of top QTL snp to plot
#' @param sig.line -log10(p) value to draw line on plot
#'
#' @returns
#' GGplot of manhattan plot with points colored by maximum R2 to markers in the qtl.df.
#' @export
#'
#' @examples
#' # Work in progress
make_panvar_manhattan <- function(gwas.res,
                                    qtl.df = NULL,
                                    pvals.in.log = TRUE,
                                    plot.r2.thresh = .2,
                                    ld.list,
                                    window,
                                    sig.line)
  # qualitative.annotation (shape)
  # quantitative.annotation (color) = LD
  # quantitative.fill.scall (from its own function)
  # qualititave.shape.scale (from its own function)
  {
  
  
  # chrom <- unique(qtl.df$CHR)
  
  gwas.sub <- gwas.res %>%
    # filter(.data$CHR == chrom) %>%
    mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
    # filter(between(physical.pos, this.pos - window * 1000, this.pos + window * 1000)) %>%
    left_join(ld.list$table, by = "marker.ID") %>%
    filter(!is.na(.data$R2))
  # for rug plot
  marker.list.in.window <- ld.list$table %>%
    mutate(POS = get_bp_from_id(.data$marker.ID)) %>%
    select("marker.ID", "POS") %>%
    distinct()
  
  
  # make manhattan
  plot.df <- gwas.sub %>%
    # alpha scale
    mutate(how.to.plot = case_when(.data$R2 > plot.r2.thresh ~ 1,
                                   TRUE ~ .4)) %>%
    # color scale
    mutate(plot.R2 = case_when(.data$R2 < plot.r2.thresh ~ NA,
                               TRUE ~ R2))
  
  # change pvalue if needed
  if(!pvals.in.log){
    plot.df <- plot.df %>%
      mutate(PVAL = -log10(.data$PVAL))
  }
  
  # how far to spread labels past ends, in percentage
  y.spread.expansion <- .1
  y.spread.factors <- c(1 + y.spread.expansion, 1 - y.spread.expansion)
  y.spread.factor.window <- (window * y.spread.expansion) * 1000
  
  # plot limits
  this.pos <- get_bp_from_id(ld.list$key.snp)
  plot.limits <- c(this.pos + window * 1000, this.pos - window * 1000)
  plot.limits.ex <- c(plot.limits[1] + y.spread.factor.window, plot.limits[2] - y.spread.factor.window)
  
  # Base plot
  man <-
    ggplot(aes(x = .data$POS, y = .data$PVAL), data = plot.df) +
    geom_point(aes(fill = .data$plot.R2, alpha = .data$how.to.plot), size = 3, shape = 21, color = "black") +
    scale_alpha(guide = "none") +
    binned_scale(aesthetics = "fill",
                 name = "R2 \n",
                 palette = function(x) c("#43638E", "#88DAA0", "#DBC32D", "#B94712"),
                 limits = c(plot.r2.thresh, 1),
                 breaks = seq(plot.r2.thresh, 1, length.out = 5)[-c(1,5)],
                 show.limits = T,
                 guide = "colorsteps",
                 na.value = "grey50") +
    # shape.scale +
    scale_x_reverse(
      limits = plot.limits.ex,
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      name = "-log10pval"
    ) +
    coord_flip() +
    scale_y_reverse() +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "grey95"),
      axis.title.y = element_blank(),
      legend.position = "left",
      legend.justification = "top",
      panel.grid = element_blank()
    ) +
    labs(y = bquote(-log[10](p - value))) +
    geom_rug(aes(x = .data$POS),
             sides = 't',
             data = marker.list.in.window,
             inherit.aes = F) +
    geom_hline(yintercept = sig.line, linetype = 'dashed')
  
  return(man)
}
