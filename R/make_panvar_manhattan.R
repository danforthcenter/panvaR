#' Make sideways manhattan plot for building locus zoom. Receives output from a single gwas model.
#'
#' @param gwas.res table of all gwas results, should contain columns (CHR, POS, PVAL), corresponding to (chromosome, physical position, and pvalue).
#' @param qtl.df table with same columns that includes only significant hits in a qtl. QTL are typically defined as hits grouped by LD by something like `plink --clump`
#' @param pvals.in.log boolean, are pvalues in input data.frames in -log10(p)?
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD
#' @param ld.list output of [luebbert::get_ld_in_window]
#' @param window kilobases on either side of top QTL snp to plot
#' @param sig.line -log10(p) value to draw line on plot
#' @param orient character, will rotate plot 90 degrees. vertical (V) or horizontal (H)
#' refers to how the "buildings" of the plot are plotted. 
#' @param qualitative.annotation character, column in gwas.res that contains qualitative annotations.
#' For example impact grades from snpeff. See [panvaR::format_snpeff_annotations].
#' Will be plotted as shapes. Only accepts up to 5 classes. 
#' @param qualitative.shape.scale ggplot scale, an object with a stored call to 
#' [ggplot2::scale_shape_manual]. More often an output of the function [panvaR::make_consistent_scale]. 
#'
#' @returns
#' GGplot of manhattan plot with points colored by R2.
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
                                  sig.line,
                                  orient = c("H", "V"),
                                  qualitative.annotation = NULL,
                                  qualitative.shape.scale = NULL)
  # qualitative.annotation (shape)
  # quantitative.annotation (color) = LD
  # quantitative.fill.scall (from its own function)
  # qualititave.shape.scale (from its own function)
  {
  
  gwas.sub <- gwas.res %>%
    as.data.frame() %>% 
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
  

  # ------------------------------------------------------------------------\
  # prepare qualitative variable --------
  # ------------------------------------------------------------------------\
  
  
  # If we're using IMPACT then use this scale
  if(is.null(qualitative.shape.scale) & is.null(qualitative.annotation)){
    # do nothing
  } else if(is.null(qualitative.shape.scale) & qualitative.annotation == "IMPACT"){
    qual.vars <-  c("HIGH", "MODERATE", "LOW", "MODIFIER")
    qualitative.shape.scale <- make_consistent_scale(values = c(24, 22, 25, 23),
                                                     vars = qual.vars,
                                                     type = "shape",
                                                     name = qualitative.annotation)
    
    plot.df$IMPACT <- factor(plot.df$IMPACT, levels = qual.vars)
  } else if(is.null(qualitative.shape.scale) & qualitative.annotation == "IMPACT_PLUS"){
    qual.vars <-  c("HIGH", "MODERATE", "LOW", "MODIFIER_CODING", "MODIFIER_INTERGENIC")
    qualitative.shape.scale <- make_consistent_scale(values = c(24, 22, 25, 23, 21),
                                                     vars = qual.vars,
                                                     type = "shape",
                                                     name = qualitative.annotation)
    
    plot.df$IMPACT_PLUS <- factor(plot.df$IMPACT_PLUS, levels = qual.vars)
    
  } else if(is.null(qualitative.shape.scale)){
    # make one using the unique values of qualitative annotation column
    qual.vars <- sort(unique(plot.df[,qualitative.annotation]))
    
    # throw an error if there are more than 5 things for now
    if(qual.vars > 5){
      stop("Qualitative variable must not have more than 5 unique values. 
           Points are plotted using the fill aesthetic to give them a color, 
           R only has 5 shapes (21-25) that can be assigned fill values.")
    }
    qualitative.shape.scale <- make_consistent_scale(values = c(21:25)[1:length(qual.vars)],
                                                     vars = qual.vars,
                                                     type = "shape",
                                                     name = qualitative.annotation)
  }
  
  # ------------------------------------------------------------------------\
  # base plot --------
  # ------------------------------------------------------------------------\
  
  man <-
    ggplot(aes(x = .data$POS, y = .data$PVAL), data = plot.df) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "grey95"),
      legend.position = "left",
      legend.justification = "top",
      panel.grid = element_blank()
    ) +
    geom_hline(yintercept = sig.line, linetype = 'dashed')
  
  # ------------------------------------------------------------------------\
  # add extra annotations to plot --------
  # ------------------------------------------------------------------------\
  
  
  if(is.null(qualitative.annotation)){
    man <- man + 
      geom_point(aes(fill = .data$plot.R2, alpha = .data$how.to.plot), size = 3, shape = 21, color = "black") +
      scale_alpha(guide = "none") +
      binned_scale(aesthetics = "fill",
                   name = "R2 \n",
                   palette = function(x) c("#43638E", "#88DAA0", "#DBC32D", "#B94712"),
                   limits = c(plot.r2.thresh, 1),
                   breaks = seq(plot.r2.thresh, 1, length.out = 5)[-c(1,5)],
                   show.limits = T,
                   guide = "colorsteps",
                   na.value = "grey50") 
  } else {
    man <- man + 
      geom_point(aes(fill = .data$plot.R2, alpha = .data$how.to.plot, shape = .data[[qualitative.annotation]]), size = 3, color = "black") +
      qualitative.shape.scale +
      scale_alpha(guide = "none") +
      binned_scale(aesthetics = "fill",
                   name = "R2 \n",
                   palette = function(x) c("#43638E", "#88DAA0", "#DBC32D", "#B94712"),
                   limits = c(plot.r2.thresh, 1),
                   breaks = seq(plot.r2.thresh, 1, length.out = 5)[-c(1,5)],
                   show.limits = T,
                   guide = "colorsteps",
                   na.value = "grey50") 
  }
  
  # add a quantitative variable
  
  
  # flip it if you want
  orient <- match.arg(orient)
  if(orient == "V"){
    man <- man +
      geom_rug(aes(x = .data$POS),
               sides = 'b',
               data = marker.list.in.window,
               inherit.aes = F) +
      scale_x_continuous(
        limits = plot.limits.ex,
        labels = scales::label_number(scale_cut = scales::cut_short_scale())
      ) +
      # bquote doesn't work with plotly
      # labs(y = bquote(-log[10](p-value))) +
      labs(y = "-log(p-value)") +
      theme(axis.title.x = element_blank())
  } else if(orient == "H"){
    man <- man +
      coord_flip() +
      scale_y_reverse() +
      geom_rug(aes(x = .data$POS),
               sides = 't',
               data = marker.list.in.window,
               inherit.aes = F) +
      scale_x_reverse(
        limits = plot.limits.ex,
        labels = scales::label_number(scale_cut = scales::cut_short_scale())
      ) +
      # bquote doesn't work with plotly
      # labs(y = bquote(-log[10](p-value))) 
      labs(y = "-log(p-value)") +
      theme(axis.title.y = element_blank())
  } else {
    stop("Orientation must be either vertical (V) or horizontal (H).")
  }

  return(man)
}
