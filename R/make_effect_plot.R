# ------------------------------------------------------------------------\
# main function --------
# ------------------------------------------------------------------------\

#' Make volcano style effect size vs pvalue plot
#'
#' @param panvar.table.list list, output from [panvaR::make_panvar_tables]. Provide either this list or both gwas.res and ld.list.
#' @param gwas.res data.frame of all gwas results, should contain columns (CHR, POS, PVAL), corresponding to (chromosome, physical position, and pvalue).
#' @param ld.list list, output of [panvaR::get_ld_in_window]
#' @param pvals.in.log boolean, are pvalues in input data.frames in -log10(p)?
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD
#' @param unplotted.alpha numeric, number from 0 to 1 to indicate alpha values of snps below the plot.r2.thresh. 
#' To not plot these snps set value to 0. 
#' @param window numeric, kilobases on either side of top QTL snp to plot
#' @param sig.line numeric, -log10(p) value to draw line on plot
#' @param orient character, will rotate plot 90 degrees. vertical (V) or horizontal (H)
#' refers to how the "buildings" of the plot are plotted. 
#' "V" places pvalue on y-axis, "H" places pvalues on x-axis. 
#' @param qualitative.annotation character, column in gwas.res that contains qualitative annotations.
#' For example impact grades from snpeff. See [panvaR::format_snpeff_annotations].
#' Will be plotted as shapes. Only accepts up to 5 classes. "IMPACT" and "IMPACT_PLUS" are special 
#' cases that will have a pre-assigned scale used if supplied here.
#' @param qualitative.shape.scale ggplot scale, an object with a stored call to 
#' [ggplot2::scale_shape_manual]. More often an output of the function [panvaR::make_consistent_scale]. 
#' @param quantitative.annotation character, column in gwas.res that contains quantitative annotations. 
#' For example, variant effect scores. Will be plotted as fill to points. 
#' @param quantitative.fill.scale character or scale object, either a character indicating the
#' `option` parameter passed to [ggplot2::scale_fill_viridis_b] that alters the color scale used.
#' Or a previous call to a ggplot2 fill scale for example [ggplot2::scale_fill_stepsn].
#' @param include.legend boolean, if TRUE, legend will be included. 
#'
#' @returns
#' GGplot object of plot. Points colored by maximum R2 to snps in qtl.df
#' @export
#'
#' @examples
#' # Work in progress
make_effect_plot <- function(panvar.table.list = NULL,
                             gwas.res = NULL,
                             ld.list = NULL,
                             pvals.in.log = TRUE,
                             plot.r2.thresh = .2,
                             unplotted.alpha = .4,
                             window,
                             sig.line,
                             orient = c("V", "H"),
                             qualitative.annotation = NULL,
                             qualitative.shape.scale = NULL,
                             quantitative.annotation = NULL,
                             quantitative.fill.scale = NULL,
                             include.legend = T){
  
  # make sure we don't use both
  if(!is.null(panvar.table.list) & (!is.null(ld.list) | !is.null(gwas.res))){
    stop("Must provide only one of panvar.table.list or ld.list and gwas.res.")
  }
  
  # use ld.list and gwas.res
  if(is.null(panvar.table.list)){
    if(is.null(ld.list) & is.null(gwas.res)){
      stop("Must provide either panvar.table.list or ld.list and gwas.res.")
    } else if(is.null(ld.list) | is.null(gwas.res)){
      stop("Must provide ld.list and gwas.res.")
    } else {
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
      
      # add middlesnp back in
      gwas.sub.mid.snp <- gwas.res %>% 
        filter(.data$marker.ID == ld.list$key.snp) %>% 
        mutate(LD = 1)
      gwas.sub <- bind_rows(gwas.sub, gwas.sub.mid.snp)
    }
    key.snp <- ld.list$key.snp
    
    # use panvar.table.list
  } else {
    key.snp <- panvar.table.list$key.snp
    gwas.sub <- panvar.table.list$gwas %>% 
      mutate(R2 = .data$LD)
    
    # for rug plot
    marker.list.in.window <- gwas.sub %>%
      mutate(POS = get_bp_from_id(.data$marker.ID)) %>%
      select("marker.ID", "POS") %>%
      distinct()
  }
  
  
  # make manhattan
  plot.df <- gwas.sub %>%
    # alpha scale
    mutate(how.to.plot = case_when(.data$R2 > plot.r2.thresh ~ 1,
                                   TRUE ~ unplotted.alpha)) %>%
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
  
  # plot only points in window
  this.pos <- get_bp_from_id(key.snp)
  plot.df <- plot.df %>% 
    filter(between(.data$POS, this.pos - window * 1000, this.pos + window * 1000)) 
  # plot.limits <- c(this.pos + window * 1000, this.pos - window * 1000)
  # plot.limits.ex <- c(plot.limits[1] + y.spread.factor.window, plot.limits[2] - y.spread.factor.window)
  
  
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
  # prepare quantitative variable --------
  # ------------------------------------------------------------------------\
  
  # check if the scale provided is a ggplot scale or not
  if("Scale" %in% class(quantitative.fill.scale)){
    user.supplied.scale.option <- "scale"
  } else {
    user.supplied.scale.option <- "not.scale"
  }
  
  if(!is.null(quantitative.annotation) & user.supplied.scale.option == "scale"){
    # Both variable and scale provided 
    # use user supplied scale
    quantitative.fill.scale <- quantitative.fill.scale
    # only plot that stuff that's above the LD threshold 
    plot.df <- plot.df %>% 
      mutate(plot.quant.var = case_when(.data$R2 < plot.r2.thresh ~ NA,
                                        TRUE ~ .data[[quantitative.annotation]]))
  } else if(!is.null(quantitative.annotation)){
    # Only variable provided 
    # set a binned color scale
    if(is.null(quantitative.fill.scale)){
      quantitative.fill.scale <- scale_fill_viridis_b(name = quantitative.annotation)
    } else {
      quantitative.fill.scale <- scale_fill_viridis_b(name = quantitative.annotation, option = quantitative.fill.scale)
    }
    # only plot that stuff that's above the LD threshold 
    plot.df <- plot.df %>% 
      mutate(plot.quant.var = case_when(.data$R2 < plot.r2.thresh ~ NA,
                                        TRUE ~ .data[[quantitative.annotation]]))
  } else if(is.null(quantitative.annotation) & !is.null(quantitative.fill.scale)){
    warning("No quantitative annotation specified. Provided quantitative.fill.scale ignored.")
  } else {
    # do nothing
  }
  
  # ------------------------------------------------------------------------\
  # base plot --------
  # ------------------------------------------------------------------------\
  
  effect.plot <-
    ggplot(aes(x = .data$EFF, y = .data$PVAL), data = plot.df) +
    theme_bw() +
    # geom_point(aes(color = .data$plot.R2, alpha = .data$how.to.plot), shape = 16, show.legend = include.legend) +
    # scale_color_stepsn(colors = my.colors, name = "R2") +
    geom_hline(yintercept = sig.line) +
    theme(panel.grid = element_blank()) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray')
  
  # ------------------------------------------------------------------------\
  # add extra annotations to plot --------
  # ------------------------------------------------------------------------\
  
  # set default LD color scale
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
  
  # Some logic for including extra variables
  if(!is.null(qualitative.annotation) & !is.null(quantitative.annotation)){
    # Use a quant and a qual
    effect.plot <- effect.plot +
      geom_point(aes(fill = .data$plot.quant.var, alpha = .data$how.to.plot, shape = .data[[qualitative.annotation]]),
                 size = 3, color = "black", show.legend = include.legend) +
      scale_alpha(guide = "none", range = c(0, 1)) +
      quantitative.fill.scale +
      qualitative.shape.scale
    
  } else if(is.null(qualitative.annotation) & !is.null(quantitative.annotation)){
    # Use just a quant 
    effect.plot <- effect.plot +
      geom_point(aes(fill = .data$plot.quant.var, alpha = .data$how.to.plot), 
                 shape = 21, size = 3, color = "black", show.legend = include.legend) +
      scale_alpha(guide = "none", range = c(0, 1)) +
      quantitative.fill.scale
    
  } else if(!is.null(qualitative.annotation) & is.null(quantitative.annotation)){
    # Use just a qual
    effect.plot <- effect.plot + 
      geom_point(aes(fill = .data$plot.R2, alpha = .data$how.to.plot, shape = .data[[qualitative.annotation]]), 
                 size = 3, color = "black", show.legend = include.legend) +
      qualitative.shape.scale +
      scale_alpha(guide = "none", range = c(0, 1)) +
      default.LD.fill.scale
    
  } else {
    # Use neither
    effect.plot <- effect.plot + 
      geom_point(aes(fill = .data$plot.R2, alpha = .data$how.to.plot), 
                 size = 3, shape = 21, color = "black", show.legend = include.legend) +
      scale_alpha(guide = "none", range = c(0, 1)) +
      default.LD.fill.scale
  }
  
  # flip it if you want
  orient <- match.arg(orient)
  if(orient == "V"){
    effect.plot <- effect.plot +
      scale_x_continuous(
        # limits = plot.limits.ex,
        labels = scales::label_number(scale_cut = scales::cut_short_scale())
      ) +
      # bquote doesn't work with plotly
      # labs(y = bquote(-log[10](p-value))) +
      labs(y = "-log(p-value)") +
      theme(axis.title.x = element_blank())
  } else if(orient == "H"){
    effect.plot <- effect.plot +
      coord_flip() +
      # scale_y_reverse() +
      # scale_x_reverse(
      #   limits = plot.limits.ex,
      #   labels = scales::label_number(scale_cut = scales::cut_short_scale())
      # ) +
      # bquote doesn't work with plotly
      # labs(y = bquote(-log[10](p-value))) 
      labs(y = "-log(p-value)") +
      theme(axis.title.y = element_blank())
  } else {
    stop("Orientation must be either vertical (V) or horizontal (H).")
  }
  
  return(effect.plot)
  
}

