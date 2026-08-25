#' Make sideways manhattan plot for building locus zoom. Receives output from a single gwas model.
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
#' @param point.shape.variable character, column in gwas.res that contains qualitative annotations to be mapped to point shapes.
#' For example impact grades from snpeff. See [panvaR::format_snpeff_annotations].
#' Only accepts up to 5 classes. "IMPACT" and "IMPACT_PLUS" are special 
#' cases that will have a pre-assigned scale used if supplied here.
#' @param point.shape.scale ggplot scale, an object with a stored call to 
#' [ggplot2::scale_shape_manual]. More often an output of the function [panvaR::make_consistent_scale]. 
#' @param point.fill.variable.c character, column in gwas.res that contains quantitative annotations to be plotted
#' as a continuous variable mapped to point fill. For example, variant effect scores. Only provide either continuous 
#' or discrete quantitative annotations.
#' @param point.fill.scale.c character or scale object, either a character indicating the
#' `option` parameter passed to [ggplot2::scale_fill_viridis_b] that alters the color scale used.
#' Or a previous call to a ggplot2 continuous fill scale for example [ggplot2::scale_fill_stepsn].
#' @param point.fill.variable.d character, column in gwas.res that contains annotations to be plotted as a discrete 
#' variable mapped to point fill. For example, Year or Trial if combining multiple gwas results. 
#' Only provide either continuous or discrete quantitative annotations.
#' @param point.fill.scale.d character or scale object, either a character indicating the
#' `option` parameter passed to [ggplot2::scale_fill_viridis_d] that alters the color scale used.
#' Or a previous call to a ggplot2 discrete fill scale for example [ggplot2::scale_fill_discrete].
#'
#' @returns
#' GGplot of manhattan plot with points colored by R2. Accepts input
#' @export
#'
#' @examples
#' # organize options
#' tag.snp <- "Chr_05-6857045"
#' gwas.df <- read.csv(system.file(
#'     "extdata",
#'     "PanvarExample_GLM_GWASresults.csv",
#'     package = "panvaR"))
#' annotation.table <- read.csv(system.file(
#'     "extdata",
#'     "Setaria_shattering_annotation.csv",
#'     package = "panvaR"))
#' plink.path <- bigsnpr::download_plink2()
#' temp.dir <- file.path(tempdir(), "panvar_ex")
#' dir.create(temp.dir, showWarnings = FALSE)
#' geno.bed.filename <- "Setaria_shattering_example_pruned.bed"
#' geno.bed.directory <- system.file("extdata", package="panvaR")
#'
#' # make input tables
#' tables <- make_panvar_tables(
#'   gwas.res = gwas.df,
#'   tag.snp = tag.snp,
#'   annotation.table = annotation.table,
#'   plink.path = plink.path,
#'   pvals.in.log = F,
#'   geno.bed.filename = geno.bed.filename,
#'   geno.bed.directory = geno.bed.directory,
#'   window = 25,
#'   temp.dir = temp.dir,
#'   compute.scores = FALSE,
#'   snp.to.gene.buffer = 0)
#'   
#' # make plot
#' plot_panvar_manhattan(
#'   panvar.table.list = tables,
#'   pvals.in.log = FALSE,
#'   window = 25,
#'   sig.line = 6)
#' 
#' # flip it vertical if you want
#' plot_panvar_manhattan(
#'   panvar.table.list = tables,
#'   pvals.in.log = FALSE,
#'   window = 25,
#'   sig.line = 6,
#'   orient = "V")
#' 
#' # clean up
#' unlink(temp.dir, recursive = TRUE)
plot_panvar_manhattan <- function(panvar.table.list = NULL,
                                  gwas.res = NULL,
                                  ld.list = NULL,
                                  pvals.in.log = TRUE,
                                  plot.r2.thresh = .2,
                                  unplotted.alpha = .4,
                                  window,
                                  sig.line = 3,
                                  orient = c("H", "V"),
                                  point.shape.variable = NULL,
                                  point.shape.scale = NULL,
                                  point.fill.variable.c = NULL,
                                  point.fill.scale.c = NULL,
                                  point.fill.variable.d = NULL,
                                  point.fill.scale.d = NULL,
                                  plot.text.size = 11,
                                  plot.legend.size = 1.2,
                                  snp.highlight.df = NULL,
                                  snp.highlight.point.size = 4,
                                  snp.highlight.shape.var = NULL,
                                  snp.highlight.shape.scale = NULL,
                                  snp.highlight.color.var = NULL,
                                  snp.highlight.color.scale = NULL){
  
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
    marker.list.in.window <- panvar.table.list$window.snps %>%
      # mutate(POS = get_bp_from_id(.data$marker.ID)) %>%
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
    if(!is.null(snp.highlight.df)){
      snp.highlight.df <- snp.highlight.df %>% 
        mutate(PVAL = -log10(.data$PVAL))
    }
  }
  
  # how far to spread labels past ends, in percentage
  y.spread.expansion <- .1
  y.spread.factors <- c(1 + y.spread.expansion, 1 - y.spread.expansion)
  y.spread.factor.window <- (window * y.spread.expansion) * 1000
  
  # plot limits
  this.pos <- get_bp_from_id(key.snp)
  plot.limits <- c(this.pos + window * 1000, this.pos - window * 1000)
  plot.limits.ex <- c(plot.limits[1] + y.spread.factor.window, plot.limits[2] - y.spread.factor.window)
  
  # filter to only in window
  plot.df <- plot.df %>% 
    filter(between(.data$POS, this.pos - window * 1000, this.pos + window * 1000)) 
  
  # ------------------------------------------------------------------------\
  # prepare qualitative variable --------
  # ------------------------------------------------------------------------\
  
  
  # If we're using IMPACT then use this scale
  if(is.null(point.shape.scale) & is.null(point.shape.variable)){
    # do nothing
  } else if(is.null(point.shape.scale) & point.shape.variable == "IMPACT"){
    qual.vars <-  c("HIGH", "MODERATE", "LOW", "MODIFIER")
    point.shape.scale <- make_consistent_scale(values = c(24, 22, 25, 23),
                                                     vars = qual.vars,
                                                     type = "shape",
                                                     name = point.shape.variable)
    
    plot.df$IMPACT <- factor(plot.df$IMPACT, levels = qual.vars)
  } else if(is.null(point.shape.scale) & point.shape.variable == "IMPACT_PLUS"){
    qual.vars <-  c("HIGH", "MODERATE", "LOW", "MODIFIER_CODING", "MODIFIER_INTERGENIC")
    point.shape.scale <- make_consistent_scale(values = c(24, 22, 25, 23, 21),
                                                     vars = qual.vars,
                                                     type = "shape",
                                                     name = point.shape.variable)
    
    plot.df$IMPACT_PLUS <- factor(plot.df$IMPACT_PLUS, levels = qual.vars)
    
  } else if(is.null(point.shape.scale)){
    # make one using the unique values of qualitative annotation column
    qual.vars <- sort(unique(as.matrix(plot.df[,point.shape.variable])))
    
    # throw an error if there are more than 5 things for now
    if(length(qual.vars) > 5){
      stop("Qualitative variable must not have more than 5 unique values. 
           Points are plotted using the fill aesthetic to give them a color, 
           R only has 5 shapes (21-25) that can be assigned fill values.")
    }
    point.shape.scale <- make_consistent_scale(values = c(21:25)[1:length(qual.vars)],
                                                     vars = qual.vars,
                                                     type = "shape",
                                                     name = point.shape.variable)
  }
  
  # ------------------------------------------------------------------------\
  # prepare quantitative variable --------
  # ------------------------------------------------------------------------\
  
  # check for only providing one quantitative scale
  if(!is.null(point.fill.variable.c) & !is.null(point.fill.variable.d)){
    stop("Only one of point.fill.variable.c or point.fill.variable.d should be provided.")
  }
  # 
  if(!is.null(point.fill.scale.c) & !is.null(point.fill.scale.d)){
    stop("Only one of point.fill.variable.c or point.fill.variable.d should be provided.")
  }
  
  # check if a scale is provided but no variable name
  if(is.null(point.fill.variable.c) & !is.null(point.fill.scale.c)){
    warning("No quantitative annotation specified. Provided quantitative.fill.scale ignored.")
  }
  if(is.null(point.fill.variable.d) & !is.null(point.fill.scale.d)){
    warning("No quantitative annotation specified. Provided quantitative.fill.scale ignored.")
  }
  
  
  
  
  # If quantitative annotation is continuous: 
  if(!is.null(point.fill.variable.c)){
    # check if the scale provided is a ggplot scale or not
    if("Scale" %in% class(point.fill.scale.c)){
      user.supplied.scale.option <- "scale"
    } else {
      user.supplied.scale.option <- "not.scale"
    }
    
    if(!is.null(point.fill.variable.c) & user.supplied.scale.option == "scale"){
      # Both variable and scale provided 
      # use user supplied scale
      quantitative.fill.scale <- point.fill.scale.c
      # only plot that stuff that's above the LD threshold 
      plot.df <- plot.df %>% 
        mutate(plot.quant.var = case_when(.data$R2 < plot.r2.thresh ~ NA,
                                          TRUE ~ .data[[point.fill.variable.c]]))
    } else if(!is.null(point.fill.variable.c)){
      # Only variable provided 
      # set a binned color scale
      if(is.null(point.fill.scale.c)){
        quantitative.fill.scale <- scale_fill_viridis_b(name = point.fill.variable.c)
      } else {
        quantitative.fill.scale <- scale_fill_viridis_b(name = point.fill.variable.c, option = point.fill.scale.c)
      }
      # only plot that stuff that's above the LD threshold 
      plot.df <- plot.df %>% 
        mutate(plot.quant.var = case_when(.data$R2 < plot.r2.thresh ~ NA,
                                          TRUE ~ .data[[point.fill.variable.c]]))
    } else {
      # do nothing
    }
    quantitative.annotation <- point.fill.variable.c
  } else if(!is.null(point.fill.variable.d)){
    # check if the scale provided is a ggplot scale or not
    if("Scale" %in% class(point.fill.scale.c)){
      user.supplied.scale.option <- "scale"
    } else {
      user.supplied.scale.option <- "not.scale"
    }
    
    if(!is.null(point.fill.variable.d) & user.supplied.scale.option == "scale"){
      # Both variable and scale provided 
      # use user supplied scale
      quantitative.fill.scale <- point.fill.scale.d
      # only plot that stuff that's above the LD threshold 
      plot.df <- plot.df %>% 
        mutate(plot.quant.var = case_when(.data$R2 < plot.r2.thresh ~ NA,
                                          TRUE ~ .data[[point.fill.variable.d]]))
    } else if(!is.null(point.fill.variable.d)){
      # Only variable provided 
      # set a binned color scale
      if(is.null(point.fill.scale.d)){
        quantitative.fill.scale <- scale_fill_viridis_d(name = point.fill.variable.d, 
                                                        na.translate = F)
      } else {
        quantitative.fill.scale <- scale_fill_viridis_d(name = point.fill.variable.d, 
                                                        option = point.fill.scale.d,
                                                        na.translate = F)
      }
      # only plot that stuff that's above the LD threshold 
      plot.df <- plot.df %>% 
        mutate(plot.quant.var = case_when(.data$R2 < plot.r2.thresh ~ NA,
                                          TRUE ~ .data[[point.fill.variable.d]]))
    } else {
      # do nothing
    }
    quantitative.annotation <- point.fill.variable.d
  
  }
  
  # make sure quantitative.annotation has a value if not provided
  if(is.null(point.fill.variable.c) & is.null(point.fill.variable.d)){
    quantitative.annotation <- NULL
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
      panel.grid = element_blank(),
      text = ggplot2::element_text(size = plot.text.size),
      legend.key.size = unit(plot.legend.size, "lines")
    ) +
    geom_hline(yintercept = sig.line, linetype = 'dashed')
  
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
    guide = guide_colorsteps(order = 1),
    na.value = "grey50"
  )
  
  # Some logic for including extra variables
  if(!is.null(point.shape.variable) & !is.null(quantitative.annotation)){
    # Use a quant and a qual
    man <- man +
      geom_point(aes(fill = .data$plot.quant.var, alpha = .data$how.to.plot, shape = .data[[point.shape.variable]]),
                 size = 3, color = "black") +
      scale_alpha_continuous(guide = "none", limits = c(0,NA), range = c(0, 1)) +
      quantitative.fill.scale +
      point.shape.scale
    
  } else if(is.null(point.shape.variable) & !is.null(quantitative.annotation)){
    # Use just a quant 
    man <- man +
      geom_point(aes(fill = .data$plot.quant.var, alpha = .data$how.to.plot), shape = 21, size = 3, color = "black") +
      scale_alpha_continuous(guide = "none", limits = c(0,NA), range = c(0, 1)) +
      quantitative.fill.scale
    
  } else if(!is.null(point.shape.variable) & is.null(quantitative.annotation)){
    # Use just a qual
    man <- man + 
      geom_point(aes(fill = .data$plot.R2, alpha = .data$how.to.plot, shape = .data[[point.shape.variable]]), size = 3, color = "black") +
      point.shape.scale +
      scale_alpha_continuous(guide = "none", limits = c(0,NA), range = c(0, 1)) +
      default.LD.fill.scale
    
  } else {
    # Use neither
    man <- man +
      geom_point(aes(fill = .data$plot.R2, alpha = .data$how.to.plot), size = 3, shape = 21, color = "black") +
      scale_alpha_continuous(guide = "none", limits = c(0,NA), range = c(0, 1)) +
      default.LD.fill.scale
  }
  
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

  if(!is.null(snp.highlight.df)){
    if(!is.null(snp.highlight.shape.var) & !is.null(snp.highlight.color.var)){
      # both shape and color provided 
      man <- man +
        geom_point(data = snp.highlight.df, size = snp.highlight.point.size, aes(shape = .data[[snp.highlight.shape.var]], color = .data[[snp.highlight.color.var]]))
    } else if(!is.null(snp.highlight.shape.var) & is.null(snp.highlight.color.var)){
      # only shape provided 
      man <- man +
        geom_point(data = snp.highlight.df, size = snp.highlight.point.size, aes(shape = .data[[snp.highlight.shape.var]], color = 'red'))
    } else if(is.null(snp.highlight.shape.var) & !is.null(snp.highlight.color.var)){
      # only color provided 
      man <- man +
        geom_point(data = snp.highlight.df, size = snp.highlight.point.size, aes(color = .data[[snp.highlight.color.var]]))
    } else {
      # neither provided 
      man <- man +
        geom_point(data = snp.highlight.df, size = snp.highlight.point.size, color = "red")
    }
  }
  
  return(man)
}
