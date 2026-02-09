# ------------------------------------------------------------------------\
# main function --------
# ------------------------------------------------------------------------\

#' Make a plot that lines up manhattan and gene locations
#'
#' @param panvar.table.list list, output from [panvaR::make_panvar_tables]. Provide either this list or both gwas.res and ld.list.
#' @param gwas.res data.frame of all gwas results, should contain columns (CHR,
#'   POS, PVAL), corresponding to (chromosome, physical position, and pvalue).
#' @param tag.snp character, marker.ID of snp around which to calculate LD. In
#'   the form 'CHR-POS'
#' @param annotation.table table with annotations with columns (geneID, CHR,
#'   start, end, annotation). start and end correspond to base-pair coordinates
#'   of start and end of gene. CHR is chromosome of gene.
#' @param pvals.in.log boolean, if TRUE PVAL column has already been converted to -log10(pvalue)
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD
#' @param unplotted.alpha numeric, number from 0 to 1 to indicate alpha values
#'   of snps below the plot.r2.thresh. To not plot these snps set value to 0.
#' @param window numeric, total window size in KB, all variants within .5 *
#'   window are calculated.
#' @param sig.line numeric, -log10(p) value to draw line on plot
#' @param orient character, will rotate plot 90 degrees. vertical (V) or
#'   horizontal (H) refers to how the "buildings" of the plot are plotted. "V"
#'   places pvalue on y-axis, "H" places pvalues on x-axis.
#' @param qualitative.annotation character, column in `gwas.res` that contains
#'   qualitative annotations. For example impact grades from snpeff. See
#'   [panvaR::format_snpeff_annotations]. Will be plotted as shapes. Only
#'   accepts up to 5 classes. "IMPACT" and "IMPACT_PLUS" are special cases that
#'   will have a pre-assigned scale used if supplied here.
#' @param qualitative.shape.scale ggplot scale, an object with a stored call to
#'   [ggplot2::scale_shape_manual]. More often an output of the function
#'   [panvaR::make_consistent_scale].
#' @param quantitative.annotation character, column in `gwas.res` that contains
#'   quantitative annotations. For example, variant effect scores. Will be
#'   plotted as fill to points.
#' @param quantitative.fill.scale character or scale object, either a character
#'   indicating the `option` parameter passed to [ggplot2::scale_fill_viridis_b]
#'   that alters the color scale used. Or a previous call to a ggplot2 fill
#'   scale for example [ggplot2::scale_fill_stepsn].
#' @param plot.title character, a title for the plot
#' @param include.gene.id boolean, if TRUE, `geneID` column will be included in
#'   annotation plot.
#' @param highlight.gene.ids character, vector of geneID's that will be
#'   highlighted in the plot.
#' @param gene.highlight.color character, a color to highlight specific geneIDs
#' @param annotation.point.variable.location character, either "anno" or "gwas"
#'   to indicate the location of the column indicated in
#'   `annotation.point.variable` in either `annotation.table` or `gwas.res`
#'   respectively. If column is located in `gwas.res`, the maximum value of all
#'   snps in a given gene will be used.
#' @param annotation.point.variable character, variable in either `gwas.res` or
#'   `annotation.table` that indicates how to color points plotted next to gene
#'   descriptions. If not supplied, no points are plotted. The input "LD" is
#'   reserved and will use LD.
#' @param annotation.point.scale ggplot2 scale object, a color scale to customize
#'   how point.color is displayed.
#' @param plot.effect boolean, if TRUE include volcano style effect vs pvalue
#'   plot as inset.
#'
#' @returns ggplot2 object of plot with manhattan plot alongside genes for a
#' given genomic window.
#' @export
#'
#' @examples
#' # work in progress
make_panvar_plot <- function(panvar.table.list = NULL, 
                             gwas.res = NULL,
                             ld.list = NULL,
                             annotation.table = NULL,
                             pvals.in.log = T,
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
                             annotation.point.variable.location = c("gwas", "anno"),
                             annotation.point.variable = "LD",
                             annotation.point.scale = NULL,
                             plot.effect = F) {
  
  # ------------------------------------------------------------------------\
  # make manhattan --------
  # ------------------------------------------------------------------------\
  
  message("Making manhattan")
  man <- make_panvar_manhattan(panvar.table.list = panvar.table.list, 
                               gwas.res = gwas.res,
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
  # make annotation --------
  # ------------------------------------------------------------------------\
  
  message("Making annotation plot")
  # setting column to pull fill variable from
  # setting the scale
  if(is.null(annotation.point.variable)){
    # do nothing
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
  }
  
  # make sure the variable is in the table if its provided 
  if(!is.null(annotation.point.variable)){
    if(is.null(panvar.table.list)){
      # Check the point variable is in annotation table if not providing panvar.table.list
      if(!annotation.point.variable %in% names(annotation.table)){
        stop(paste0("annotation.point.variable: ", annotation.point.variable, 
                    " not found in annotation.table."))
      }
    } else {
      if(!annotation.point.variable %in% names(panvar.table.list$anno)){
        stop(paste0("annotation.point.variable: ", annotation.point.variable, 
                    " not found in panvar.table.list$anno"))
      }
    }
  }
  
  # if(is.null(annotation.point.variable)){
  #   point.color.option <- NULL
  # } else if(annotation.point.variable == "LD"){
  #   default.LD.fill.scale <- binned_scale(
  #     aesthetics = "fill",
  #     name = "R2 \n",
  #     palette = function(x)
  #       c("#43638E", "#88DAA0", "#DBC32D", "#B94712"),
  #     limits = c(plot.r2.thresh, 1),
  #     breaks = seq(plot.r2.thresh, 1, length.out = 5)[-c(1, 5)],
  #     show.limits = T,
  #     guide = "colorsteps",
  #     na.value = "grey50"
  #   )
  #   annotation.point.scale <- default.LD.fill.scale
  #   # anno.in <- anno.in %>% 
  #   #   rename("LD" = "maximum.value")
  #   point.color.option <- "LD"
  # } else if(!is.null(annotation.point.variable)) {
  #   point.color.option <- annotation.point.variable
  # } 
  
  
  anno <- make_gene_annotation_plot(panvar.table.list = panvar.table.list,
                                    annotation.table = annotation.table,
                                    middle.snp = ld.list$key.snp,
                                    window = window,
                                    include.id = include.gene.id,
                                    highlight.ids = highlight.gene.ids,
                                    highlight.color = gene.highlight.color,
                                    use.arrows = F,
                                    point.color = annotation.point.variable,
                                    point.fill.scale = annotation.point.scale)
  
  # ------------------------------------------------------------------------\
  # make effect --------
  # ------------------------------------------------------------------------\
  
  # if (plot.effect) {
  #   effect.plot <- make_effect_plot(gwas.res,
  #                                   qtl.df,
  #                                   pvals.in.log,
  #                                   plot.r2.thresh,
  #                                   ld.list,
  #                                   window,
  #                                   sig.line)
  #   
  # }
  
  # ------------------------------------------------------------------------\
  # final plot --------
  # ------------------------------------------------------------------------\
  
  # chrom <- unique(qtl.df$CHR)
  
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

