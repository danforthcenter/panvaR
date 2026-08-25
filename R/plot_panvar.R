# ------------------------------------------------------------------------\
# main function --------
# ------------------------------------------------------------------------\

#' Make a plot that lines up manhattan and gene locations
#'
#' @param panvar.table.list list, output from [panvaR::make_panvar_tables]. Provide either this list or both gwas.res and ld.list.
#' @param gwas.res data.frame of all gwas results, should contain columns (CHR,
#'   POS, PVAL), corresponding to (chromosome, physical position, and pvalue).
#' @param annotation.table table with annotations with columns (geneID, CHR,
#'   start, end, annotation). start and end correspond to base-pair coordinates
#'   of start and end of gene. CHR is chromosome of gene.
#' @param ld.list list, output of [panvaR::get_ld_in_window]
#' @param pvals.in.log boolean, if TRUE PVAL column has already been converted to -log10(pvalue)
#' @param plot.r2.thresh minimum LD with qtl snps to plot snps colored by LD
#' @param unplotted.alpha numeric, number from 0 to 1 to indicate alpha values
#'   of snps below the plot.r2.thresh. To not plot these snps set value to 0.
#' @param window numeric, total window size in KB, all variants within .5 *
#'   window are calculated.
#' @param sig.line numeric, -log10(p) value to draw line on plot
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
#' @param include.gene.id boolean, if TRUE, `geneID` column will be included in
#'   annotation plot.
#' @param highlight.gene.ids character, vector of geneID's that will be
#'   highlighted in the plot.
#' @param gene.highlight.color character, a color to highlight specific geneIDs
#' @param annotation.point.variable character, variable in  `annotation.table` that indicates how to color points plotted next to gene
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
#' plot_panvar(
#'   panvar.table.list = tables,
#'   pvals.in.log = FALSE,
#'   window = 25,
#'   sig.line = 6)
#'   
#' # clean up 
#' unlink(temp.dir, recursive = TRUE)
plot_panvar <- function(panvar.table.list = NULL,
                        gwas.res = NULL,
                        ld.list = NULL,
                        annotation.table = NULL,
                        pvals.in.log = T,
                        plot.r2.thresh = .2,
                        unplotted.alpha = .4,
                        window,
                        sig.line,
                        point.shape.variable = NULL,
                        point.shape.scale = NULL,
                        point.fill.variable.c = NULL,
                        point.fill.scale.c = NULL,
                        point.fill.variable.d = NULL,
                        point.fill.scale.d = NULL,
                        plot.title = "",
                        include.gene.id = F,
                        highlight.gene.ids = NULL,
                        gene.highlight.color = "red",
                        annotation.point.variable = "LD",
                        annotation.point.scale = NULL,
                        plot.effect = F,
                        plot.text.size = 11,
                        plot.legend.size = 1.2) {
  
  # ------------------------------------------------------------------------\
  # make manhattan --------
  # ------------------------------------------------------------------------\
  
  message("Making manhattan")
  man <- plot_panvar_manhattan(panvar.table.list = panvar.table.list, 
                               gwas.res = gwas.res,
                               pvals.in.log = pvals.in.log,
                               plot.r2.thresh = plot.r2.thresh,
                               unplotted.alpha = unplotted.alpha,
                               ld.list = ld.list,
                               window = window,
                               sig.line = sig.line,
                               orient = "H",
                               point.shape.variable = point.shape.variable,
                               point.shape.scale = point.shape.scale,
                               point.fill.variable.c = point.fill.variable.c,
                               point.fill.scale.c = point.fill.scale.c,
                               point.fill.variable.d = point.fill.variable.d,
                               point.fill.scale.d = point.fill.scale.d,
                               plot.text.size = plot.text.size,
                               plot.legend.size = plot.legend.size)
  
  # ------------------------------------------------------------------------\
  # make annotation --------
  # ------------------------------------------------------------------------\
  
  message("Making annotation plot")
  # setting column to pull fill variable from
  # setting the scale
  if(is.null(annotation.point.variable)){
    # do nothing
  } else if(annotation.point.variable == "LD"){
    annotation.point.scale <- default.panvar.LD.scale(type = "color", plot.r2.thresh = plot.r2.thresh)
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
  
  # if(is.null(ld.list$key.snp)){
  #   middle.snp <- NULL
  # } else {
  #   middle.snp <- paste0(get_chrom_from_id(ld.list$key.snp),
  #                        "_",
  #                        get_bp_from_id(ld.list$key.snp))
  # }

  
  
  anno <- plot_gene_annotation(
    panvar.table.list = panvar.table.list,
    annotation.table = annotation.table,
    middle.snp = ld.list$key.snp,
    window = window,
    include.id = include.gene.id,
    highlight.ids = highlight.gene.ids,
    highlight.color = gene.highlight.color,
    use.arrows = F,
    point.color = annotation.point.variable,
    point.fill.scale = annotation.point.scale,
    plot.text.size = plot.text.size,
    plot.legend.size = plot.legend.size
  )
  
  # ------------------------------------------------------------------------\
  # make effect --------
  # ------------------------------------------------------------------------\
  
  if (plot.effect) {
    effect.plot <- plot_effect(
      panvar.table.list = panvar.table.list,
      gwas.res = gwas.res,
      pvals.in.log = pvals.in.log,
      plot.r2.thresh = plot.r2.thresh,
      unplotted.alpha = unplotted.alpha,
      ld.list = ld.list,
      window = window,
      sig.line = sig.line,
      orient = "V",
      point.shape.variable = point.shape.variable,
      point.shape.scale = point.shape.scale,
      point.fill.variable.c = point.fill.variable.c,
      point.fill.scale.c = point.fill.scale.c,
      point.fill.variable.d = point.fill.variable.d,
      point.fill.scale.d = point.fill.scale.d,
      include.legend = F
    )
    
  }
  
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

