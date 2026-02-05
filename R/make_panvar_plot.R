

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

#' Make a plot that lines up manhattan and gene locations
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
#' @param pvals.in.log
#' @param geno.bed.filename character, prefix of genotype files in plink
#'   (bed/bim/fam) format. Do not include ".bed" extension.
#' @param geno.bed.directory character, directory where genotype files are
#'   located
#' @param temp.dir character, where to output some temporary files.
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
                             annotation.point.variable.location = c("gwas", "anno"),
                             annotation.point.variable = "LD",
                             annotation.point.scale = NULL,
                             compute.scores = F,
                             score.vars = NULL,
                             score.dirs = NULL,
                             score.weights = NULL,
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
  # subset gwas  --------
  # ------------------------------------------------------------------------\
  
  # filter gwas df to just in window and join LD
  gwas.sub <- gwas.res %>%
    as.data.frame() %>% 
    mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
    left_join(ld.list$table, by = "marker.ID") %>%
    filter(!is.na(.data$R2)) %>% 
    # rename this here to work with special case "LD" given to point variable
    rename("LD" = "R2")
  
  # add middlesnp back in
  gwas.sub.mid.snp <- gwas.res %>% 
    filter(marker.ID == middle.snp) %>% 
    mutate(LD = 1)
  gwas.sub <- bind_rows(gwas.sub, gwas.sub.mid.snp)
  
  # ------------------------------------------------------------------------\
  # make scores --------
  # ------------------------------------------------------------------------\
  
  # make_scores <- function(input.df, cols, directions, weights = NULL){
  if(compute.scores){
    message("Computing scores")
    # add dist column 
    score.in <- gwas.sub %>% 
      mutate(key.snp.pos = get_bp_from_id(middle.snp)) %>% 
      mutate(DIST = abs(POS - key.snp.pos)) 
    
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
  
  annotation.point.variable.location <- match.arg(annotation.point.variable.location)
  
  if(annotation.point.variable.location == "gwas"){
    if(!is.null(annotation.point.variable)){
      message("Generating snp to gene correspondence")
      
      # filter anno to just window
      this.chrom <- get_chrom_from_id(ld.list$key.snp)
      this.pos <- get_bp_from_id(ld.list$key.snp)
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
  } else {
    anno.in <- annotation.table
  }
  
  # ------------------------------------------------------------------------\
  # make annotation --------
  # ------------------------------------------------------------------------\
  
  # setting column to pull fill variable from
  # setting the scale
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
    if(annotation.point.variable.location == "anno"){
      point.color.option <- annotation.point.variable
    } else {
      anno.in <- anno.in %>% 
        rename_with(~annotation.point.variable, "maximum.value")
      point.color.option <- annotation.point.variable
      # point.color.option <- "maximum.value"
    }
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

