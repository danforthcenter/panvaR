# make annotated gene list plot

# ------------------------------------------------------------------------\
# unexported functions --------
# ------------------------------------------------------------------------\

# get the minmum distance from a gene to a snp
get.gene.dist.from.snp <- function(s, start, end){
  a <- abs(s - start)
  b <- abs(s - end)
  return(min(a, b))
}

# make evenly spaced values across a range for a given number of items
get.gene.y.pos <- function(from, to, length.out) {
  length.out <- length.out + 2
  result <- seq(from, to, length.out = length.out)
  result <- result[-c(1,length(result))]
  return(result)
}

# ------------------------------------------------------------------------\
# main function --------
# ------------------------------------------------------------------------\


#' Plot genes and their locations 
#' 
#' Use a table of genes, descriptions and their start and end points to plot their locations 
#' and annotations along the genome. 
#'
#' @param annotation.table table with annotations with columns (geneID, CHR, start, end, annotation). start and end correspond to base-pair coordinates of start and end of gene. CHR is chromosome of gene.
#' @param middle.snp character, SNP name in form "CHR-POS" center of window. Often the key.snp output of [panvaR::get_ld_in_window]
#' @param window integer, kilobases on either side of middle.snp to plot
#' @param include.id boolean, include geneID in gene annotations or not
#' @param highlight.ids character, optional, vector of ids to highlight
#' @param highlight.color character, optional, color to highlight ids
#' @param use.arrows boolean, if TRUE, use [gggenes::geom_gene_arrow] to draw representations of genes. 
#' If the direction is encoded in another variable, supply that to <NEW ARGUMENT NAME>
#' @param point.color character, variable in annotation.table that indicates how to color 
#' points plotted next to gene descriptions. If not supplied, no points are plotted. 
#' The input "LD" is reserved to give functionality to [panvaR::make_panvar_plot].
#' If used, legend will not be displayed. 
#' @param point.fill.scale ggplot2 scale object, a fill scale to customize how point.color 
#' is displayed. 
#'
#' @returns
#' GGplot object 
#' @export
#'
#' @examples
#' # Work in progress
make_gene_annotation_plot <- function(annotation.table,
                                      middle.snp,
                                      window,
                                      include.id = F,
                                      gene.color = "blue",
                                      highlight.ids = NULL,
                                      highlight.color = "red",
                                      use.arrows = F,
                                      point.color = NULL,
                                      point.fill.scale = NULL){
  
  
  # get pos and chrom
  this.chrom <- stringr::str_extract(middle.snp, "^(.*?)-", group = 1)
  this.pos <- get_bp_from_id(middle.snp)
  
  # how far to spread labels past ends, in percentage
  y.spread.expansion <- .1
  y.spread.factors <- c(1 + y.spread.expansion, 1 - y.spread.expansion)
  y.spread.factor.window <- (window * y.spread.expansion) * 1000
  
  # plot limits
  plot.limits <- c(this.pos + window * 1000, this.pos - window * 1000)
  plot.limits.ex <- c(plot.limits[1] + y.spread.factor.window, plot.limits[2] - y.spread.factor.window)
  
  
  # make gene plot
  anno.sub <- annotation.table %>%
    # select("geneID", "CHR", "start", "end", "annotation") %>%
    filter(.data$CHR == this.chrom) %>%
    rowwise() %>%
    mutate(dist.from.snp = get.gene.dist.from.snp(this.pos, .data$start, .data$end)) %>%
    filter(.data$dist.from.snp <= window * 1000) %>%
    mutate(id.plus.anno = paste0(.data$geneID, ", ", .data$annotation),
           plot.label = ifelse(include.id, .data$id.plus.anno, .data$annotation)) %>% 
    mutate(gene.mid = median(c(.data$start, .data$end))) %>% 
    arrange(.data$gene.mid)
  
  anno.spread <- anno.sub %>%
    tibble::add_column(y.pos = get.gene.y.pos(min(plot.limits.ex),
                                      max(plot.limits.ex),
                                      nrow(anno.sub)))
  
  if(!is.null(highlight.ids)){
    anno.spread <- anno.spread %>%
      mutate(gene.label.color = case_when(geneID %in% highlight.ids ~ "A",
                                          TRUE ~ "B"))
  }
  
  mid <- this.pos
  breaks.anno <- seq(from = mid - window * 1000,
                     to = mid + window * 1000,
                     length.out = 9)
  if(use.arrows){
    anno <-
      ggplot(data = anno.spread) +
      # geom_segment(aes(y = .data$start, yend = .data$end, x = .5, xend = .5), linewidth = 2, color = "red") +
      gggenes::geom_gene_arrow(aes(xmin = .data$start, xmax = .data$end, y = .5), fill = gene.color) +
      scale_x_reverse(limits = plot.limits.ex,
                      labels = function(x) paste0((x - this.pos) / 1000, " KB"),
                      breaks = breaks.anno) +
      ylim(.5, .875) +
      geom_segment(aes(y = .5, yend = .55, x = .data$gene.mid, xend = .data$y.pos)) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            #axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid = element_blank(),
            axis.ticks.x = element_blank()) +
      coord_flip()
  } else {
    anno <-
      ggplot(data = anno.spread) +
      geom_segment(aes(y = .data$start, yend = .data$end, x = .5, xend = .5), linewidth = 2, color = gene.color) +
      scale_y_reverse(limits = plot.limits.ex,
                      labels = function(x) paste0((x - this.pos) / 1000, " KB"),
                      breaks = breaks.anno) +
      xlim(.5, .875) +
      geom_segment(aes(x = .5, xend = .55, y = .data$gene.mid, yend = .data$y.pos)) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            #axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  # Add some points if you want
  if(!is.null(point.color)){
    anno <- anno +
      geom_point(aes(x = .555, y = .data$y.pos, fill = .data[[point.color]]), shape = 21, size = 2) +
      theme(legend.position = "right",
            legend.justification = "top")
    if(!is.null(point.fill.scale)){
      anno <- anno + 
        point.fill.scale
    } else {
      anno <- anno + 
        scale_fill_viridis_b()
    }
    # move the text over a little to accomodate point
    text.x.start <- .56
  } else {
    text.x.start <- .55
  }
  
  # make a special case for LD where we hide the legend
  if(!is.null(point.color)){
    if(point.color == "LD"){
      anno <- anno +
        theme(legend.position = "none")
    }
  }
  
  if(!is.null(highlight.ids)){
    if(use.arrows){
      anno <- anno +
        ggfittext::geom_fit_text(aes(ymin = text.x.start, ymax = .85, x = .data$y.pos, label = .data$plot.label, color = .data$gene.label.color),
                                 place = "left",
                                 #grow = TRUE,
                                 hjust = 0,
                                 padding.y = grid::unit(.1, "lines"),
                                 min.size = 4,
                                 show.legend = F) +
        scale_color_manual(values = c(highlight.color, "black"))
    } else {
      anno <- anno +
        ggfittext::geom_fit_text(aes(xmin = text.x.start, xmax = .85, y = .data$y.pos, label = .data$plot.label, color = .data$gene.label.color),
                                 place = "left",
                                 #grow = TRUE,
                                 hjust = 0,
                                 padding.y = grid::unit(.1, "lines"),
                                 min.size = 4,
                                 show.legend = F) +
        scale_color_manual(values = c(highlight.color, "black"))
    }
  } else {
    if(use.arrows){
      anno <- anno +
        ggfittext::geom_fit_text(aes(ymin = text.x.start, ymax = .85, x = .data$y.pos, label = .data$plot.label),
                                 place = "left",
                                 #grow = TRUE,
                                 hjust = 0,
                                 padding.y = grid::unit(.1, "lines"),
                                 min.size = 4) 
    } else {
      anno <- anno +
        ggfittext::geom_fit_text(aes(xmin = text.x.start, xmax = .85, y = .data$y.pos, label = .data$plot.label),
                                 place = "left",
                                 #grow = TRUE,
                                 hjust = 0,
                                 padding.y = grid::unit(.1, "lines"),
                                 min.size = 4) 
    }
  }
  
  return(anno)
  
}

