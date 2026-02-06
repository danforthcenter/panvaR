
# set global options 
set_plink_path("/home/cluebbert/bin/plink2")
set_panvar_prefix("SetShattering")
set_out_dir("~/scratch/panvar_test/setaria_shatter_full")


# ------------------------------------------------------------------------\
# test sub modules --------
# ------------------------------------------------------------------------\

# ld
out_ld <- 
  panvaR::get_ld_in_window(tag.snp = "5-6857045",
                           window = 500,
                           in.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full/",
                           out.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full",
                           geno.bed = "SetShattering_PlinkQC_maf0.05_missing0.1")

gwas.df <- data.table::fread(file.path(options()$panvar_outdir, "SetShattering_GLM_GWASresults.csv"))

# manhattan
make_panvar_manhattan2(gwas.res = gwas.df,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      window = 10,
                      sig.line = 6, orient = "V")

# annotation
anno <- read.csv("~/scratch/setaria_biomart.txt") %>% 
  filter(str_detect(Chromosome.Name, "scaffold", negate = T)) %>% 
  mutate(CHR = as.numeric(str_replace(Chromosome.Name, "Chr_", ""))) %>% 
  select(CHR, geneID = Gene.Name, 
         start = Gene.Start..bp.,
         end = Gene.End..bp.,
         annotation = Description) %>% 
  mutate(annotation = case_when(annotation == "" ~ "No gene description.",
                                TRUE ~ annotation)) %>% 
  distinct()

make_gene_annotation_plot(annotation.table = anno,
                          middle.snp = out_ld$key.snp,
                          window = 250, 
                          include.id = T,
                          use.arrows = F)

# ------------------------------------------------------------------------\
# test full plot --------
# ------------------------------------------------------------------------\

# p <- 
make_panvar_plot(
  gwas.res = gwas.df,
  tag.snp = "5-6857045",
  annotation.table = anno,
  plink.path = options()$plink_path,
  pvals.in.log = F,
  geno.bed.filename = "SetShattering_PlinkQC_maf0.05_missing0.1",
  geno.bed.directory = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full/",
  plot.r2.thresh = .2,
  unplotted.alpha = .4,
  window = 100,
  sig.line = 6,
  orient = "H",
  qualitative.annotation = NULL,
  qualitative.shape.scale = NULL,
  quantitative.annotation = NULL,
  quantitative.fill.scale = NULL,
  plot.title = "",
  include.gene.id = T,
  highlight.gene.ids = NULL,
  gene.highlight.color = "red",
  annotation.point.scale = ggplot2::scale_fill_viridis_b(option = "plasma"),
  plot.effect = F
)

# doesn't work :(
ggplotly(p)

# ------------------------------------------------------------------------\
# test full plot with some other values --------
# ------------------------------------------------------------------------\

gwas.df <- data.table::fread("~/scratch/panvar_test/setaria_shatter_full/Shatter_GLM_res_with_snpeff_and_cad.csv",
                             data.table = F)

make_panvar_plot(
  gwas.res = gwas.df,
  tag.snp = "5-6857045",
  annotation.table = anno,
  plink.path = options()$plink_path,
  pvals.in.log = F,
  geno.bed.filename = "SetShattering_PlinkQC_maf0.05_missing0.1",
  geno.bed.directory = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full/",
  plot.r2.thresh = .2,
  unplotted.alpha = .4,
  window = 100,
  sig.line = 6,
  orient = "H",
  qualitative.annotation = "IMPACT_PLUS",
  qualitative.shape.scale = NULL,
  # quantitative.annotation = "zero_shot_positive",
  quantitative.fill.scale = NULL,
  plot.title = "Setaria Shattering",
  include.gene.id = T,
  highlight.gene.ids = NULL,
  gene.highlight.color = "red",
  annotation.point.variable = "snp.score",
  annotation.point.scale = ggplot2::scale_color_viridis_b(option = "plasma"),
  plot.effect = F,
  compute.scores = T
)

luebbert::memesave("~/OneDrive/vm_transfer/SetariaShattering_panvar_snp.score.png")

# ------------------------------------------------------------------------\
# test make input to points in anno --------
# ------------------------------------------------------------------------\

gwas.res <- gwas.df 
ld.list <- out_ld
middle.snp <- out_ld$key.snp
annotation.table <- anno
window <- 25

# filter gwas df to just in window
gwas.sub <- gwas.res %>%
  as.data.frame() %>% 
  mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
  left_join(ld.list$table, by = "marker.ID") %>%
  filter(!is.na(.data$R2))

# filter anno to just window
this.chrom <- get_chrom_from_id(middle.snp)
this.pos <- get_bp_from_id(middle.snp)
anno.sub <- annotation.table %>%
  filter(.data$CHR == this.chrom) %>%
  rowwise() %>%
  mutate(dist.from.snp = get.gene.dist.from.snp(this.pos, .data$start, .data$end)) %>%
  filter(.data$dist.from.snp <= window * 1000) %>%
  mutate(gene.mid = median(c(.data$start, .data$end))) %>% 
  arrange(.data$gene.mid)

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

# marker in 2 genes "5-6867941"
library(dplyr)

test <- gwas.sub %>% 
  rowwise() %>% 
  mutate(snp.in.gene = get.gene.from.snp(POS, anno.sub)) %>% 
  filter(!is.null(snp.in.gene)) %>% 
  unnest_longer(snp.in.gene) %>% 
  group_by(snp.in.gene) %>% 
  summarize(maximum.value = max(R2))

# make it work with any number of variables, always do LD...
variables <- c("LOGPVAL", "R2")
test <- gwas.sub %>% 
  rowwise() %>% 
  mutate(snp.in.gene = get.gene.from.snp(POS, anno.sub)) %>% 
  filter(!is.null(snp.in.gene)) %>% 
  unnest_longer(snp.in.gene) %>% 
  group_by(snp.in.gene)

x <- test %>% 
  summarize(across(all_of(variables), ~ max(.x, na.rm = T)))  


# overlapping genes
bp <- gwas.sub$POS[9188]

between(bp, anno.sub$start, anno.sub$end)

