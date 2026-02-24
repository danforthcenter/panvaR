# ------------------------------------------------------------------------\
# test with tables list input --------
# ------------------------------------------------------------------------\

snpeffann <- read.csv("~/scratch/panvar_test/setaria_annotated_chrsnumeric_snpeffinfo_1persnp.csv") %>% 
  select(-CHROM, -POS)

# set global options 
set_plink_path("/home/cluebbert/bin/plink2")
set_panvar_prefix("SetShattering")
set_out_dir("~/scratch/panvar_test/setaria_shatter_full")

gwas.df <- data.table::fread(file.path(options()$panvar_outdir, "SetShattering_GLM_GWASresults.csv")) %>% 
  left_join(snpeffann, by = "marker.ID") %>% 
  filter(str_detect(CHR, "scaffold", negate = T))
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

qtl.df.test <- gwas.df %>% 
  dplyr::slice_max(LOGPVAL, n = 3)

tables <- make_panvar_tables(gwas.res = gwas.df,
                             tag.snp = "5-6857045",
                             # qtl.df = qtl.df.test, 
                             annotation.table = anno,
                             plink.path = options()$plink_path,
                             pvals.in.log = F,
                             geno.bed.filename = "SetShattering_PlinkQC_maf0.05_missing0.1",
                             geno.bed.directory = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full/",
                             window = 500,
                             compute.scores = T,
                             snp.to.gene.vars = c("LD", "snp.score"),
                             snp.to.gene.buffer = 5)

# ------------------------------------------------------------------------\
# test plotly compatible --------
# ------------------------------------------------------------------------\

# make manhattan
man <- 
plot_panvar_manhattan(panvar.table.list = tables,
                      pvals.in.log = F,
                      window = 20, 
                      sig.line = 6)

man

anno <- 
  plotly_gene_annotation(panvar.table.list = tables,
                            window = 20, 
                       point.color = "snp.score")
anno
ggplotly(anno)
ggp <- ggplotly(anno,tooltip = c("text", "snp.score"))
ggp$x$data[[4]]$textposition <- "right"
ggp


plotly::subplot(man, ggp, nrows = 1, shareY = T)

# ------------------------------------------------------------------------\
# example --------
# ------------------------------------------------------------------------\

# from plotly.subplots import make_subplots
# import plotly.graph_objects as go
# 
# fig = make_subplots(rows=1, cols=2)
# 
# fig.add_trace(
#   go.Scatter(x=[1, 2, 3], y=[4, 5, 6]),
#   row=1, col=1
# )
# 
# fig.add_trace(
#   go.Scatter(x=[20, 30, 40], y=[50, 60, 70]),
#   row=1, col=2
# )
# 
# fig.update_layout(height=600, width=800, title_text="Side By Side Subplots")
# fig.show()