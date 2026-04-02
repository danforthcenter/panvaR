

# ------------------------------------------------------------------------\
# look at annotation file --------
# ------------------------------------------------------------------------\

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

anno_ld <- anno %>% 
  mutate(ld_sim = runif(nrow(anno)))

make_gene_annotation_plot(annotation.table = anno_ld,
                          middle.snp = out_ld$key.snp,
                          window = 25, 
                          include.id = T,
                          use.arrows = F,
                          point.color = "LD",
                          point.fill.scale = NULL)

# ------------------------------------------------------------------------\
# horizontal --------
# ------------------------------------------------------------------------\

make_gene_annotation_plot_horiz(anno,
                          middle.snp = out_ld$key.snp,
                          window = 25, 
                          include.id = T,
                          use.arrows = F)

# ------------------------------------------------------------------------\
# setaria --------
# ------------------------------------------------------------------------\

anno <- read.csv("~/scratch/panvar_test/setaria_shatter_full/Anno_setaria_shattering.csv")

anno_cleanup <- anno %>% 
  mutate(annotation = str_replace(annotation, "^.+ - ", ""))


panvaR::plot_gene_annotation(annotation.table = anno_cleanup,
                             middle.snp = "9-43025000",
                             window = 100,
                             include.id = T)
