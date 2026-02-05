

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

make_gene_annotation_plot(anno_ld,
                          middle.snp = out_ld$key.snp,
                          window = 25, 
                          include.id = T,
                          use.arrows = F,
                          point.color = "ld_sim",
                          point.fill.scale = NULL)

# ------------------------------------------------------------------------\
# horizontal --------
# ------------------------------------------------------------------------\

make_gene_annotation_plot_horiz(anno,
                          middle.snp = out_ld$key.snp,
                          window = 25, 
                          include.id = T,
                          use.arrows = F)

