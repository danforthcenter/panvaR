

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

make_gene_annotation_plot(anno,
                          middle.snp = out_ld$key.snp,
                          window = 25, 
                          include.id = T,
                          use.arrows = F)

# ------------------------------------------------------------------------\
# horizontal --------
# ------------------------------------------------------------------------\

make_gene_annotation_plot_horiz(anno,
                          middle.snp = out_ld$key.snp,
                          window = 25, 
                          include.id = T,
                          use.arrows = F)

