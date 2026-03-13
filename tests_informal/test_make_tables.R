# make_panvar_tables <- function(gwas.res,
#                                qtl.df = NULL,
#                                tag.snp = NULL,
#                                annotation.table,
#                                plink.path,
#                                pvals.in.log = T,
#                                geno.bed.filename,
#                                geno.bed.directory = "/.",
#                                temp.dir = tempdir(),
#                                window,
#                                compute.scores = F,
#                                score.vars = NULL,
#                                score.dirs = NULL,
#                                score.weights = NULL)

snpeffann <- read.csv("~/scratch/panvar_test/setaria_annotated_chrsnumeric_snpeffinfo_1persnp.csv") %>% 
    select(-CHROM, -POS)

# set global options 
set_plink_path("/home/cluebbert/bin/plink2")
set_panvar_prefix("SetShattering")
set_out_dir("~/scratch/panvar_test/setaria_shatter_full")

gwas.df <- data.table::fread(file.path(options()$panvar_outdir, "SetShattering_GLM_GWASresults.csv")) %>% 
  left_join(snpeffann, by = "marker.ID")
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

impact_score <- data.frame(IMPACT_PLUS = c("MODIFIER_INTERGENIC", "MODIFIER_CODING", "LOW", "MODERATE", "HIGH"),
                           impact_score = c(1:5))

gwas.df <- gwas.df %>% 
  left_join(impact_score, by = "IMPACT_PLUS")
  

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
                             snp.to.gene.vars = c("LD", "snp.score", "impact_score"),
                             snp.to.gene.buffer = 5)


head(tables$anno)
head(tables$gwas)

x <- tables$gwas %>% 
  filter(str_detect(genes_near_snp, "Sevir.5G085400"))

out.anno <- tables$anno

# ------------------------------------------------------------------------\
# manhattan with tables --------
# ------------------------------------------------------------------------\
make_panvar_manhattan(gwas.res = gwas.df,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      window = 250,
                      sig.line = 6, orient = "V")

make_panvar_manhattan(panvar.table.list = tables,
                       pvals.in.log = F,
                       window = 25,
                       sig.line = 6,
                       orient = "H")

# ------------------------------------------------------------------------\
# annotation with tables --------
# ------------------------------------------------------------------------\


make_gene_annotation_plot(panvar.table.list = tables,
                          window = 20, 
                          include.id = T,
                          use.arrows = F)

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
# test snp to gene --------
# ------------------------------------------------------------------------\

gwas.df <- data.table::fread(file.path(options()$panvar_outdir, "SetShattering_GLM_GWASresults.csv")) %>% 
  left_join(snpeffann, by = "marker.ID")
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

gene <- "Sevir.5G085400"
gene.info <- anno %>% 
  filter(geneID == gene)

snps.in.gene <- gwas.df %>% 
  filter(CHR == 5) %>% 
  filter(between(POS, gene.info$start, gene.info$end))

gwas.sub <- gwas.df %>% 
  filter(CHR == 5) %>% 
  filter(between(POS, gene.info$start - 25e4, gene.info$end + 25e4))

snp.to.gene.buffer <- 0

gwas.sub_with.genes <- gwas.sub %>% 
  rowwise() %>% 
  mutate(snp.in.gene_list = get.gene.from.snp(.data$POS, anno, snp.to.gene.buffer)) %>% 
  mutate(snp.in.gene = paste0(snp.in.gene_list, collapse = "|")) 

point.color.stat <- gwas.sub_with.genes %>% 
  filter(!is.null(.data$snp.in.gene_list)) %>% 
  unnest_longer(.data$snp.in.gene_list) %>% 
  group_by(.data$snp.in.gene_list) %>% 
  # summarize(maximum.value = max(.data[[annotation.point.variable]])) %>% 
  summarize(across(all_of("LOGPVAL"), ~ max(.x, na.rm = T)))  %>% 
  rename("geneID" = "snp.in.gene_list")

x <- gwas.sub_with.genes %>% 
  select(-snp.in.gene_list) %>%
  rename("genes_near_snp" = "snp.in.gene")
  

# temp 
point.color.stat <- gwas.sub %>% 
  rowwise() %>% 
  mutate(snp.in.gene = get.gene.from.snp(.data$POS, anno.sub, snp.to.gene.buffer)) %>% 
  filter(!is.null(.data$snp.in.gene)) %>% 
  unnest_longer(.data$snp.in.gene) %>% 
  group_by(.data$snp.in.gene) %>% 
  # summarize(maximum.value = max(.data[[annotation.point.variable]])) %>% 
  summarize(across(all_of(snp.to.gene.vars), ~ max(.x, na.rm = T)))  %>% 
  rename("geneID" = "snp.in.gene")
