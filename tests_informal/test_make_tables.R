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

gwas.df <- data.table::fread(file.path(options()$panvar_outdir, "SetShattering_GLM_GWASresults.csv"))
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
                             snp.to.gene.vars = c("LD", "snp.score"))


head(tables$anno)
head(tables$gwas)


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
                       window = 250,
                       sig.line = 6,
                       orient = "H")
