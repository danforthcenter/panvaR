setwd("~/Projects/Sorghum13CMash/")

# tagsnps
tagsnps <- read.csv("data/SorghumD13_sig_markers.csv")
this.tag.snp <- tagsnps$marker.ID[tagsnps$name == "sig1"]

gwas.df <- data.table::fread("results/panvar_res/sig1/sig1_Sorg13C_GLM_GWASresults.csv", data.table = F)

# prepare annotation, chose arabi as annotation for now
annotation.table <- read.csv("data/SorghumV3_annotation.csv") %>% 
  mutate(CHR = as.numeric(str_replace(CHR, "Chr", ""))) %>% 
  mutate(annotation = arabi.anno)

test <- gwas.df %>%
  filter(is.na(PVAL))
# luebbert::ritecsv(test, "results/sig1_NA_pvalues_fromMVPGLMgwas.csv")
# these all are strange snps, lets remove them
gwas.df <- gwas.df %>%
  filter(!is.na(PVAL))

# add cad scores now
cad_score <- read.table("data/plantCAD_Crawford_BAP_d13Csighits_2MBwindow.txt", header = F)
names(cad_score) <- c("CHR", "POS", "PlantCad_zeroshot")
cad_score <- cad_score %>%
  mutate(CHR = as.numeric(str_replace(CHR, "Chr", ""))) %>%
  mutate(marker.ID = paste0(CHR, "-", POS)) %>%
  mutate(PlantCad_zeroshot = str_replace(PlantCad_zeroshot, "plantCAD_zero_shot=", "")) %>%
  # relocate(marker.ID) %>%
  select(-CHR, -POS) %>% 
  mutate(PlantCad_zeroshot = as.numeric(PlantCad_zeroshot))

gwas.df_cad <- gwas.df %>%
  left_join(cad_score, by = "marker.ID")

tables <-
  make_panvar_tables(gwas.res = gwas.df_cad,
                     tag.snp = this.tag.snp,
                     annotation.table = annotation.table,
                     pvals.in.log = F,
                     plink.path = NULL,
                     geno.bed.filename = "sig1_Sorg13C_PlinkQC_maf0.05_missing0.1.bed",
                     geno.bed.directory = "results/panvar_res/sig1",
                     window = 500,
                     snp.to.gene.vars = c("PlantCad_zeroshot", "snp.score", "LD"),
                     snp.to.gene.buffer = 5,
                     compute.scores = T,
                     score.vars = c("LD", "LOGPVAL", "PlantCad_zeroshot"),
                     score.dirs = c(1, 1, -1))

# get bonferroni

plot_panvar(tables,
            window = 50,
            sig.line = 6,
            pvals.in.log = F,
            annotation.point.variable = "PlantCad_zeroshot",
            plot.effect = T)


# ------------------------------------------------------------------------\
# testing with clumps --------
# ------------------------------------------------------------------------\

setwd("~/Projects/Sorghum13CMash/")

# tagsnps
# tagsnps <- read.csv("data/SorghumD13_sig_markers.csv")
# this.tag.snp <- tagsnps$marker.ID[tagsnps$name == "sig1"]
# try to use a clump 
clumps <- read.csv("results/1.SorghumD13C_MashBayesFactors_clumped.csv") %>% 
  filter(clump_num == 2)


gwas.df <- data.table::fread("results/panvar_res/sig1/sig1_Sorg13C_GLM_GWASresults.csv", data.table = F)

this.qtl.df <- gwas.df %>% 
  filter(marker.ID %in% clumps$marker.ID)

# prepare annotation, chose arabi as annotation for now
annotation.table <- read.csv("data/SorghumV3_annotation.csv") %>% 
  mutate(CHR = as.numeric(str_replace(CHR, "Chr", ""))) %>% 
  mutate(annotation = arabi.anno)

test <- gwas.df %>%
  filter(is.na(PVAL))
# luebbert::ritecsv(test, "results/sig1_NA_pvalues_fromMVPGLMgwas.csv")
# these all are strange snps, lets remove them
gwas.df <- gwas.df %>%
  filter(!is.na(PVAL))

# add cad scores now
cad_score <- read.table("data/plantCAD_Crawford_BAP_d13Csighits_2MBwindow.txt", header = F)
names(cad_score) <- c("CHR", "POS", "PlantCad_zeroshot")
cad_score <- cad_score %>%
  mutate(CHR = as.numeric(str_replace(CHR, "Chr", ""))) %>%
  mutate(marker.ID = paste0(CHR, "-", POS)) %>%
  mutate(PlantCad_zeroshot = str_replace(PlantCad_zeroshot, "plantCAD_zero_shot=", "")) %>%
  # relocate(marker.ID) %>%
  select(-CHR, -POS) %>% 
  mutate(PlantCad_zeroshot = as.numeric(PlantCad_zeroshot))

gwas.df_cad <- gwas.df %>%
  left_join(cad_score, by = "marker.ID")

tables <-
  make_panvar_tables(gwas.res = gwas.df_cad,
                     # tag.snp = this.tag.snp,
                     qtl.df = this.qtl.df,
                     annotation.table = annotation.table,
                     pvals.in.log = F,
                     plink.path = NULL,
                     geno.bed.filename = "sig1_Sorg13C_PlinkQC_maf0.05_missing0.1.bed",
                     geno.bed.directory = "results/panvar_res/sig1",
                     window = 500,
                     snp.to.gene.vars = c("PlantCad_zeroshot", "snp.score", "LD"),
                     snp.to.gene.buffer = 5,
                     compute.scores = T,
                     score.vars = c("LD", "LOGPVAL", "PlantCad_zeroshot"),
                     score.dirs = c(1, 1, -1))



plot_panvar_manhattan(panvar.table.list = tables, 
                      pvals.in.log = F,
                      plot.r2.thresh = .2,
                      window = 75,
                      sig.line = cutoff.line)

plot_panvar(tables,
            window = 200,
            pvals.in.log = F,
            sig.line = 6,
            annotation.point.variable = "PlantCad_zeroshot",
            annotation.point.scale = scale_color_viridis_b(option = "inferno"),
            plot.effect = T,
            include.gene.id = T)


# ------------------------------------------------------------------------\
# look at clump1 --------
# ------------------------------------------------------------------------\

library(tidyverse)
setwd("~/Projects/Sorghum13CMash/")
# prepare annotation, choose arabi as annotation for now
annotation.table <- read.csv("data/Sorghum_V3_Annotation_best_gene_descriptions.csv") %>%
  rename(annotation = BestDescription)

# set up inputs
top.gwas.phenos <- read.csv("results/2.top_gwas_phenos.csv")

clumps.df <- read.csv("results/1.SorghumD13C_MashBayesFactors_clumped.csv")

loop.df <- top.gwas.phenos %>%
  rowwise() %>%
  mutate(gwas.file.path = list.files(path = paste0("results/panvar_res/", pheno),
                                     pattern = "GWASresults",
                                     full.names = T)[1],
         bed.file.name = list.files(path = paste0("results/panvar_res/", pheno),
                                    pattern = "bed$"),
         bed.file.directory = paste0("results/panvar_res/", pheno)) %>%
  filter(clump_num <= 10)

temp <- clumps.df %>%
  group_by(clump_num) %>%
  mutate(key.snp = marker.ID[which(bf == max(bf))]) %>%
  select(clump_num, key.snp) %>%
  distinct()

loop.df <- loop.df %>%
  left_join(temp, by = "clump_num")

bonf <- -log10(.05/ 4837900)
cutoff.line <- bonf
window <- 100

i <- 1

this.clump.num <- loop.df$clump_num[i]
this.gwas.path <- loop.df$gwas.file.path[i]
this.bed.file <- loop.df$bed.file.name[i]
this.bed.dir <- loop.df$bed.file.directory[i]
this.outdir <- paste0("results/panvar_res/", "Clump", this.clump.num)
this.key.snp <- loop.df$key.snp[i]

set_out_dir(this.outdir)
# set_panvar_prefix(paste0(this.dataset, "_Sorg13C"))

# tagsnps
# this.tag.snp <- tagsnps$marker.ID[tagsnps$name == this.dataset]

# gwas.df <- data.table::fread("results/panvar_res/sig1/sig1_Sorg13C_GLM_GWASresults.csv", data.table = F)
gwas.df <- data.table::fread(this.gwas.path, data.table = F)

gwas.df <- gwas.df %>%
  filter(!is.na(PVAL))

# add cad scores now
cad_score <- read.table("data/plantCAD_Crawford_BAP_top10BFclumps_2MBwindow.txt", header = F)
names(cad_score) <- c("CHR", "POS", "PlantCad_zeroshot")
cad_score <- cad_score %>%
  mutate(CHR = as.numeric(str_replace(CHR, "Chr", ""))) %>%
  mutate(marker.ID = paste0(CHR, "-", POS)) %>%
  mutate(PlantCad_zeroshot = str_replace(PlantCad_zeroshot, "plantCAD_zero_shot=", "")) %>%
  relocate(marker.ID) %>%
  select(-CHR, -POS) %>%
  mutate(PlantCad_zeroshot = as.numeric(PlantCad_zeroshot)) %>%
  mutate(neg_PlantCad_zeroshot = -PlantCad_zeroshot)

gwas.df_cad <- gwas.df %>%
  left_join(cad_score, by = "marker.ID")

this.clump.markers <- clumps.df %>%
  filter(clump_num == this.clump.num) %>%
  pull(marker.ID)

this.qtl.df <- gwas.df_cad %>%
  filter(marker.ID %in% this.clump.markers)

tables <-
  make_panvar_tables(gwas.res = gwas.df_cad,
                     qtl.df = this.qtl.df,
                     annotation.table = annotation.table,
                     pvals.in.log = F,
                     plink.path = "plink2",
                     geno.bed.filename = this.bed.file,
                     geno.bed.directory = this.bed.dir,
                     window = window,
                     snp.to.gene.vars = c("neg_PlantCad_zeroshot", "snp.score", "LD","LOGPVAL"),
                     snp.to.gene.buffer = 5,
                     compute.scores = T,
                     score.vars = c("LD", "LOGPVAL", "PlantCad_zeroshot"),
                     score.dirs = c(1, 1, -1))

this.qtl.df$marker.ID %in% tables$gwas$marker.ID

p <-
  plot_panvar(tables,
              window = window,
              pvals.in.log = F,
              sig.line = cutoff.line,
              unplotted.alpha = .1,
              annotation.point.variable = "neg_PlantCad_zeroshot",
              annotation.point.scale = scale_color_viridis_b(option = "inferno"),
              plot.effect = T,
              include.gene.id = T,
              plot.title = this.key.snp)

x <- filter(tables$gwas, marker.ID == this.key.snp)
x

memesave(myfilename = file.path(this.outdir, paste0("Clump", this.clump.num, "_PanvarPlot_bestdescription.png")),
         plot = p,
         width = 20)

saveRDS(p, file.path(this.outdir, paste0("Clump", this.clump.num, "_PanvarPlotObject_", window, "KBwindow.rds")))
ritecsv(tables$anno, file.path(this.outdir, paste0("Clump", this.clump.num, "_PanvarGeneTable_", window, "KBwindow.csv")))
ritecsv(tables$gwas, file.path(this.outdir, paste0("Clump", this.clump.num, "_PanvarGWASTable_", window, "KBwindow.csv")))


# ------------------------------------------------------------------------\
# snp2gene --------
# ------------------------------------------------------------------------\

setwd("~/Projects/panvaR")

gwas.sub <- read.csv("tests_informal/files/snp2gene_gwas.sub.csv")
anno.sub <- read.csv("tests_informal/files/snp2gene_anno.sub.csv")
my.gene.anno <- anno.sub %>% 
  filter(geneID == "Sobic.009G144400")
test <- gwas.sub %>% 
  filter(between(POS, 50173181 - 15000, 50173181 + 15000))

# get.gene.from.snp <- function(bp, gene.df, gene.buffer = 0){
#   gene.buffer * 1000
#   check.vec <- data.table::between(bp, gene.df$start - gene.buffer, gene.df$end + gene.buffer)
#   if(any(check.vec)){
#     gene.id.out <- list(gene.df$geneID[check.vec])
#     # gene.id.out <- paste(gene.df$geneID[check.vec], collapse = "|")
#   } else {
#     gene.id.out <- NA
#   }
#   return(gene.id.out)
# }

test.bp <- 50168394
gene.df <- anno.sub
gene.buffer = 5
gene.buffer <- gene.buffer * 1000
check.vec <- data.table::between(test.bp, gene.df$start - gene.buffer, gene.df$end + gene.buffer)
any(check.vec)
gene.id.out <- list(gene.df$geneID[check.vec])
gene.id.out
test.bp <- 50168394

get.gene.from.snp(test.bp, gene.df = anno.sub, gene.buffer = 5)

point.color.stat <- test %>% 
  rowwise() %>% 
  mutate(within.5kb.gene = between(POS, 50173181 - 5000, 50173181 + 5000))
  mutate(snp.in.gene = get.gene.from.snp(.data$POS, anno.sub, 5)) %>% 
  filter(!is.null(.data$snp.in.gene)) %>% 
  unnest_longer(.data$snp.in.gene) %>% 
  group_by(.data$snp.in.gene) %>% 
  # summarize(maximum.value = max(.data[[annotation.point.variable]])) %>% 
  summarize(across(all_of(snp.to.gene.vars), ~ max(.x, na.rm = T)))  %>% 
  rename("geneID" = "snp.in.gene")

anno.out <- anno.sub %>% 
  left_join(point.color.stat, by = "geneID")
