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
                     snp.to.gene.vars = c("PlantCad_zeroshot", "snp.score"),
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
