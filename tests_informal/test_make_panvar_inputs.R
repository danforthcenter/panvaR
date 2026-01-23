


pheno <- read.delim("inst/extdata/Setaria_shattering_example_phenotype.tsv")

make_panvar_inputs(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
                   phenotype.table = pheno,
                   out.prefix = "INPUTSPHENO",
                   out.dir = "~/scratch/temp")

pheno <- read.csv("~/scratch/gwas_in_F_pgrp/6.BG.bothftcovs_all_experiments_BLUPs_scaled_withPCtraits.csv")
pheno <- pheno %>% 
  filter(FT_cov == "noFT") %>% 
  select(Genotype, GDU_to_Anthesis_50_AllExps)

make_panvar_inputs(genotype.path = "~/scratch/gwas_in_F_pgrp/1.WiDiv_942g_AGPv4_imputed_maf0.05_maxmaf0.95_masMissing0.1_onlyphenotypedlines.bed",
                   phenotype.table = pheno,
                   min.maf = 0,
                   max.missing.snp = 1)
