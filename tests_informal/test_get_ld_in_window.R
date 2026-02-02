# still haven't fixed the non numeric chromsome ids in the example data so this is the id
# "Chr_05-6857045"
# Setaria_shattering_example_pruned

set_plink_path("/home/cluebbert/bin/plink2")

# test tag snp
out_ld <- 
panvaR::get_ld_in_window(tag.snp = "Chr_05-6857045",
                 window = 500*2,
                 in.dir = "inst/extdata/",
                 out.dir = "~/scratch/panvar_test/",
                 geno.bed = "Setaria_shattering_example_pruned")

head(out_ld)

# test qtl.df
qtl.df <- read.csv("~/scratch/gwas_in_F_pgrp/6F.gwas_res_above.eff.bonf_with.clump_noFT.model_clean.csv")
sub <- qtl.df %>% 
  filter(trait == "GDU_to_Anthesis_50") %>% 
  arrange(desc(pval))
sub <- qtl.df %>% 
  filter(clump_num == 1156) %>% 
  rename(CHR = chromosome, 
         POS = physical.pos,
         LOGPVAL = logpval)

panvaR::get_ld_in_window(qtl.df = sub,
                         window = 500*2,
                         in.dir = "~/scratch/gwas_in_F_pgrp/",
                         out.dir = "~/scratch/gwas_in_F_pgrp/",
                         geno.bed = "1.WiDiv_942g_AGPv4_imputed_maf0.05_maxmaf0.95_masMissing0.1_onlyphenotypedlines_fixnames", 
                         verbose = T)


# ------------------------------------------------------------------------\
# sorghum --------
# ------------------------------------------------------------------------\

set_out_dir("~/Projects/Sorghum13CMash/results/")
set_panvar_prefix("Sorg13C")
set_plink_path("~/bin/plink2")

ldlist <- get_ld_in_window(tag.snp = "1-21909576",
                           window = 500,
                           geno.bed = "Sorg13C_PlinkQC_maf0.05_missing0.1.bed",
                           in.dir = "~/Projects/Sorghum13CMash/results/")

gwas.res <- fread("~/Projects/Sorghum13CMash/results/Sorg13C_GLM_GWASresults.csv", data.table = F)
make_panvar_manhattan(gwas.res = gwas.res,
                      ld.list = ldlist,
                      window = 500,
                      sig.line = 8,
                      pvals.in.log = F)
