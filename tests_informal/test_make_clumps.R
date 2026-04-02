set_plink_path("/home/cluebbert/bin/plink2")

bfs <- read.csv("~/Projects/Sorghum13CMash/data/SorghumD13C_MashBayesFactors.csv")
df.in <- bfs %>% 
  mutate(marker.ID = paste0(Chr, "-", Pos)) %>% 
  rename(CHR = Chr,
         POS = Pos,
         PVAL = bf)


snp_make_clumps(geno.bed.filename = "Crawford_BAP_new.recode_chr1to10",
                geno.bed.dir = "~/Projects/Sorghum13CMash/data/",
                gwas.res = df.in,
                window = 250,
                ld.thresh = .6)



make_ld(plink.path = NULL,
        snp.name = "9-50221901",
        window = 250,
        bedfile = "Crawford_BAP_new.recode_chr1to10",
        in.dir = "/home/cluebbert/Projects/Sorghum13CMash/data",
        out.dir = tempdir())
        
# ------------------------------------------------------------------------\
# another dataset --------
# ------------------------------------------------------------------------\

setwd("~/Projects/Sorghum13CMash/")

bf.df_sub <- read.csv("results/setaria_d13C_mash/MashResults_BF_over5.csv")

df.in <- bf.df_sub %>%
  mutate(marker.ID = paste0(Chr, "-", Pos)) %>%
  rename(CHR = Chr,
         POS = Pos,
         PVAL = bf) %>% 
  mutate(CHR = paste0("Chr_0", CHR),
         marker.ID = paste0(CHR, "_", POS))

clumps <-
  panvaR::snp_make_clumps(geno.bed.filename = "8.2.IB007_maf.1.maxmaf.9.hetsIMP.filteredSNPs.hetFilter0.25.recode",
                          geno.bed.dir = "results/setaria_d13C_mash/",
                          gwas.res = df.in,
                          pvals.in.log = T,
                          window = 50,
                          ld.thresh = .6,
                          plink.path = normalizePath("~/bin/plink2"))

x <- data.table::fread("results/setaria_d13C_mash/8.2.IB007_maf.1.maxmaf.9.hetsIMP.filteredSNPs.hetFilter0.25.recode.bim")
head(x)
