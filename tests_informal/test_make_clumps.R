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
        
