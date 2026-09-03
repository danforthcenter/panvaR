

library(tidyverse)

clump.glm.res <- read.csv("C:/Users/cluebbert/OneDrive - DDPSC/~Rprojects~/PGRP_both_pops/results_gwas/testing/TEST_bg_clump1_glm.res.csv")
this.clump <- read.csv("C:/Users/cluebbert/OneDrive - DDPSC/~Rprojects~/PGRP_both_pops/results_gwas/testing/TEST_bg_clump1_qtl.df.csv")
annotation.table <- read.csv("C:/Users/cluebbert/OneDrive - DDPSC/~Rprojects~/PGRP_both_pops/results_gwas/testing/ZeaMaysB73_V5_allgenes_Uniprot_descriptions_with_coordinates.csv")%>% 
  mutate(annotation = case_when(GeneNameShort != "N/A" ~ paste0(GeneNameShort, ", ", Description),
                                TRUE ~ Description)) %>% 
  mutate(CHR = as.numeric(CHROM)) %>% 
  rename(geneID = GeneID)

clump.glm.res <- clump.glm.res %>% 
  rename(EFF = Effect)
this.clump <- this.clump %>% 
  rename(EFF = Effect)

panvar.res <-
  panvaR::make_panvar_tables(gwas.res = clump.glm.res,
                             qtl.df = this.clump,
                             pvals.in.log = F,
                             annotation.table = annotation.table,
                             geno.bed.filename = "1.BG_schnable23_filtered_maf0.05.maxMissing0.1_fixnames",
                             geno.bed.directory = "C:/Users/cluebbert/OneDrive - DDPSC/~Rprojects~/PGRP_both_pops/results_gwas/testing/",
                             window = 500)


trial.scale <- make_consistent_scale(values = c(15, 16, 17), vars = c("AllExps", "PGRP18", "PGRP19"), type = "shape")


plot_panvar_manhattan(panvar.res,
                      window = 500,
                      sig.line = 6,
                      pvals.in.log = F,
                      point.fill.variable.d = "trait",
                      snp.highlight.df = this.clump,
                      snp.highlight.shape.var = "trial",
                      snp.highlight.shape.scale = trial.scale,
                      snp.highlight.color.var = "trait_only",
                      snp.highlight.color.scale = khroma::scale_color_bright())

plot_panvar_manhattan(panvar.res,
                      window = 500,
                      sig.line = 6,
                      pvals.in.log = F,
                      point.shape.variable = "trait",
                      point.fill.variable.d = "trait")

# Might want to try doing both ways, by trait is interesting as it shows effect direction and could indicate 2 hits 
# color by trait
plot_panvar(panvar.res,
            window = 500,
            sig.line = 6, 
            pvals.in.log = F,
            plot.effect = T,
            point.fill.variable.d = "trait",
            snp.highlight.df = this.clump,
            snp.highlight.shape.var = "trial",
            snp.highlight.shape.scale = trial.scale,
            snp.highlight.color.var = "trait_only",
            snp.highlight.color.scale = khroma::scale_color_bright())

# color by ld
plot_panvar(panvar.res,
            window = 500,
            sig.line = 6, 
            pvals.in.log = F,
            plot.effect = T,
            # point.fill.variable.d = "trait",
            snp.highlight.df = this.clump,
            snp.highlight.shape.var = "trial",
            snp.highlight.shape.scale = trial.scale,
            snp.highlight.color.var = "trait_only",
            snp.highlight.color.scale = khroma::scale_color_bright())
