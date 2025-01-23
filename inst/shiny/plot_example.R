# library(panvaR)
library(tidyverse)
library(luebbert)
library(scales)


bonf.cut <- -log10(.05 / 5566446) 
gene.snp <- 6849363


# gene locations
geno.loc <- read.delim("data/setaria_gene_locations.tsv") 

# panvar results
results <- read.csv("data/setaria_shattering_0r2_allimpacts_2MB_highestsnp.csv")

# get only high ld ones
res.ld.sub <- results %>% 
  filter(LD >= .5) 

# region start/end
qtl.start <- min(res.ld.sub$BP)
qtl.end <- max(res.ld.sub$BP)

# get all regardless of significance in region
res.ld.sub_plot <- results %>% 
  filter(BP > qtl.start, BP < qtl.end)


min(res.ld.sub_plot$BP)

# only impactful ones
res.sub <- results %>% 
  filter(IMPACT != "tag_snp") %>% 
  filter(IMPACT != "MODIFIER") %>% 
  filter(IMPACT != "LOW")

# only ones over ld threshold
res.sub_over5 <- res.sub %>% 
  filter(LD > .5) %>% 
  filter(Pvalues > bonf.cut)

tag_df <- results %>% 
  filter(IMPACT == "tag_snp")

# color palette
pal <- khroma::color("incandescent")
my.colors <- pal(5)

# plotting df for genes
# geno.loc_res.sub <- geno.loc %>% 
#   filter(Chrom == unique(results$CHROM)) %>% 
#   filter(Ext_Start >= min(res.sub$BP) & Ext_Stop <= max(res.sub$BP)) %>% 
#   filter(Gene %in% temp$gene_only) %>% 
#   rowwise() %>% 
#   mutate(mid.point = median(c(Ext_Start, Ext_Stop)))


res.sub_over5 <- res.sub %>% 
  filter(LD > .5) %>% 
  filter(Pvalues > bonf.cut)


# make plot
ggplot(aes(x = BP, y = Pvalues), data = res.sub_over5) + 
  geom_point(aes(shape = IMPACT, fill = LD),
             size = 3,
             color = "black") + 
  #geom_point(data = tag_df, color = "black", size = 4) +
  geom_vline(data = tag_df, aes(xintercept = BP), linewidth = 1.5) +
  geom_vline(xintercept = gene.snp, color = "red", linetype = "dashed", linewidth = 1.5) +
  geom_hline(aes(yintercept = bonf.cut), color = "black", linetype = "dashed") +
  scale_shape_manual(values = c(25, 21, 22)) +
  scale_fill_stepsn(colors = my.colors, name = "R2") +
  theme_bw() +
  theme(text = element_text(size = 16),
        legend.position = "left",
        legend.text = element_text(size = 10)) +
  scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = "Position", 
       y = bquote(-log[10](p-value))) 
  # geom_segment(aes(x = Ext_Start, xend = Ext_Stop, y = max(res.ld.sub_plot$Pvalues) + 1),
  #              linewidth = 10,
  #              color = "grey44",
  #              data = geno.loc_res.sub) +
  # geom_point(aes(x = mid.point, y = max(res.ld.sub_plot$Pvalues) + 1),
  #            color = "grey22",
  #            data = geno.loc_res.sub) 
