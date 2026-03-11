

# test <- data.table::fread("~/scratch/setaria_vcf_chrsnumeric/setaria_annotated_chrsnumeric_CHR9.vcf.gz", nrows = 1000)
library(dplyr)

path <-
execute_snpsift2("~/scratch/setaria_vcf_chrsnumeric/setaria_annotated_chrsnumeric_CHR9.vcf.gz")

x <- fread(path)

# sub <- x[1:100,]
# test <-
# sub %>% 
#   mutate(testcol = strsplit(x = `ANN[*].FEATUREID`, split = ","))

dat.long <- x %>% 
  mutate(across(starts_with("ANN"), ~ strsplit(.x, split = ","))) %>% 
  tidyr::unnest_longer(col = starts_with("ANN"))
names(dat.long) <- stringr::str_replace_all(names(dat.long), "ANN\\[\\*\\]\\.", "")

# get some stats 
info <- dat.long %>% 
  mutate(marker.ID = paste(CHROM, POS, sep = "-")) %>% 
  group_by(marker.ID)

# how many genes per marker?
gene.num <- info %>% 
  summarize(num.genes = length(unique(GENE))) 
hist(gene.num$num.genes)
mean(gene.num$num.genes)
# same but only those labeled as 'protein coding'
gene.num.prot <- info %>% 
  filter(BIOTYPE == "protein_coding") %>% 
  summarize(num.genes = length(unique(GENE))) 
hist(gene.num.prot$num.genes)
mean(gene.num.prot$num.genes)
# number effects per marker
effect.num <- info %>% 
  summarize(num.eff = length(unique(EFFECT)))
hist(effect.num$num.eff)

# plot what effects and impacts are in the set
plot.df <- info %>% 
  select(marker.ID, EFFECT, IMPACT, BIOTYPE) %>% 
  distinct() %>% 
  mutate(dummy = "A")
ggplot(aes(x = dummy, fill = EFFECT), data = plot.df) + 
  geom_bar() +
  facet_wrap(BIOTYPE ~IMPACT, scales = "free_y")

plot.df2 <- plot.df %>% 
  filter(IMPACT == "MODIFIER") 
ggplot(aes(x = BIOTYPE, fill = EFFECT), data = plot.df2) +
  geom_bar() +
  labs(title = "Modifier only",
       x = "protein coding or no?")

# based on this we can see that modifier splits into intergenic vs upstream/downstream variants
# I think it makes sense to use this extra information to split this category
# First, how many markers that are only modifier have intergenic and upstream/downstream designations?
plot.df <- info %>% 
  mutate(only.modifier = all(IMPACT == "MODIFIER")) %>% 
  filter(only.modifier) %>% 
  mutate(num.biotypes = length(unique(BIOTYPE))) %>% 
  filter(num.biotypes == 1) %>% 
  mutate(dummy = "A")
# in only modifier markers with just 1 
ggplot(aes(x = BIOTYPE, fill = EFFECT), data = plot.df) +
  geom_bar() +
  labs(title = "All Modifier sites, 1 biotype")

# some genes get designated as ONLY upstream/downstream 


# make a decision tree to choose the most useful impact
sort <- info %>% 
  ungroup() %>% 
  # First make a new impact category
  mutate(IMPACT_PLUS = case_when(IMPACT == "MODIFIER" & BIOTYPE == "protein_coding" ~ "MODIFIER_CODING",
                                 IMPACT == "MODIFIER" & BIOTYPE == "." ~ "MODIFIER_INTERGENIC",
                                 TRUE ~ IMPACT)) 

# do any markers have multiple low/med/high?
test <- sort %>% 
  filter(IMPACT_PLUS %in% c("LOW", "MODERATE", "HIGH")) %>% 
  group_by(marker.ID) %>% 
  mutate(num.impacts = length(unique(IMPACT_PLUS))) %>% 
  filter(num.impacts > 1)
max(test$num.impacts) 
# yeah some do

# lets make a numeric scale to rank our preference for choosing these impacts
scoring.key <- data.frame(IMPACT_PLUS = c("MODIFIER_INTERGENIC", "MODIFIER_CODING", "LOW", "MODERATE", "HIGH"),
           IMPACT_score = c(1:5))

sort_score <- sort %>% 
  left_join(scoring.key, by = "IMPACT_PLUS") %>% 
  group_by(marker.ID) %>% 
  mutate(is.max.impact.score = IMPACT_score == max(IMPACT_score)) %>% 
  filter(is.max.impact.score)

# still have multiple rows per marker
test <- sort_score %>% 
  mutate(n = n()) %>% 
  filter(n > 1)
# multiple transcripts
test <- sort_score %>% 
  ungroup() %>% 
  select(-FEATUREID) %>% 
  distinct() %>% 
  group_by(marker.ID) %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>% 
  arrange(desc(n))
max(table(test$marker.ID))

# lets keep 2 tables
# table with 1 impact (row) per marker
# table with 1 row per transcript
final_snp_impact_grades <- sort_score %>% 
  select(CHROM, POS, marker.ID, IMPACT_PLUS, IMPACT) %>% 
  distinct() 
max(table(final_snp_impact_grades$marker.ID))
names(table(finalsnp))

# can't retain info about which gene this max impact is because there are ties!
# :( lots of snps as both "upstream" and "downstream" variants. What do we do with these?????? 
test <- sort_score %>% 
  select(CHROM, POS, GENE, marker.ID, IMPACT_PLUS, IMPACT, IMPACT_score) %>% 
  distinct() 
max(table(test$marker.ID))
these <- names(table(test$marker.ID))[which(table(test$marker.ID) > 1)]
# this is really slow... ? 
test2 <- test %>% 
  filter(marker.ID %in% sample(these, 100)) %>% 
  arrange(marker.ID)
# to-do output one row per gene with most impactful impact listed
x <- test %>% 
  group_by(GENE) %>% 
  mutate(is.max.impact.score = IMPACT_score == max(IMPACT_score)) %>% 
  filter(is.max.impact.score) %>% 
  select(GENE, IMPACT_PLUS, IMPACT, IMPACT_score) %>% 
  distinct()

y <- test %>% 
  mutate(GENE_sep = strsplit(GENE, split = "-"))%>% 
  tidyr::unnest_longer(GENE_sep) %>% 
  group_by(GENE_sep) %>% 
  mutate(is.max.impact.score = IMPACT_score == max(IMPACT_score)) %>% 
  filter(is.max.impact.score) %>% 
  select(GENE = GENE_sep, IMPACT_PLUS, IMPACT, IMPACT_score) %>% 
  distinct()

non.inter.genes <- x$GENE[str_detect(x$GENE, "-", negate = T)]
length(non.inter.genes)
length(unique(y$GENE))
# 3 genes had only intergenic snps! 

formatted_snpeff_annotations <- sort 
  
# make this a function:
out <- 
format_snpeff_annotations("~/scratch/setaria_vcf_chrsnumeric/setaria_annotated_chrsnumeric_CHR9.vcf.gz")

# want to find a way to use this table to make snp to gene correspondence
# big file (3GB)
# 1) filter to only snps you care about (in window)
#     - can we read in only those rows? 
# 2) join your snp based variable to this table
# 3) group by gene and get maximum 
test <- out$formatted_snpeff_annotations

out <- format_snpeff_annotations("~/scratch/setaria_annotated_chrsnumeric.vcf.gz")


# ------------------------------------------------------------------------\
# further testing --------
# ------------------------------------------------------------------------\

path <-
  execute_snpsift("~/scratch/setaria_vcf_chrsnumeric/setaria_annotated_chrsnumeric_CHR9.vcf.gz")

dat.long <- fread(path) %>% 
  mutate(across(starts_with("ANN"), ~ strsplit(.x, split = ","))) %>% 
  tidyr::unnest_longer(col = starts_with("ANN"))
names(dat.long) <- stringr::str_replace_all(names(dat.long), "ANN\\[\\*\\]\\.", "")

# make a new impact category that splits modifier
sort <- dat.long %>% 
  ungroup() %>% 
  mutate(marker.ID = paste(.data$CHROM, .data$POS, sep = "-")) %>% 
  mutate(IMPACT_PLUS = case_when(.data$IMPACT == "MODIFIER" & .data$BIOTYPE == "protein_coding" ~ "MODIFIER_CODING",
                                 .data$IMPACT == "MODIFIER" & .data$BIOTYPE == "." ~ "MODIFIER_INTERGENIC",
                                 TRUE ~ .data$IMPACT)) 

#  make a numeric scale to rank our preference for choosing these impacts if a snp has multiple
scoring.key <- data.frame(IMPACT_PLUS = c("MODIFIER_INTERGENIC", "MODIFIER_CODING", "LOW", "MODERATE", "HIGH"),
                          IMPACT_score = c(1:5))

# retain our favorite impact
sort_score <- sort %>% 
  left_join(scoring.key, by = "IMPACT_PLUS") %>% 
  group_by(.data$marker.ID) %>% 
  mutate(is.max.impact.score = .data$IMPACT_score == max(.data$IMPACT_score)) %>% 
  filter(.data$is.max.impact.score) 

a <- sort_score %>% 
  select("CHROM", "POS", "marker.ID", "IMPACT", "IMPACT_PLUS") %>% 
  distinct() 

b <- sort_score %>% 
  mutate(Genes.impacted = paste0(unique(GENE), collapse = "|")) %>% 
  select("Genes.impacted", "CHROM", "POS", "marker.ID", "IMPACT", "IMPACT_PLUS") %>% 
  distinct()
