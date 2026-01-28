

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

formatted_snpeff_annotations <- sort 
  
# make this a function:
out <- 
format_snpeff_annotations("~/scratch/setaria_vcf_chrsnumeric/setaria_annotated_chrsnumeric_CHR9.vcf.gz")

out <- format_snpeff_annotations("~/scratch/setaria_annotated_chrsnumeric.vcf.gz")
