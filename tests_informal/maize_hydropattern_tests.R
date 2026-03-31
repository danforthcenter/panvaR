
# ------------------------------------------------------------------------\
# testing scores error --------
# ------------------------------------------------------------------------\

# had an error where NA's were leading to all NA when making scores,
# fixed this in the min/max function by adding na.rm
setwd("~/Projects/Hydropattern_finemap/")


window <- 50

# prepare annotation
annotation.table <- read.csv("data/ZeaMaysB73_V5_allgenes_Uniprot_descriptions_with_coordinates.csv") %>% 
  mutate(annotation = case_when(GeneNameShort != "N/A" ~ paste0(GeneNameShort, ", ", Description),
                                TRUE ~ Description)) %>% 
  mutate(CHR = as.numeric(CHROM)) %>% 
  rename(geneID = GeneID)

# prepare gwas results
clumps.df <- read.csv("results/1.Paper_sig_hits_clumped_by_closest_marker.csv") 

clump.phenos <- clumps.df %>% 
  select(clump_num, Trait) %>% 
  distinct()

bed.file <- "hydropattern_Schanble2023__PlinkQC_maf0.05_missing0.1.bed"
bed.dir <- "results/shared_snp"


all.gwas.df <- read.csv("results/clump1_gwasres_fortesting.csv")

tables <- 
make_panvar_tables(gwas.res = all.gwas.df,
                   # qtl.df = this.qtl.df,
                   tag.snp = "5-209562356",
                   annotation.table = annotation.table,
                   pvals.in.log = F,
                   plink.path = "plink2",
                   geno.bed.filename = bed.file,
                   geno.bed.directory = bed.dir,
                   window = window,
                   snp.to.gene.vars = c("neg_PlantCad_zeroshot", "snp.score", "LD","LOGPVAL"),
                   snp.to.gene.buffer = 5,
                   compute.scores = T,
                   score.vars = c("LD", "LOGPVAL", "neg_PlantCad_zeroshot"),
                   score.dirs = c(1, 1, 1))

head(tables$gwas)
x <- tables$anno


# ------------------------------------------------------------------------\
# annotation error --------
# ------------------------------------------------------------------------\

# having an issue with this table plotting gene annotations
tables <- readRDS("../Hydropattern_finemap/results/clump5_tables_fortesting.rds")

plot_gene_annotation(tables,
                     window = 50,
                     point.color = "LOGPVAL")


# turned out this was some weird behavior with geom_fit_text
# the y height of the text bounding box was getting automatically set strangely.
# didn't figure out exactly how. 
# Fixed by setting a "height" parameter it geom_fit_text if there is just one row
# and null if more than one to allow the program to automatically pick these when 
# its not broken. 

seg.df <- data.frame(ystart = 0, yend = 1, label = tables$anno$annotation)

ggplot(data = seg.df) +
  geom_segment(aes(y = ystart, yend = yend, x = 1)) +
  geom_fit_text(aes(xmin = 1.5, xmax = 2, y = .5, label = label),
                place = "left",
                #grow = TRUE,
                hjust = 0,
                padding.y = grid::unit(.1, "lines"),
                min.size = 4)


anno +
  # geom_segment(aes(x = text.x.start, xend = .85, y = .data$y.pos), color = "red") +
  ggfittext::geom_fit_text(aes(xmin = text.x.start, xmax = .85, 
                               ymin = .data$y.pos - text.box.width.half, ymax = .data$y.pos + text.box.width.half, 
                               label = .data$plot.label),
                           min.size = 4,
                           place = "left",
                           hjust = 0) 