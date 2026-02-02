
# set global options 
set_plink_path("/home/cluebbert/bin/plink2")
set_panvar_prefix("SetShattering")
set_out_dir("~/scratch/panvar_test/setaria_shatter_full")

# ------------------------------------------------------------------------\
# 1) Standardize inputs  --------
# ------------------------------------------------------------------------\

pheno <- read.delim("inst/extdata/Setaria_shattering_example_phenotype.tsv")

make_panvar_inputs(genotype.path = "~/scratch/setaria_annotated_chrsnumeric.vcf.gz",
                   phenotype.path = "inst/extdata/Setaria_shattering_example_phenotype.tsv",
                   out.prefix = "SetShattering",
                   out.dir = "~/scratch/panvar_test/setaria_shatter_full")


# ------------------------------------------------------------------------\
# 2) Run gwas --------
# ------------------------------------------------------------------------\


panvar_mvp_gwas(inputs.dir = options()$panvar_outdir,
                npcs = 3,
                gwas.model = "GLM",
                output.manhattan = T)




# For each SNP ------------------------------------------------------------


# ------------------------------------------------------------------------\
# 1) Calculate LD --------
# ------------------------------------------------------------------------\


out_ld <- 
  panvaR::get_ld_in_window(tag.snp = "5-6857045",
                           window = 500,
                           in.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full/",
                           out.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full",
                           geno.bed = "SetShattering_PlinkQC_maf0.05_missing0.1")

# ------------------------------------------------------------------------\
# 2) Plot  --------
# ------------------------------------------------------------------------\

gwas.df <- data.table::fread(file.path(options()$panvar_outdir, "SetShattering_GLM_GWASresults.csv"))

make_panvar_manhattan(gwas.res = gwas.df,
                      qtl.df = NULL,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      window = 500,
                      sig.line = 6, orient = "V")

# add snpeff annotations
snpeffann <- read.csv("~/scratch/panvar_test/setaria_annotated_chrsnumeric_snpeffinfo_1persnp.csv") %>% 
  select(-CHROM, -POS)
gwas.df_ann <- gwas.df %>% 
  left_join(snpeffann, by = "marker.ID")

myshapescale <- make_consistent_scale(values = c(21:25), vars = unique(gwas.df_ann$IMPACT_PLUS), type = "shape", show.example = T)

p <-
make_panvar_manhattan(gwas.res = gwas.df_ann,
                      qtl.df = NULL,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      window = 500,
                      sig.line = 6, 
                      orient = "V",
                      qualitative.annotation = "IMPACT_PLUS",
                      qualitative.shape.scale = NULL)


p
ggplotly(p)

# add something that approximates zero shot score
gwas.df_ann <- gwas.df_ann %>% 
  mutate(zero.shot.sim = -rnorm(nrow(gwas.df_ann), -1, sd = 2))

# p <- 
make_panvar_manhattan(gwas.res = gwas.df_ann,
                      qtl.df = NULL,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      plot.r2.thresh = .6,
                      unplotted.alpha = 0,
                      window = 500,
                      sig.line = 6, 
                      orient = "V",
                      qualitative.annotation = "IMPACT_PLUS",
                      qualitative.shape.scale = NULL,
                      quantitative.annotation = "zero.shot.sim")


cadscore <- arrow::read_parquet("~/scratch/setaria_plantCAD_scores.parquet", col_select = -annotation) %>% 
  mutate(CHR = as.numeric(str_replace(CHR, "Chr_", ""))) %>% 
  mutate(marker.ID = paste0(CHR, "-", POS)) %>% 
  select(marker.ID, zero_shot) %>% 
  mutate(zero_shot = as.numeric(zero_shot))

gwas.df_ann <- left_join(gwas.df_ann, cadscore, by = "marker.ID") %>% 
  mutate(zero_shot_positive = -zero_shot)

make_panvar_manhattan(gwas.res = gwas.df_ann,
                      qtl.df = NULL,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      plot.r2.thresh = .8,
                      unplotted.alpha = 0,
                      window = 50,
                      sig.line = 6, 
                      orient = "V",
                      qualitative.annotation = "IMPACT_PLUS",
                      qualitative.shape.scale = NULL,
                      quantitative.annotation = "zero_shot_positive",
                      quantitative.fill.scale = "magma")




# ------------------------------------------------------------------------\
# test consistent scale --------
# ------------------------------------------------------------------------\

qual.vars <-  c("HIGH", "MODERATE", "LOW", "MODIFIER")
values <- c(24, 22, 25, 23)
vars <- qual.vars
names(values) <- levels(factor(vars, levels = vars))

out <- scale_shape_manual(name = "TEST", values = values)

x <- data.frame(shape = values, label = names(values), x.coord = seq(1, length(values), by = 1))


p <-
  ggplot(aes(x = .data$x.coord, y = 1), data = x) +
  geom_point(aes(shape = x$label), show.legend =  F, size = 20, fill = "orange") +
  # scale_shape_manual(values = x$shape) +
  out +
  geom_text(aes(label = .data$label, y = 2), size = 8) +
  theme_void() +
  theme(plot.margin = unit(rep(2,4), "cm")) +
  coord_cartesian(clip = "off") +
  ylim(0,3)
print(p)



qualitative.shape.scale <- make_consistent_scale(values = c(24, 22, 25, 23),
                                                 vars = c("HIGH",  "LOW","MODERATE", "MODIFIER"),
                                                 type = "shape",
                                                 name = ,show.example = T)


