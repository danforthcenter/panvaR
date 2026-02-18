# ------------------------------------------------------------------------\
# testing PCs --------
# ------------------------------------------------------------------------\

pheno <- read.delim("inst/extdata/Setaria_shattering_example_phenotype.tsv")

make_panvar_inputs(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
                   phenotype.path = "inst/extdata/Setaria_shattering_example_phenotype.tsv",
                   out.prefix = "INPUTSPHENO",
                   out.dir = "~/scratch/temp",
                   extra.plink.options = c("--max-maf", ".95"))

pcs <- rMVP::attach.big.matrix("~/scratch/temp/INPUTSPHENO.pc.desc")


# pcadapt package
library(pcadapt)

mat <- pcadapt::bed2matrix("~/scratch/temp/INPUTSPHENO_PlinkQC_maf0.05_missing0.1.bed")
mat <- pcadapt::bed2matrix(path.expand("~/scratch/temp/INPUTSPHENO_PlinkQC_maf0.05_missing0.1.bed"))

bedfile.path <- read.pcadapt(path.expand("~/scratch/temp/INPUTSPHENO_PlinkQC_maf0.05_missing0.1.bed"),
                             type = "bed")
x <- pcadapt(input = bedfile.path, K = 20) 
plot(x, option = "scree")

max.pc <- 20

plot.df <- data.frame(PC = 1:max.pc, prop.var = x$singular.values^2 / sum(x$singular.values ^ 2), eigenvalues = x$singular.values^2)
ggplot(aes(x = PC, y = prop.var), data = plot.df) +
  geom_point() +
  geom_line() +
  theme_bw(13) +
  labs(title = "Scree Plot",
       y = "Proportion of Variance Explained (approximate)")

test <- plot.df %>% 
  mutate(cumulative = cumsum(prop.var))

ggplot(aes(x = PC, y = cumulative), data = test) +
  geom_point() +
  geom_line() +
  theme_bw(13) +
  labs(title = "Scree Plot - Cumulative",
       y = "Cumulative Proportion of Variance Explained (approximate)")


plot.df <- data.frame(pc.num = c(1:10), pve = round(100 * cumsum(var_exp), 3))


# ------------------------------------------------------------------------\
# test function --------
# ------------------------------------------------------------------------\

plot_pc_scree("~/scratch/temp/INPUTSPHENO_PlinkQC_maf0.05_missing0.1.bed",
              10,
              cumulative = T)
