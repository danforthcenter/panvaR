pheno <- read.delim("inst/extdata/Setaria_shattering_example_phenotype.tsv")

make_panvar_inputs(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
                   phenotype.table = pheno,
                   out.prefix = "GWAStest",
                   out.dir = "~/scratch/temp")

panvar_mvp_gwas(inputs.dir = "~/scratch/temp",
                in.prefix = "GWAStest",
                phenotype.table = pheno,
                npcs = 3,
                gwas.model = "GLM")


inputs.dir = "~/scratch/panvar_test/snpqc/"
in.prefix = "PanvarIN"

matched.files <- list.files(path = inputs.dir, pattern = paste0("^", in.prefix))

