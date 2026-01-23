pheno <- read.delim("inst/extdata/Setaria_shattering_example_phenotype.tsv")

# right now there are non-numeric chromsomes in the example files I think, throws a warning from MVP
make_panvar_inputs(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
                   phenotype.table = pheno,
                   out.prefix = "GWAStest",
                   out.dir = "~/scratch/temp")
panvar_mvp_gwas(inputs.dir = "~/scratch/panvar_test/gwas/",
                in.prefix = "GWAStest",
                phenotype.table = pheno,
                npcs = 3,
                out.dir = "~/scratch/temp",
                gwas.model = "GLM")

# run full genome, unfixed chromosomes. does it break anything besides throwing a warning?
# outputs warning if plotting
make_panvar_inputs(genotype.path = "~/scratch/setaria_annotated_vcf.vcf.vcf",
                   phenotype.table = pheno,
                   out.prefix = "GWAStest",
                   out.dir = "~/scratch/panvar_test/gwas2/")
panvar_mvp_gwas(inputs.dir = "~/scratch/panvar_test/gwas2/",
                in.prefix = "GWAStest",
                phenotype.table = pheno,
                npcs = 3,
                output.manhattan = T,
                out.dir = "~/scratch/panvar_test/gwas2/",
                gwas.model = "GLM")

# test full genome
make_panvar_inputs(genotype.path = "~/scratch/setaria_annotated_chrsnumeric.vcf.gz",
                   phenotype.table = pheno,
                   out.prefix = "GWAStest",
                   out.dir = "~/scratch/panvar_test/gwas/")

panvar_mvp_gwas(inputs.dir = "~/scratch/panvar_test/gwas/",
                in.prefix = "GWAStest",
                phenotype.table = pheno,
                npcs = 3,
                out.dir = "~/scratch/panvar_test/gwas/",
                gwas.model = "GLM")

# test mlm
make_panvar_inputs(genotype.path = "~/scratch/setaria_annotated_chrsnumeric.vcf.gz",
                   phenotype.table = pheno,
                   calc.kinship = T,
                   out.prefix = "GWAStest",
                   out.dir = "~/scratch/panvar_test/gwas/")

panvar_mvp_gwas(inputs.dir = "~/scratch/panvar_test/gwas/",
                in.prefix = "GWAStest",
                phenotype.table = pheno,
                npcs = 3,
                out.dir = "~/scratch/panvar_test/gwas/",
                gwas.model = "MLM")


inputs.dir = "~/scratch/panvar_test/snpqc/"
in.prefix = "PanvarIN"

matched.files <- list.files(path = inputs.dir, pattern = paste0("^", in.prefix))

