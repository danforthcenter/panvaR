# gwas.res <-
#   panvar_gwas(genotype_data = "inst/extdata/Setaria_shattering_example_pruned.bed",
#               phenotype_input = "inst/extdata/Setaria_shattering_example_phenotype.tsv",
#               maf = 0,
#               missing_rate = 1,
#               pc_min = 1)
# 
out <-
panvar_func(vcf_file_path = "/home/cluebbert/scratch/setaria_annotated_chrsnumeric.vcf.gz",
            phenotype_data = "inst/extdata/Setaria_shattering_example_phenotype.tsv",
            tag_snps = "5:6857045")


out$plot
