# gwas.res <-
#   panvar_gwas(genotype_data = "inst/extdata/Setaria_shattering_example_pruned.bed",
#               phenotype_input = "inst/extdata/Setaria_shattering_example_phenotype.tsv",
#               maf = 0,
#               missing_rate = 1,
#               pc_min = 1)
# 
# test setting out dir
# works with either trailing slash or not, need full path though (doesn't expand "~/")
# set_out_dir("/home/cluebbert/scratch/")

out <-
panvar_func(vcf_file_path = "/home/cluebbert/scratch/setaria_annotated_chrsnumeric.vcf.gz",
            phenotype_data = "inst/extdata/Setaria_shattering_example_phenotype.tsv",
            tag_snps = "5:6857045")


out$plot


out <-
  panvar_func(vcf_file_path = "/home/cluebbert/scratch/setaria_annotated_chrsnumeric.vcf.gz",
              phenotype_data = "inst/extdata/Setaria_shattering_example_phenotype.tsv",
              tag_snps = "5:6857045",
              pc_rds = "inst/extdata/Setaria_shattering_example_PCobj_6pcs.rds")


x <- readRDS("inst/extdata/Setaria_shattering_example_PCobj_6pcs.rds")
y <- matrix(c(1,1))
# test custom pcs
cleaned_up_bed_file <- "inst/extdata/Setaria_shattering_example_pruned.bed"
rds_file_base <- base_name_func(cleaned_up_bed_file,include_dir = TRUE)
# the path to the rds file made by `rds_file_base`
rds_file_path <- paste0(c(rds_file_base,".rds"),collapse = "")
# Read the genotype using `snp_attach`
genotype_rds_data <- snp_attach("~/scratch/panvar_test/Setaria_shattering_example_pruned.rds")
genotype_info <- genotype_rds_data$fam
geno.order <- genotype_info$sample.ID


length(match(rownames(x), geno.order))
