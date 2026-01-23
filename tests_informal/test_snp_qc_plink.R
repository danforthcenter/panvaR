
# test output logic working
set_out_dir("~/scratch/panvar_test/temp2")
# if I don't supply, should give set_out_dir directory
snp_qc_plink(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
             out.prefix = "TESTING")
# if I do supply should over ride and give that one
snp_qc_plink(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
             out.prefix = "TESTING",
             out.dir = "~/scratch/temp")
options(panvar_outdir = NULL)
snp_qc_plink(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
             out.prefix = "TESTING")



# test plink executable logic working 
set_plink_path("/home/cluebbert/bin/plink2")
snp_qc_plink(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
             out.prefix = "TESTING",
             plink.path = "~/",
             out.dir = "~/scratch/temp")
set_plink_path(NULL)
snp_qc_plink(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
             out.prefix = "TESTING",
             plink.path = "/home/cluebbert/bin/plink2",
             out.dir = "~/scratch/temp")


# test with sample list
snp_qc_plink(genotype.path = "~/scratch/setaria_annotated_chrsnumeric.vcf.gz",
             out.name = "TESTING",
             sample.list.path = "~/scratch/Panvar_list.of.samples.with.phenotype_shattering.txt")


# test prefix logic
set_out_dir("~/scratch/panvar_test/temp3")
set_plink_path("/home/cluebbert/bin/plink2")
set_panvar_prefix("Awesome")
snp_qc_plink(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed")
set_panvar_prefix(NULL)
snp_qc_plink(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed")
