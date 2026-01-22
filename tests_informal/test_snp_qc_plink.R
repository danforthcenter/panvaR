set_out_dir("~/scratch/panvar_test/snpqc")

snp_qc_plink(genotype.path = "inst/extdata/Setaria_shattering_example_pruned.bed",
             out.name = "TESTING")

snp_qc_plink(genotype.path = "~/scratch/setaria_annotated_chrsnumeric.vcf.gz",
             out.name = "TESTING",
             sample.list.path = "~/scratch/Panvar_list.of.samples.with.phenotype_shattering.txt")
