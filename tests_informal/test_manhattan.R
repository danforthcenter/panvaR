
pheno <- read.delim("inst/extdata/Setaria_shattering_example_phenotype.tsv")

make_panvar_inputs(genotype.path = "~/scratch/setaria_annotated_chrsnumeric.vcf.gz",
                   phenotype.path = "inst/extdata/Setaria_shattering_example_phenotype.tsv",
                   out.prefix = "SetShattering",
                   out.dir = "~/scratch/panvar_test/setaria_shatter_full")

set_plink_path("/home/cluebbert/bin/plink2")
set_panvar_prefix("SetShattering")
set_out_dir("~/scratch/panvar_test/setaria_shatter_full")


panvar_mvp_gwas(inputs.dir = options()$panvar_outdir,
                npcs = 3,
                gwas.model = "GLM",
                output.manhattan = T)

out_ld <- 
  panvaR::get_ld_in_window(tag.snp = "5-6857045",
                           window = 500,
                           in.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full/",
                           out.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full",
                           geno.bed = "SetShattering_PlinkQC_maf0.05_missing0.1")

x <- out_ld$table
x <- out_ld$table %>% 
  mutate(pos = get_bp_from_id(marker.ID))
range(x$pos)

list.files(options()$panvar_outdir)
gwas.df <- data.table::fread(file.path(options()$panvar_outdir, "SetShattering_GLM_GWASresults.csv"))

gwas.df[which.max(gwas.df$LOGPVAL), ]

make_panvar_manhattan(gwas.res = gwas.df,
                      qtl.df = NULL,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      window = 500,
                      sig.line = 6, orient = "V")
